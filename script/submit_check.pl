#!/usr/bin/perl -w
# Submit real-data afterburner jobs (MuDst -> picoDst), skipping runs whose
# picoDst output already exists (resumable version of the afterburner submit).
#
# Usage:
#   submit_check.pl [submit|submitmerge] [align] [maxjobs] [key=value ...]
#
#   (no option)  dry run — just reports how many jobs would be queued/skipped
#   submit       submit to condor
#   submitmerge  submit to condor, same as submit (merge step not yet wired up)
#   align        also run StFwdAlignmentMaker (unbiased hit-removed residuals,
#                see proposal_alignment_path.txt) for each file in this batch,
#                and copy the resulting FwdAlignment_*.root ntuples back to
#                $outdir/alignment/. Off (0) if omitted -- alignment costs
#                real extra walltime per job (~20% in interactive testing), so
#                this is opt-in per batch, not the routine picoDst-production
#                default. When given, the per-batch job cap defaults to a much
#                smaller $maxAlignJobs instead of $maxjob (see below), since
#                this is meant for a focused stats campaign, not the full
#                3800+-file production list -- pass an explicit maxjobs as the
#                3rd argument to override that cap either way.
#   maxjobs      optional explicit cap on jobs queued this invocation,
#                overriding either default below.
#
# Dataset selection (key=value, any order, after the positional args above).
# Defaults reproduce the old hardcoded behaviour exactly, so existing
# invocations are unchanged:
#
#   list=FILE      input MuDst list, one path per line.
#                  Relative names resolve against $WORKDIR.
#                  default: prod_mudst.txt  (-> the standard FwdStreamTest
#                  production_pp500_2022_DEV_st_fwd_test_mudst.txt)
#   out=DIR        output directory for picoDst + residual/ + alignment/ + log/.
#                  Relative names resolve against $DATADIR.
#                  default: $DATADIR/pico
#   magfield=N     magnetic field mode (see runab / fwd_afterburner_db.C):
#                    1  from the MuDst's own magneticField() stamp (default)
#                    0  forced zero
#                   -1  legacy, no field at all -- what every afterburner
#                       production before 2026-09-09 actually did
#                  Nothing needs setting per dataset: mode 1 gets field-on and
#                  field-off runs both right, straight out of the file.
#                  default: 1
#   gapfix=1       applyFstGapFix: undo the FST outer-sensor 1 deg gap sign at
#                  load time. Needed for every MuDst written before
#                  StFstHitMaker is fixed, i.e. all current productions.
#                  default: 0
#   fstmirror=1    applyFstMirror: per-wedge FST mirror. TEMPORARY -- the real
#                  fix starts in StarVMC/Geometry/FstmGeo/FstmGeo.xml and has to
#                  include the per-wedge z. default: 0
#   fttmirror=1    applyFttXYMirror: exchange which FTT strip orientation
#                  measures X and which measures Y. TEMPORARY, and note an
#                  equivalent 2026-07 swap LOWERED purity 14.9% -> 7.0%.
#                  default: 0
#   ftttimecut=N   StFttClusterMaker time-cut mode: 2 = calibrated time
#                  (default, production), 1 = AcceptAll. Use 1 for sparse
#                  data like zf2022/zf2024: per-VMM time calibration needs
#                  >200 distinct dbcid values per VMM and essentially never
#                  warms up there, so mode 2 keeps almost no FTT hits.
#                  default: 2
#   submitdir=DIR  where to submit from (sets condor Initialdir). The job script
#                  resolves StRoot/, fGeom.root and the macro from $cwd, so this
#                  must be the working copy. default: $WORKDIR
#
# Named datasets, just shorthand for the matching list=:
#   zf2022         list=zeroFieldAlignment_2022_DEV_mudst.txt
#   zf2024         list=zeroFieldAlignment_2024_DEV_mudst.txt
# (out= is still required with these -- the output directory is deliberately not
#  guessed, so a campaign can never silently land on top of another one.)
#
# NOTE: the "skip if picoDst exists" check below only looks at picoDst -- if
# you run plain `submit` first and later want `align` for the SAME files,
# they will be skipped here (picoDst already exists) and won't get an
# alignment ntuable produced. Point at a subset of prod_mudst.txt (or clear
# those picoDst outputs) if you need to add alignment to already-processed
# files.
#
# Reads:
#   $WORKDIR/prod_mudst.txt   — one MuDst path per line
# Writes:
#   $DATADIR/pico/*.picoDst.root
#   $DATADIR/pico/residual/*.FwdDetResidual_*.root
#   $DATADIR/pico/alignment/*.FwdAlignment_*.root   (only with align)
#   condor/submit.txt
#
# Examples:
#   submit_check.pl submit                  # normal picoDst+residual production, all files
#   submit_check.pl submit align            # + alignment ntuple, capped at $maxAlignJobs files
#   submit_check.pl submit align 100        # + alignment ntuple, capped at 100 files
#
#   # zero-field 2022 -> its own output dir (dry run first, then submit)
#   submit_check.pl zf2022 out=pico_zf2022_20260909
#   submit_check.pl submit zf2022 out=pico_zf2022_20260909
#
#   # equivalent, spelled out
#   submit_check.pl submit list=zeroFieldAlignment_2022_DEV_mudst.txt \
#                          out=pico_zf2022_20260909
#
# Arguments are matched by CONTENT, not position, so order does not matter and
# no placeholder is ever needed to reach a later one.

use strict;
use warnings;
use File::Basename;

# === CONFIGURATION ===
my $WORKDIR      = "/star/u/akio/fcstrk11/star-sw-fwd";
my $DATADIR      = "/gpfs01/star/pwg_tasks/FwdCalib/akio";
my $CONTAINER    = "/cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif";
my $CONDORDIR    = "$WORKDIR/condor";
my $SGL          = "singularity exec -e -B /direct -B /star -B /afs -B /gpfs -B /sdcc/lustre02 $CONTAINER";
my $maxjob       = 10000; # picoDst-only production: no practical cap (list has ~3846 files)
my $maxAlignJobs = 30;    # alignment campaigns: conservative default -- see usage note above
# =====================

# --- argument parsing -------------------------------------------------------
# Matched by content, not position: an explicit maxjobs no longer requires
# passing "align"/"none" placeholders just to reach the third slot.
my $opt       = "none";
my $doAlign   = 0;
my $maxjobArg = undef;      # resolved after we know $doAlign
my $list      = "prod_mudst.txt";
my $out       = "pico";
my $magfield  = 1;
my $gapfix    = 0;
my $fstmirror = 0;
my $fttmirror = 0;
my $ftttimecut = 2;
my $fttwindow  = "-40,100";   # dbcid ticks relative to the per-VMM anchor
my $diag       = "-1,0,0";    # FST-blind diagnostic: TYPE,NOADD,MIX (see runab)
my $submitdir = $WORKDIR;
my $outGiven  = 0;
my $listGiven = 0;

foreach my $a (@ARGV) {
    if    ($a eq "submit" || $a eq "submitmerge") { $opt = $a; }
    elsif ($a eq "none")                          { $opt = "none"; }
    elsif ($a eq "align")                         { $doAlign = 1; }
    elsif ($a =~ /^\d+$/)                         { $maxjobArg = $a + 0; }
    elsif ($a eq "zf2022") { $list = "zeroFieldAlignment_2022_DEV_mudst.txt"; $listGiven = 1; }
    elsif ($a eq "zf2024") { $list = "zeroFieldAlignment_2024_DEV_mudst.txt"; $listGiven = 1; }
    elsif ($a =~ /^list=(.+)$/)          { $list = $1; $listGiven = 1; }
    elsif ($a =~ /^out=(.+)$/)           { $out = $1; $outGiven = 1; }
    elsif ($a =~ /^magfield=(-?\d+)$/)   { $magfield = $1 + 0; }
    elsif ($a =~ /^gapfix=(\d+)$/)       { $gapfix = $1 + 0; }
    elsif ($a =~ /^fstmirror=(\d+)$/)    { $fstmirror = $1 + 0; }
    elsif ($a =~ /^fttmirror=(\d+)$/)    { $fttmirror = $1 + 0; }
    elsif ($a =~ /^ftttimecut=(\d+)$/)   { $ftttimecut = $1 + 0; }
    elsif ($a =~ /^fttwindow=(-?\d+,-?\d+)$/) { $fttwindow = $1; }
    elsif ($a =~ /^diag=(-?\d+,[01],[01])$/)  { $diag = $1; }
    elsif ($a =~ /^submitdir=(.+)$/)     { $submitdir = $1; }
    else { die "Unrecognized argument '$a'. See the usage block at the top of $0.\n"; }
}
$maxjobArg = ($doAlign ? $maxAlignJobs : $maxjob) unless defined $maxjobArg;

# Resolve relative names; absolute paths pass through untouched.
my $listfile = ($list =~ m{^/}) ? $list : "$WORKDIR/$list";
my $outdir   = ($out  =~ m{^/}) ? $out  : "$DATADIR/$out";
my $logdir   = "$outdir/log";

# Refuse to write a non-default dataset into the default output directory: that
# is how one campaign silently lands on top of another, and the "skip if pico
# exists" check below would then read the WRONG files as already done.
if ($listGiven && !$outGiven) {
    die "a non-default list= was given without out=: refusing to write it into the\n"
      . "default $outdir. Pass an explicit out=DIR for this campaign.\n";
}
die "magfield must be 1, 0 or -1 (got $magfield)\n" unless grep { $_ == $magfield } (1, 0, -1);
die "ftttimecut must be 1 or 2 (got $ftttimecut)\n" unless $ftttimecut == 1 || $ftttimecut == 2;
{   my ($lo, $hi) = split /,/, $fttwindow;
    die "fttwindow must be LO,HI with LO < HI (got $fttwindow)\n" unless defined $hi && $lo < $hi;
    # Per-run anchors live in ./fttDataWindow and are shipped by runab; without
    # them a tight window is meaningless (the DB entry is a single global one).
    if ($hi - $lo < 40 && ! -d "$WORKDIR/fttDataWindow") {
        die "fttwindow=$fttwindow is tight but $WORKDIR/fttDataWindow does not exist\n";
    }
}
die "Input list not found: $listfile\n" unless -e $listfile;
die "Output dir does not exist: $outdir\n(create it first -- not auto-created, so a\n"
  . "typo cannot silently start a fresh campaign in the wrong place.)\n" unless -d $outdir;

my $exe    = "$WORKDIR/script/runab";
my $sgl    = "$SGL $exe";

print "Option = $opt, align = $doAlign, maxjobs this invocation = $maxjobArg\n";
print "list       = $listfile\n";
print "outdir     = $outdir\n";
my %mfdesc = (1 => "from MuDst magneticField()", 0 => "forced zero", -1 => "LEGACY, no field at all");
print "magfield   = $magfield   ($mfdesc{$magfield})\n";
print "geofix     = gapfix=$gapfix fstmirror=$fstmirror fttmirror=$fttmirror\n";
print "ftttimecut = $ftttimecut   (" . ($ftttimecut == 1 ? "AcceptAll" : "calibrated time") . ")\n";
print "fttwindow  = [$fttwindow]   (dbcid ticks; ignored when ftttimecut=1)\n";
{ my ($dt,$dn,$dm) = split /,/, $diag;
  print "diag       = type $dt (" . ($dt<0 ? "all types" : "one track type") . "), noAdd $dn, mix $dm\n"; }
print "submitdir  = $submitdir\n";
if (!grep { $_ eq "submit" || $_ eq "submitmerge" } @ARGV) {
    print "No 'submit' given: dry run only. Add 'submit' to actually queue jobs.\n";
}

print "WORKDIR=$WORKDIR\n";

system("/bin/mkdir -p $CONDORDIR") unless -d $CONDORDIR;
system("/bin/mkdir -p $logdir")    unless -d $logdir;

# Name the submit file after the campaign so two datasets can be prepared
# without one overwriting the other's queue.
my $tag = basename($outdir);
my $condor = "$CONDORDIR/submit_$tag.txt";
print("Creating $condor\n");
unlink $condor if -e $condor;

open(my $OUT, '>', $condor) or die "Cannot write $condor: $!\n";
print $OUT "Executable   = /bin/env\n";
print $OUT "Universe     = vanilla\n";
print $OUT "notification = never\n";
print $OUT "getenv       = True\n";
# runab resolves StRoot/, fGeom.root and the macro from $cwd. Condor defaults
# $cwd to the submit directory, so this was previously implicit and correct only
# while submitting from $WORKDIR. Stating it makes the job independent of where
# condor_submit happens to be run.
print $OUT "Initialdir   = $submitdir\n";
print $OUT "\n";

my $njob  = 0;
my $nskip = 0;
open(my $IN, '<', $listfile) or die "Cannot read $listfile: $!\n";
while (my $file = <$IN>) {
    chomp($file);
    next if $file =~ /^\s*$/;
    my $base      = basename($file);
    my $inputdir  = dirname($file);
    my $picobase  = $base;
    $picobase     =~ s/\.MuDst\.root$/.picoDst.root/;
    my $picofile  = "$outdir/$picobase";

    if (-e $picofile) {
        print "SKIP (pico exists): $picobase\n";
        $nskip++;
        next;
    }

    my $log = "$logdir/$base";
    $log =~ s/\.MuDst\.root$/.log/;
    my $alignArg = $doAlign ? "1" : "0";
    print $OUT "Arguments = \"$sgl $inputdir $base $outdir $alignArg $magfield $gapfix,$fstmirror,$fttmirror $ftttimecut $fttwindow $diag\"\n";
    print $OUT "Log    = $log\n";
    print $OUT "Output = $log\n";
    print $OUT "Error  = $log\n";
    print $OUT "Queue\n\n";
    $njob++;
    last if $njob >= $maxjobArg;
}
close($IN);
close($OUT);
print "$nskip jobs skipped (pico already exists)\n";
print "$njob jobs queued\n";

if ($opt eq "submit" || $opt eq "submitmerge") {
    print("Submitting ${condor}\n");
    system("condor_submit ${condor}");
    system("$WORKDIR/running200.pl");
}
#if($opt eq "merge" || $opt eq "submitmerge"){
#    system("merge.pl\n");
#}
