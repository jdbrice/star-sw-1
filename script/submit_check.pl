#!/usr/bin/perl -w
# Submit real-data afterburner jobs (MuDst -> picoDst), skipping runs whose
# picoDst output already exists (resumable version of the afterburner submit).
#
# Usage:
#   submit_check.pl [submit|submitmerge] [align] [maxjobs]
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

my $opt      = "none";
my $doAlign  = 0;
if (@ARGV >= 1) { $opt = $ARGV[0]; }
if (@ARGV >= 2 && $ARGV[1] eq "align") { $doAlign = 1; }
my $maxjobArg = ($doAlign ? $maxAlignJobs : $maxjob);
if (@ARGV >= 3) { $maxjobArg = $ARGV[2] + 0; }
print "Option = $opt, align = $doAlign, maxjobs this invocation = $maxjobArg\n";
if (@ARGV < 1) {
    print "No option given: dry run only. Add 'submit' to actually queue jobs.\n";
}

my $outdir = "$DATADIR/pico";
my $logdir = "$DATADIR/pico/log";
my $exe    = "$WORKDIR/script/runab";
my $sgl    = "$SGL $exe";

print "WORKDIR=$WORKDIR\n";

system("/bin/mkdir -p $CONDORDIR") unless -d $CONDORDIR;
system("/bin/mkdir -p $logdir")    unless -d $logdir;

my $condor = "$CONDORDIR/submit.txt";
print("Creating $condor\n");
unlink $condor if -e $condor;

open(my $OUT, '>', $condor) or die "Cannot write $condor: $!\n";
print $OUT "Executable   = /bin/env\n";
print $OUT "Universe     = vanilla\n";
print $OUT "notification = never\n";
print $OUT "getenv       = True\n";
print $OUT "\n";

my $njob  = 0;
my $nskip = 0;
open(my $IN, '<', "$WORKDIR/prod_mudst.txt") or die "Cannot read $WORKDIR/prod_mudst.txt: $!\n";
while (my $file = <$IN>) {
    chomp($file);
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
    print $OUT "Arguments = \"$sgl $inputdir $base $outdir $alignArg\"\n";
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
