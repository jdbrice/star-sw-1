#!/usr/bin/perl -w
# Submit simulation jobs (particle gun or Pythia) to condor.
#
# Usage:
#   submit_sim.pl type=TYPE [nev=N] [r1=N] [r2=N] [submit]
#
#   TYPE detection:
#     contains "pythia"  → Pythia MC    e.g. pythia.JPsi.vz0
#     otherwise          → particle gun e.g. ele.pt1.5.vz0  pos.e5.vz10
#   nev   events per job  (default: 1000 particle, 200 Pythia)
#   r1/r2 run number range (default: 1..100)
#   submit  actually submit to condor (default: just create job file)
#
# Examples:
#   submit_sim.pl type=ele.pt1.5.vz0     nev=1000 r1=1 r2=100 submit
#   submit_sim.pl type=ele.e5.vz0        nev=1000 r1=1 r2=50  submit
#   submit_sim.pl type=pythia.JPsi.vz0   nev=200  r1=1 r2=100 submit
#   submit_sim.pl type=pythia.dy.vz0     nev=200  r1=1 r2=100 submit
#
#   NOTE: as of 2026-07-13, runPythia.C's vz argument is the GAUSSIAN
#   SIGMA of the vertex (vertex mean fixed at 0,0,0 with sigmaX=sigmaY=
#   0.1cm), not the vertex center -- e.g. vz40 means a realistic
#   sigma_z=40cm spread around z=0, not a vertex offset to z=40. Used
#   for realistic-vertex J/psi and DYBG (Drell-Yan background) MC
#   samples (see todo.txt #10):
#
#   submit_sim.pl type=pythia.JPsi.vz40  nev=200  r1=1 r2=200  submit
#   submit_sim.pl type=pythia.dybg.vz40  nev=500  r1=1 r2=2000 submit

use strict;
use warnings;

# === CONFIGURATION ===
my $WORKDIR   = "/star/u/akio/fcstrk11/star-sw-fwd";
my $DATADIR   = "/gpfs01/star/pwg_tasks/FwdCalib/akio";
my $CONTAINER = "/cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif";
my $CONDORDIR = "$WORKDIR/condor";
my $SGL       = "singularity exec -e -B /direct -B /star -B /afs -B /gpfs -B /sdcc/lustre02 $CONTAINER";
# MC production (submitted here) is bulk/backlog work that should never
# compete with interactive real-data analysis batches (submit_ana.pl,
# submit_check.pl -- left at HTCondor's default priority 0). Negative
# priority = scheduled after default-priority jobs whenever both have idle
# jobs waiting for this user's slots. Doesn't preempt already-running jobs,
# but gets real-data batches to the front of the queue as slots turn over.
my $MC_PRIORITY = -10;
# =====================

my %args;
my $do_submit = 0;
for my $arg (@ARGV) {
    if    ($arg eq 'submit') { $do_submit = 1; }
    elsif ($arg =~ /^([\w.]+)=(.+)$/) { $args{$1} = $2; }
    else  { die "Unknown argument: $arg\n"; }
}

my $type = $args{type} or die "type= required\n";

if ($type =~ /pythia/i) {
    run_pythia($type, \%args, $do_submit);
} elsif ($type eq 'data') {
    die "type=data is not valid for submit_sim.pl — use submit_check.pl for real data\n";
} else {
    run_particle($type, \%args, $do_submit);
}

# --------------------------------------------------------------------------
sub run_particle {
    my ($type, $args, $submit) = @_;

    my ($pid, $spec, $vz);
    if    ($type =~ /^(\w+)\.(pt[\d.]+)\.vz(-?\d+)$/) { ($pid,$spec,$vz) = ($1,$2,$3); }
    elsif ($type =~ /^(\w+)\.(e[\d.]+)\.vz(-?\d+)$/)  { ($pid,$spec,$vz) = ($1,$2,$3); }
    else  { die "Cannot parse particle TYPE='$type'.\nExpected: {pid}.pt{PT}.vz{VZ} or {pid}.e{E}.vz{VZ}\n"; }

    my $nev = $args->{nev} // 1000;
    my $r1  = $args->{r1}  // 1;
    my $r2  = $args->{r2}  // 100;
    printf "particle  type=%s  nev=%d  runs=%d..%d\n", $type, $nev, $r1, $r2;

    my $outdir  = "$DATADIR/particle/$type";
    my $logdir  = "$outdir/log";
    my $exe     = "$WORKDIR/script/runsim";
    my $condor  = "$CONDORDIR/submit_sim_${type}.txt";

    mkdir_p($CONDORDIR); mkdir_p($outdir); mkdir_p($logdir);

    open(my $fh, '>', $condor) or die "Cannot write $condor: $!\n";
    print $fh "Executable   = /bin/env\nUniverse     = vanilla\n";
    print $fh "notification = never\ngetenv       = True\npriority     = $MC_PRIORITY\n\n";

    my $njob = 0;
    for my $r ($r1 .. $r2) {
        my $log = "$logdir/${type}.run${r}.log";
        print $fh "Arguments = \"$SGL $exe $nev $r $pid $spec $vz $outdir $WORKDIR\"\n";
        print $fh "Log = $log\nOutput = $log\nError = $log\nQueue\n\n";
        $njob++;
    }
    close $fh;
    printf "%d jobs → %s\n", $njob, $condor;
    submit_or_not($condor, $submit);
}

# --------------------------------------------------------------------------
sub run_pythia {
    my ($type, $args, $submit) = @_;

    my ($proc, $vz);
    if ($type =~ /^pythia\.(\w+)\.vz(-?\d+)$/i) { ($proc,$vz) = ($1,$2); }
    else { die "Cannot parse Pythia TYPE='$type'.\nExpected: pythia.{PROC}.vz{VZ}\n"; }

    my $nev = $args->{nev} // 200;
    my $r1  = $args->{r1}  // 1;
    my $r2  = $args->{r2}  // 100;
    printf "pythia  type=%s  nev=%d  runs=%d..%d\n", $type, $nev, $r1, $r2;

    my $outdir  = "$DATADIR/pythia/$type";
    my $logdir  = "$outdir/log";
    my $exe     = "$WORKDIR/runpythia";
    my $condor  = "$CONDORDIR/submit_sim_${type}.txt";

    mkdir_p($CONDORDIR); mkdir_p($outdir); mkdir_p($logdir);

    open(my $fh, '>', $condor) or die "Cannot write $condor: $!\n";
    print $fh "Executable   = /bin/env\nUniverse     = vanilla\n";
    print $fh "notification = never\ngetenv       = True\npriority     = $MC_PRIORITY\n\n";

    my $njob = 0;
    for my $r ($r1 .. $r2) {
        my $log = "$logdir/${type}.run${r}.log";
        # runpythia: NEV RUN PARTICLE VZ OUTDIR GEANT BFC
        print $fh "Arguments = \"$SGL $exe $nev $r $proc $vz $outdir 1 1\"\n";
        print $fh "Log = $log\nOutput = $log\nError = $log\nQueue\n\n";
        $njob++;
    }
    close $fh;
    printf "%d jobs → %s\n", $njob, $condor;
    submit_or_not($condor, $submit);
}

# --------------------------------------------------------------------------
sub submit_or_not {
    my ($condor, $submit) = @_;
    if ($submit) { print "Submitting $condor\n"; system("condor_submit $condor"); }
    else         { print "  (add 'submit' to actually submit)\n"; }
}

sub mkdir_p { system("/bin/mkdir -p $_[0]") unless -d $_[0]; }
