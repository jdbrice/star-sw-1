#!/usr/bin/perl -w
# Submit pico analysis jobs (picoMatch + picoDilepton) for all three passes.
#
# Usage:
#   submit_ana.pl type=TYPE outdir=OUTDIR [option=NN] [maxjob=N] [run|submit|merge]
#
#   TYPE detection:
#     "data"           → real data  (reads runs2.txt + picolist/{RUN}.{N}.lis)
#     contains "pythia"→ Pythia MC  (reads picolist/{TYPE}.nsets.txt)
#     otherwise        → particle gun (reads picolist/{TYPE}.nsets.txt)
#   outdir = output label under DataDisk/  (e.g. 202606 or ele.pt1.5.vz0)
#   option = analysis bits  default 11 (bit0=picoMatch, bit1=picoDilepton)
#   maxjob = condor job cap  default 10000
#   run    = execute locally instead of condor
#   submit = submit to condor
#   merge  = run merge_ana.pl after submit/run
#
# Examples:
#   submit_ana.pl type=ele.pt1.5.vz0   outdir=ele.pt1.5.vz0   submit
#   submit_ana.pl type=pythia.JPsi.vz0 outdir=pythia.JPsi.vz0 submit
#   submit_ana.pl type=data             outdir=202606           submit merge

use strict;
use warnings;

# === CONFIGURATION ===
my $WORKDIR   = "/star/u/akio/fcstrk11/star-sw-fwd";
my $DATADIR   = "/gpfs01/star/pwg_tasks/FwdCalib/akio";
my $CONTAINER = "/cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif";
my $CONDORDIR = "$WORKDIR/condor";
my $SGL       = "singularity exec -e -B /direct -B /star -B /afs -B /gpfs -B /sdcc/lustre02 $CONTAINER";
# =====================

my %args;
my ($do_submit, $do_run, $do_merge) = (0, 0, 0);
for my $arg (@ARGV) {
    if    ($arg eq 'submit') { $do_submit = 1; }
    elsif ($arg eq 'run')    { $do_run    = 1; }
    elsif ($arg eq 'merge')  { $do_merge  = 1; }
    elsif ($arg =~ /^([\w.]+)=(.+)$/) { $args{$1} = $2; }
    else  { die "Unknown argument: $arg\n"; }
}

my $type   = $args{type}   or die "type= required\n";
my $outdir = $args{outdir} or die "outdir= required\n";
my $option = $args{option} // 11;
my $maxjob = $args{maxjob} // 10000;

my $exe     = "$WORKDIR/script/runpico_uni";
my $fullout = "$DATADIR/$outdir";
my $logdir  = "$fullout/log";
mkdir_p($CONDORDIR); mkdir_p($fullout);
mkdir_p("$fullout/hist_match"); mkdir_p("$fullout/hist_dilep"); mkdir_p($logdir);

if ($type eq 'data') {
    submit_data($type, $option, $maxjob, $fullout, $logdir, $exe, $do_submit, $do_run);
} else {
    submit_sim($type, $option, $maxjob, $fullout, $logdir, $exe, $do_submit, $do_run);
}

if ($do_merge) {
    system("$WORKDIR/script/merge_ana.pl outdir=$outdir");
}

# --------------------------------------------------------------------------
sub submit_sim {
    my ($type, $option, $maxjob, $fullout, $logdir, $exe, $submit, $run_local) = @_;

    my $nsets_file = "$WORKDIR/picolist/${type}.nsets.txt";
    die "picolist/${type}.nsets.txt not found.\nRun: find_pico.pl type=$type  first.\n"
        unless -e $nsets_file;
    open(my $fh, '<', $nsets_file) or die;
    my $nsets = <$fh> + 0;  close $fh;
    die "nsets=0 — nothing to do\n" unless $nsets > 0;

    my $condor = "$CONDORDIR/submit_ana_${type}.txt";
    open(my $cfh, '>', $condor) or die;
    print $cfh "Executable = /bin/env\nUniverse = vanilla\n";
    print $cfh "notification = never\ngetenv = True\n\n";

    my $njob = 0;
    for my $iset (0 .. $nsets - 1) {
        last if $njob >= $maxjob;
        my $run = 1;   # readpico3 dispatch: run<10^7 → picolist_pythia/{iset}.list
        my $log = "$logdir/${type}.${iset}.log";
        # runpico_uni TYPE RUN ISET OPTION WORKDIR OUTDIR
        print $cfh "Arguments = \"$SGL $exe $type $run $iset $option $WORKDIR $fullout\"\n";
        print $cfh "Log = $log\nOutput = $log\nError = $log\nQueue\n\n";
        system("$exe $type $run $iset $option $WORKDIR $fullout") if $run_local;
        $njob++;
    }
    close $cfh;
    printf "%d jobs → %s\n", $njob, $condor;
    submit_or_not($condor, $submit, $run_local);
}

# --------------------------------------------------------------------------
sub submit_data {
    my ($type, $option, $maxjob, $fullout, $logdir, $exe, $submit, $run_local) = @_;

    my $runs2 = "$WORKDIR/runs2.txt";
    die "runs2.txt not found.\nRun: find_pico.pl type=data  first.\n" unless -e $runs2;

    my $condor = "$CONDORDIR/submit_ana_data.txt";
    open(my $cfh, '>', $condor) or die;
    print $cfh "Executable = /bin/env\nUniverse = vanilla\n";
    print $cfh "notification = never\ngetenv = True\n\n";

    my $njob = 0;
    my $done = 0;
    open(my $rfh, '<', $runs2) or die;
    while (!$done && (my $line = <$rfh>)) {
        chomp $line;
        my ($run, $nsets) = split /\s+/, $line;
        next unless $run =~ /^\d+$/ && $nsets >= 1;
        for my $iset (0 .. $nsets - 1) {
            if ($njob >= $maxjob) { $done = 1; last; }
            my $log = "$logdir/$run.$iset.log";
            # runpico_uni TYPE RUN ISET OPTION WORKDIR OUTDIR
            print $cfh "Arguments = \"$SGL $exe $type $run $iset $option $WORKDIR $fullout\"\n";
            print $cfh "Log = $log\nOutput = $log\nError = $log\nQueue\n\n";
            system("$exe $type $run $iset $option $WORKDIR $fullout") if $run_local;
            $njob++;
        }
    }
    close $rfh; close $cfh;
    printf "%d jobs → %s\n", $njob, $condor;
    submit_or_not($condor, $submit, $run_local);
}

# --------------------------------------------------------------------------
sub submit_or_not {
    my ($condor, $submit, $run_local) = @_;
    if ($submit) {
        print "Submitting $condor\n";
        system("condor_submit $condor");
        system("$WORKDIR/running200.pl");
    } elsif (!$run_local) {
        print "  (add 'submit' or 'run' to execute)\n";
    }
}

sub mkdir_p { system("/bin/mkdir -p $_[0]") unless -d $_[0]; }
