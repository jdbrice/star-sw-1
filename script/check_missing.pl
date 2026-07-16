#!/usr/bin/perl
# Check for missing picoDst.root output from submit_sim.pl MC production,
# and optionally resubmit just the missing runs.
#
# Usage:
#   check_missing.pl type=TYPE r1=N r2=N [nev=N] [resubmit]
#
#   TYPE detection (same as find_pico.pl / submit_sim.pl):
#     contains "pythia" → Pythia MC   (e.g. pythia.dybg.vz40)
#     otherwise         → particle gun (e.g. ele.pt1.5.vz0)
#   r1,r2    = run range to check (inclusive)
#   nev      = events/run, required only with 'resubmit' (passed to submit_sim.pl)
#   resubmit = submit_sim.pl r1=r2=run for each missing run
#
# Examples:
#   check_missing.pl type=pythia.JPsi.vz40 r1=1 r2=200
#   check_missing.pl type=pythia.dybg.vz40 r1=1 r2=20000 nev=500 resubmit

use strict;
use warnings;

# === CONFIGURATION ===
my $WORKDIR = "/star/u/akio/fcstrk11/star-sw-fwd";
my $DATADIR = "/gpfs01/star/pwg_tasks/FwdCalib/akio";
# =====================

my %args;
my $do_resubmit = 0;
for my $arg (@ARGV) {
    if    ($arg eq 'resubmit') { $do_resubmit = 1; }
    elsif ($arg =~ /^([\w.]+)=(.+)$/) { $args{$1} = $2; }
    else  { die "Unknown argument: $arg\n"; }
}

my $type = $args{type} or die "type= required\n";
my $r1   = $args{r1}   or die "r1= required\n";
my $r2   = $args{r2}   or die "r2= required\n";
my $nev  = $args{nev};
die "nev= required with 'resubmit'\n" if $do_resubmit && !defined $nev;

my $subdir = ($type =~ /pythia/i) ? "pythia/$type" : "particle/$type";
my $srcdir = "$DATADIR/$subdir";

my @missing;
for my $r ($r1 .. $r2) {
    my $f = "$srcdir/${type}.run${r}.picoDst.root";
    push @missing, $r unless -e $f;
}

my $nmiss = scalar @missing;
printf "Checked runs %d..%d in %s\n", $r1, $r2, $srcdir;
printf "Missing: %d / %d\n", $nmiss, $r2 - $r1 + 1;
print "@missing\n" if $nmiss;

exit 0 unless $nmiss;

if ($do_resubmit) {
    for my $r (@missing) {
        printf "Resubmitting run %d ...\n", $r;
        system("$WORKDIR/script/submit_sim.pl", "type=$type", "nev=$nev", "r1=$r", "r2=$r", "submit");
    }
} else {
    print "  (add 'resubmit nev=N' to resubmit the missing runs)\n";
}
