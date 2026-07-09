#!/usr/bin/perl -w
# Merge analysis histogram files with hadd.
#
# Usage:
#   merge_ana.pl outdir=OUTDIR [types=match,dilep]
#
#   outdir = directory label under $DATADIR (e.g. 202606)
#   types  = comma-separated list of histogram types (default: match,dilep)
#
# Creates:
#   $DATADIR/$OUTDIR/hist_match/all.match.root
#   $DATADIR/$OUTDIR/hist_dilep/all.dilep.root
#
# Example:
#   merge_ana.pl outdir=202606
#   merge_ana.pl outdir=202606 types=match

use strict;
use warnings;

# === CONFIGURATION ===
my $DATADIR = "/gpfs01/star/pwg_tasks/FwdCalib/akio";
# =====================

my %args;
for my $arg (@ARGV) {
    if ($arg =~ /^(\w+)=(.+)$/) { $args{$1} = $2; }
    else { die "Unknown argument: $arg\n"; }
}

my $outdir_label = $args{outdir} or die "outdir= required\n";
my $types_str    = $args{types}  // "match,dilep";
my @types        = split(/,/, $types_str);

my $outdir = "$DATADIR/$outdir_label";
die "Output directory not found: $outdir\n" unless -d $outdir;

for my $hist (@types) {
    my $histdir  = "$outdir/hist_$hist";
    my $allfile  = "$histdir/all.$hist.root";
    my $pattern  = "$hist.root";

    unless (-d $histdir) {
        warn "WARNING: $histdir not found, skipping $hist\n";
        next;
    }

    opendir(my $dh, $histdir) or die "Cannot open $histdir: $!\n";
    my @files = grep { /\.\Q$hist\E\.root$/ && $_ !~ /^all\./ } readdir($dh);
    closedir $dh;

    if (@files == 0) {
        warn "WARNING: no *.$hist.root files found in $histdir\n";
        next;
    }

    printf "Merging %d files → %s\n", scalar(@files), $allfile;
    unlink $allfile if -e $allfile;

    my $inputs = join(' ', map { "$histdir/$_" } sort @files);
    my $cmd = "hadd -k $allfile $inputs";
    print "$cmd\n";
    my $rc = system($cmd);
    if ($rc != 0) {
        warn "WARNING: hadd returned non-zero status for $hist\n";
    }
}
