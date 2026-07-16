#!/usr/bin/perl
# Catalog picoDst files into list shards for readpico_uni / readpico3.C.
#
# Usage:
#   find_pico.pl type=TYPE [maxfile=N]
#
#   TYPE detection:
#     "data"          → real data (searches DataDisk/pico/)
#     contains "pythia" → Pythia MC  (e.g. pythia.JPsi.vz0)
#     otherwise        → particle gun (e.g. ele.pt1.5.vz0)
#
#   maxfile = max pico files per list shard (default 100)
#
# Output (all under picolist/ in $WORKDIR):
#   Particle / Pythia:
#     picolist/{TYPE}.{N}.list      — list shards (N = 0, 1, 2, ...)
#     picolist/{TYPE}.nsets.txt     — number of shards produced
#   Data:
#     picolist/{RUN}.{N}.lis        — one shard per run per set
#     runs2.txt                     — "RUN NSETS" pairs for submit_ana.pl
#
# Examples:
#   find_pico.pl type=ele.pt1.5.vz0
#   find_pico.pl type=pythia.JPsi.vz0
#   find_pico.pl type=data

use strict;
use warnings;

# === CONFIGURATION ===
my $WORKDIR = "/star/u/akio/fcstrk11/star-sw-fwd";
my $DATADIR = "/gpfs01/star/pwg_tasks/FwdCalib/akio";
# =====================

my %args;
for my $arg (@ARGV) {
    if ($arg =~ /^(\w[\w.]*?)=(.+)$/) { $args{$1} = $2; }
    else { die "Unknown argument: $arg\n"; }
}

my $type    = $args{type}    or die "type= required\n";
my $maxfile = $args{maxfile} // 100;

my $picolistdir = "$WORKDIR/picolist";
system("/bin/mkdir -p $picolistdir") unless -d $picolistdir;

if ($type eq 'data') {
    find_data($maxfile, $picolistdir);
} elsif ($type =~ /pythia/i) {
    find_sim($type, "pythia/$type", $maxfile, $picolistdir);
} else {
    find_sim($type, "particle/$type", $maxfile, $picolistdir);
}

# --------------------------------------------------------------------------
sub find_sim {
    my ($type, $subdir, $maxfile, $outdir) = @_;

    my $srcdir = "$DATADIR/$subdir";
    printf "Searching %s/*.picoDst.root\n", $srcdir;
    # Use find (not a shell glob) -- large MC sets (10000s of runs) can blow
    # past the shell's argument-list limit and silently return nothing.
    my $lsdata = `find $srcdir -maxdepth 1 -name '*.picoDst.root' 2>/dev/null`;
    my @files = sort split(/\n/, $lsdata);
    my $nfile = scalar @files;
    if ($nfile == 0) {
        warn "WARNING: no picoDst files found in $srcdir\n";
        warn "  Run submit_sim.pl first and wait for jobs to finish.\n";
        exit 0;
    }
    printf "Found %d files\n", $nfile;

    # Remove old shards for this type
    system("rm -f $outdir/${type}.*.list $outdir/${type}.nsets.txt");

    my $nset = 0;
    for (my $i = 0; $i < $nfile; $i += $maxfile) {
        open(my $fh, '>', "$outdir/${type}.${nset}.list") or die;
        for (my $j = $i; $j < $i + $maxfile && $j < $nfile; $j++) {
            print $fh "$files[$j]\n";
        }
        close $fh;
        $nset++;
    }

    open(my $fh, '>', "$outdir/${type}.nsets.txt") or die;
    print $fh "$nset\n";
    close $fh;

    printf "Created %d list shards: picolist/%s.{0..%d}.list\n", $nset, $type, $nset-1;
    printf "  Next: submit_ana.pl type=%s outdir=%s submit\n", $type, $type;
}

# --------------------------------------------------------------------------
sub find_data {
    my ($maxfile, $outdir) = @_;

    my $srcpat = "$DATADIR/pico/st_*.picoDst.root";
    printf "Searching %s\n", $srcpat;
    my $lsdata = `ls $srcpat 2>/dev/null`;
    my @all_files = split /\n/, $lsdata;
    my $nfile = scalar @all_files;
    if ($nfile == 0) {
        warn "WARNING: no picoDst files found in $DATADIR/pico/\n";
        exit 0;
    }
    printf "Found %d files\n", $nfile;

    system("mv $WORKDIR/runs2.txt $WORKDIR/runs2.txt.old 2>/dev/null");

    # Group files by run number
    my (%runfiles, @runorder);
    for my $f (@all_files) {
        my $run = '';
        my $raw_idx = rindex($f, "_raw_");
        $run = substr($f, $raw_idx - 8, 8) if $raw_idx >= 8;
        next unless $run =~ /^\d+$/;
        unless (exists $runfiles{$run}) {
            push @runorder, $run;
            $runfiles{$run} = [];
        }
        push @{$runfiles{$run}}, $f;
    }

    my @sorted_runs = sort @runorder;
    open(my $rfh, '>', "$WORKDIR/runs2.txt") or die "Cannot write runs2.txt: $!\n";

    my ($nrun, $nset_total) = (0, 0);
    for my $run (@sorted_runs) {
        my @files = @{$runfiles{$run}};
        my $count = scalar @files;
        my $nsets = int(($count + $maxfile - 1) / $maxfile) || 1;
        for my $s (0 .. $nsets - 1) {
            open(my $fh, '>', "$outdir/$run.$s.lis") or die;
            for (my $i = $s*$maxfile; $i < ($s+1)*$maxfile && $i < $count; $i++) {
                print $fh "$files[$i]\n";
            }
            close $fh;
            $nset_total++;
        }
        print $rfh "$run $nsets\n";
        $nrun++;
    }
    close $rfh;

    printf "Processed %d runs, %d total list shards → picolist/ and runs2.txt\n",
        $nrun, $nset_total;
    printf "  Next: submit_ana.pl type=data outdir=OUTDIR submit\n";
}
