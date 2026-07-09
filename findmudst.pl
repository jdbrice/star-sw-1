#!/usr/bin/perl

$name="/star/dataXX/reco/production_pp500_2022/ReversedFullField/P24ia/2022/*/*/st_fwd_*.MuDst.root";
$out="list";

system("rm -rf runs.txt.old");
system("mv runs.txt runs.txt.old");

system("rm -rf $out.old/*.lis");
system("rm -rf $out.old");
system("mv $out $out.old");
system("mkdir $out");

@runs=();
$nfile = 0;
$nrun = 0;
for (my $disk=19; $disk<116; $disk++){
    $search = $name;
    $search =~ s/XX/$disk/;
    printf("Searching $search\n");
    $lsdata = `ls $search`;
    @spl = split(/\n/, $lsdata);
    foreach $line (@spl){
	$run=substr($line,rindex($line,"_raw_")-8,8);
	$new=1;
	foreach $r (@runs){
	    if($run == $r) {$new=0; break;}
	}
	if($new==1){
	    push(@runs,$run);
	    $nrun++;
	}
	$nfile++;
	open(OUT, "+>> $out/$run.lis");
	print OUT "$line\n";
	close(OUT);
	printf("$line $run $nfile $nrun\n");
    }
}

@sort=sort(@runs);
open(RUN, "+> runs.txt");
foreach $run (@sort){
    print RUN "$run\n";
}
close(RUN);

printf("nFile=$nfile nrun=$nrun\n");
