#!/usr/bin/perl
$minjobs=0;
$speedup=100;
$wait=300;

$numjobs=$minjobs;
while (1){
    $status = `condor_q $ENV{'USER'} | grep "runpico_uni data" | wc -l`;
    print "running=",$status;
    if ($status ne "0") {
        $numjobs = `echo "$status"`;
        # print $numjobs;
        if ($numjobs<=$speedup) {$wait=60;}
        if ($numjobs<=$minjobs) {exit;}
	sleep $wait;     
    } else {
	exit;
    }
}
