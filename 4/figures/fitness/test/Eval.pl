#!/usr/bin/perl
#operating system dependant

$POP = 590;    # population
$outf="./new_fit.fit";
$pat=100; #trainign pats
$maxF=4.3*65688;
if ($^O =~ /Win32/ ) #not checked
{

}
else #checked
{	
	open (F, $outf) || die "Could not open $filename: $!";
	while ($line = <F>) {
		chomp($line);
		next if $line =~ /^\s*$/; # skip blank lines
		if ($line =~ /^([A-Za-z]\w*)/) {
			$out = $1;
		} 
		else {
			my (@row) = split (/\s+/, $line);
			push (@{$out}, \@row); # insert the row-array into
        	}
	}
	close(F);
	
	open (MYFILE, './task.dat');
	$nit=0;
	while (<MYFILE>) {
		chomp;
		$task->[$nit]="$_";
		#print  $task->[$nit]."\n";
		$nit=$nit+1;
	}
	close (MYFILE);
	
#	open (F, "task.dat") || die "Could not open $filename: $!";
#	while ($line = <F>) {
#		chomp($line);
#		next if $line =~ /^\s*$/; # skip blank lines
#		if ($line =~ /^([A-Za-z]\w*)/) {
#			$task = $1;
#		} 
#		else {
#			my (@row) = split (/\s+/, $line);
#			push (@{$task}, \@row); # insert the row-array into
#        	}
#	}
#	close(F);

        
		open OUTPUT, ">task.res" or die $!;
	for($count10=0;$count10<$pat ;$count10++) #POP loop
	{
		$mindist=10000000000000; 	
		$res1=1000;
		$pos=-10;
		for($count=0;$count<$POP ;$count++) #POP loop
		{
			$dist=0;
			for($count1=0; $count1<2; $count1++){ 
			 $dist=$dist+($out->[$count][$count1]-$task->[$count1+1])*($out->[$count][$count1]-$task->[$count1+1]);
			}
			#print  $count." ".$dist."\n";
			if($dist<$mindist){
				$mindist=$dist;
				$res1=$out->[$count][2];
				$res2=$out->[$count][3];
				$pos=$count;
			}
		}
		#PCA
		$xrot=0.961455*($out->[$pos][0]-0.2283)+0.274961*($out->[$pos][1]-0.449131);
		$yrot=-0.274961*($out->[$pos][0]-0.2283)+0.961455*($out->[$pos][1]-0.449131);		
		print  $xrot." ".$yrot." ".$res1." ".$res2."\n" ;
		print OUTPUT $xrot." ".($res1/$maxF)."\n" ;
		#NoPCA
		#print  $out->[$pos][0]." ".$out->[$pos][1]." ".$res1." ".$res2."\n" ;
		#print OUTPUT $out->[$pos][0]." ".$out->[$pos][1]." ".$res1." ".$res2."\n" ;
		$out->[$pos][0]=1000000000;
		$out->[$pos][1]=1000000000;
	}
	
	$xrot=0.961455*($task->[1]-0.2283)+0.274961*($task->[2]-0.449131);
	$yrot=-0.274961*($task->[1]-0.2283)+0.961455*($task->[2]-0.449131);		
	print  "Rotated task.dat".$xrot." ".$yrot."\n" ;
	
	close OUTPUT;
}

###################################################################################
# everything worked fine (i hope), so exit with code 0
##################################################################################
