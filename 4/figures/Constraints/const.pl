#!/usr/bin/perl
#operating system dependant

$nx1 = 80;    # population
$nx2 = 80;    # population
$x1min = 0.2;    # population
$x1max = 8.0;    # population
$x2min = 0.2;    # population
$x2max = 8.0;    # population
$outf="./task.res";

if ($^O =~ /Win32/ ) #not checked
{

}
else #checked
{	
	open OUTPUT1, ">all.res" or die $!;
	open OUTPUT2, ">const_tx.fail" or die $!;
	open OUTPUT3, ">const_tx.feasible" or die $!;
	open OUTPUT4, ">const_sx.fail" or die $!;
	open OUTPUT5, ">const_sx.feasible" or die $!;
	open OUTPUT6, ">constALL.fail" or die $!;
	open OUTPUT7, ">constALL.feasible" or die $!;
	open OUTPUT8, ">const_P.fail" or die $!;
	open OUTPUT9, ">const_P.feasible" or die $!;
	open OUTPUT10, ">const_Price.fail" or die $!;
	open OUTPUT11, ">const_Price.feasible" or die $!;
	open OUTPUT12, ">const_Dx.fail" or die $!;
	open OUTPUT13, ">const_Dx.feasible" or die $!;
	for($count=0;$count<$nx1 ;$count++) #POP loop
	{
		$x1=$x1min+($x1max-$x1min)*$count/($nx1-1); 
		for($count1=0;$count1<$nx2 ;$count1++) #POP loop
		{
			$x2=$x2min+($x2max-$x2min)*$count1/($nx2-1);
			system("rm task.dat"); 
			system("rm task.res"); 
			system("rm task.cns");
 
			open OUTPUT, ">task.dat" or die $!;
			print OUTPUT "3 \n 0.6 \n".$x1."\n".$x2."\n"; 		
			close OUTPUT;
			
			system("./weldbeam.exe");
			
			open (MYFILE, './task.cns');
			$nit=0;
			while (<MYFILE>) {
				chomp;
				$task->[$nit]="$_";
				#print  $task->[$nit]."\n";
				$nit=$nit+1;
			}
			close (MYFILE);
			print OUTPUT1 $x1." ".$x2." ";
			$check=-1;
			for($cc=0;$cc<$nit ;$cc++) #POP loop
			{
				if($task->[$cc]>0){$check=1;}
				print OUTPUT1  $task->[$cc]." ";
			}
			print OUTPUT1 "\n";
			#Const1
			if($task->[0] > 0) {
			 	print OUTPUT2 $x1." ".$x2." \n";
			}
			else{
			 	print OUTPUT3 $x1." ".$x2." \n";
			}
			#Const2
			if($task->[1] > 0) {
			 	print OUTPUT4 $x1." ".$x2." \n";
			}
			else{
			 	print OUTPUT5 $x1." ".$x2." \n";
			}
			#ConstALL
			if($check > 0) {
			 	print OUTPUT6 $x1." ".$x2." \n";
			}
			else{
			 	print OUTPUT7 $x1." ".$x2." \n";
			}
			#Const7
			if($task->[6] > 0) {
			 	print OUTPUT8 $x1." ".$x2." \n";
			}
			else{
			 	print OUTPUT9 $x1." ".$x2." \n";
			}
			#Const4
			if($task->[3] > 0) {
			 	print OUTPUT10 $x1." ".$x2." \n";
			}
			else{
			 	print OUTPUT11 $x1." ".$x2." \n";
			}
			#Const6
			if($task->[5] > 0) {
			 	print OUTPUT12 $x1." ".$x2." \n";
			}
			else{
			 	print OUTPUT13 $x1." ".$x2." \n";
			}

		} 		
		print OUTPUT1 "\n";

	}
	close OUTPUT1;
	close OUTPUT2;
	close OUTPUT3;
	close OUTPUT4;
	close OUTPUT5;
	close OUTPUT6;
	close OUTPUT7;
	close OUTPUT8;
	close OUTPUT9;
	close OUTPUT10;
	close OUTPUT11;
	close OUTPUT12;
	close OUTPUT13;
}

###################################################################################
# everything worked fine (i hope), so exit with code 0
##################################################################################
