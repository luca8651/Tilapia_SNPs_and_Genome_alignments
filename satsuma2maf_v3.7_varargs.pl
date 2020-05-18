#!/usr/bin/perl -w
use strict;
use warnings;
## line 180 has changed:
#if ( $cont % 2 == 0) {  $seq1 = $seq1 . $line ; }    #0 rather than 1 (target and then query!)
### condition "$length11==$length21" added for an alignment to be printed 

my $lines =shift;   # satsuma alignment file
my $lengthfile1=shift;
my $lengthfile2=shift;
my $queryname=shift;
my $targetname=shift;
my $output_name=shift;

if ((!$lines) || (!$lengthfile1)){
print STDERR "$0: A script to convert Satsuma2 alignment format to MAF\n";
print STDERR "$0 <MergeXCorrMatches.chained.out> <query_chr_lengths_file> <target_chr_lengths_file> <query_species_name> <target_species_name> [<output_file_name>]\n";
print STDERR "$0 chr_lengths_file format: tab delimited; column 1 indicates chromosome names, column 2 the corresponding lengths\n";
exit();
}




my $cont=0;

my @array;
my @array2;

my $chr;
my $strand;
my $start;
my $end;
#my $strand;
my $plus="+";
my $start1;
my $start1copy;
my $start2;
my $end1;
my $length1;
my $length11;
my $length2;
my $length21;
my $seq1;
my $seq2;


my $IDcont=0;
my $line2;
my $line3;
my $query;
my $target;
my $strand2="+";
my $strand1;
my %Larray1;
my %Larray2;
my $name1;
my $name2;
my $l1;
my $score;
my $AMBcont=0;
my $scorecont=0;

my $test1;
my $test2;
my $test3;
my $test4;
#Lines starting with "s" contain: the sequence name, the start coordinate of the alignment, the number of bases spanned by the alignment, the strand, the sequence length, and the aligned bases.
open(FILE1, "$lengthfile1");

if ((!$output_name)){
#open my $output1, '>', "$queryname.$targetname.correct" or die $!;
#open my $output2, '>', "$queryname.$targetname.missing_start" or die $!;
#open my $output3, '>', "$queryname.$targetname.out_of_range" or die $!;
$output_name="$queryname.$targetname";
}
#else                {
open my $output1, '>', "$output_name.correct" or die $!;
open my $output2, '>', "$output_name.missing_start" or die $!;
open my $output3, '>', "$output_name.out_of_range" or die $!;
#                    }

#open my $output1, '>', "$queryname.$targetname.correct" or die $!;
#open my $output2, '>', "$queryname.$targetname.missing_start" or die $!;
#open my $output3, '>', "$queryname.$targetname.out_of_range" or die $!;



while ( my $line11=<FILE1>)  {

@array=split(/\t/, $line11);
$name1=$array[0];
$l1=$array[1];
$l1=substr($l1,0,length($l1)-1);
$l1 =~ s/ //g;
$Larray1{$name1} = $l1;
#print "query $Larray1{$name1}\n";

                          }


close($lengthfile1);

open(FILE2, "$lengthfile2");

while ( my $line12=<FILE2>)  {

@array=split(/\t/, $line12);
$name2=$array[0];
$l1=$array[1];
$l1=substr($l1,0,length($l1)-1);
$l1 =~ s/ //g;
$Larray2{$name2} = $l1;
#print "target $Larray2{$name2}\n";
 
                             }


close($lengthfile2);

print $output1 "##maf version=1 scoring=Satsuma2_perc_identity\n# Satsuma2 $lines\n";

open(FILE, "$lines");

while ( my $line=<FILE>)  {
#while ( my $line ~= $lines )  {

if ( $line=~/AMB/ || $line=~/dentity/ || $line=~/overlapping/ || $line=~/ERROR/) { $AMBcont=0;  }   ##the 3rd and 5th line after finding "AMB" should contain alignments. Matching "dentity" makes sure that $AMBcont is zero at the beginning of every alignment
else { $AMBcont++ ; }

if ( $line=~/dentity/  ) { $scorecont=0;  }   ##the 3rd and 5th line after finding "dentity" should contain alignments. These cannot be identified using AMBcont
else { $scorecont++ ; }




if ( $line=~/Query/ && $line=~/target/ && $line=~/check/ )                    {          #Query scaffold_10 [62-434] vs target lg1 [11671135-11671504] + length 369 check 369




 if ( $IDcont>0  )                                          		                   	{
  $line2=$seq1;
  $line2 =~ s/-//g; 
  $length1=length($line2);
  $line2=$seq2;
  $line2 =~ s/-//g; 
  $length2=length($line2);

  $length11=length($seq1);
  $length21=length($seq2);

  #$test1=0;
  #$test1 =~ s/-//g;
  #$test2=$seq2;
  #$test2 =~ s/-//g;
  

  $test2=($start2+$length2-1)<=$Larray2{$target};
  $test3=1;  

  if (  length($start1copy)==0  ) {
  $test1=-1;               }
  else			  {
  $test1=$start1copy+$length1-1<=$Larray1{$query}; #$start2<=$Larray1{$query};                              #start1 gives problems
  $test3=$start1copy>0; 
			  }
  #$test3=$start1; #=~/^[0-9]+$/;
  #$test3 =~ s/[0-9]//g;
  #print "ssstart1 is $start1, array entry is $Larray1{$query}; $start1+$length1-1<=$Larray1{$query} \n";
  #$test1=($start1+$length1-1)<=$Larray1{$query};
  #$test3=0; #$start1+$length1-1>0;
  #$test3=$start1+$length1-1>0;
  #$test4=$Larray1{$query}>0;
  #$test2=$start2+$length2-1;
  #$test1=$start1+$length1-1;
  #print "length11 $length11, length21 $length21, start1 is $start1, length1 is $length1, start2 is $start2, length2 is $length2, $Larray1{$query} and $Larray2{$target}\n";
  print "heey $start1copy $query $start1 boh $test3\n";

   if   (   $test1==0 || $test2==0 || $start1<0 || $start2<0 )                                                   {  #(  $start1copy<=0 ||   $test1==0 || $test2==0 )

  $query=~ s/ /_/g;
  $target=~ s/ /_/g;
  print $output3 "a score=$score\n";
  print $output3 "s $queryname.$query $start1 $length1 $strand1 $Larray1{$query} $seq1\n";
  print $output3 "s $targetname.$target $start2 $length2 $strand2 $Larray2{$target} $seq2\n\n";


   
											}
elsif ( length($start1copy)==0 || length($start2)==0   ) 	              		        { #problematic alignment    

  #print "neu $test1 and $test2 and $test3 and $test4; $length11 and $length21, start1 is $start1, length1 is $length1, start2 is $start2, length2 is $length2, $Larray1{$query} and $Larray2{$target}\n";
  print $output2 "a score=$score\n";
  print $output2 "s $queryname.$query $start1 $length1 $strand1 $Larray1{$query} $seq1\n";
  print $output2 "s $targetname.$target $start2 $length2 $strand2 $Larray2{$target} $seq2\n\n";



											}
 
  elsif ( $length11==$length21 && $seq1 =~ /\A[-ATGCNatgcn]+\Z/ && $seq2 =~ /\A[-ATGCNatgcn]+\Z/  && ( $start1copy+$length1-1<=$Larray1{$query} ) && ( $start2+$length2-1<=$Larray2{$target} ) )                                             			{      #requires the alignment coords to be in range with chr sizes
 
 
  if (length($start1copy)>0 && length($start2)>0)  	{

  $query=~ s/ /_/g;
  $target=~ s/ /_/g;
  print $output1 "a score=$score\n";
  print $output1 "s $queryname.$query $start1 $length1 $strand1 $Larray1{$query} $seq1\n";
  print $output1 "s $targetname.$target $start2 $length2 $strand2 $Larray2{$target} $seq2\n\n";


			   			}
											}
			   



                                                                               }

@array=split(/ /, $line);
$query=$array[1];    #??
$target=$array[5];      #??
$strand1=$array[7];#$strand2=$array[7];

$start1=$array[2];
@array2=split(/-/, $start1);
$start1copy=$array2[0];
$start1copy =~ s/\[//g;
$start1copy =~ s/\]//g; 

if ($strand1=~/\+/  )	{
$start1=$array2[0];
$start1 =~ s/\[//g;
$start1 =$start1-1;		#in a 0-based coord system, the start of a fwd alignment is the original start-1 
			}
else                    {
$start1=$array2[1];
$start1 =~ s/\[//g;
$start1 =~ s/\]//g;            #
#$start1 =$start1-1;
$start1 =$Larray1{$query}-$start1; #+1-1
			}

#print "heyhey $strand1 $start1copy $start1 $array2[0] $array2[1]\n";
$start2=$array[6];                                           
@array2=split(/-/, $start2);
$start2=$array2[0];
$start2 =~ s/\[//g;
$start2 =$start2-1;


										


      


 

                                                      }

if ( $line=~/dentity/   )		{

@array=split(/ /, $line);
$score=$array[4];
$cont=1;
$IDcont++;
#print $line;
if  ($IDcont>1 )   				{

#print "$Larray1{$query}\n";
#print "$Larray2{$target}\n";
#print "a score=$score\n";

                   				}



$seq1="";
$seq2="";
                                      	}

$line=substr($line,0,length($line)-1); #remove newline
$line2=$line;
$line2 =~ s/[ATGCN]//g;
$line2 =~ s/-//g;
$line3 =$line; 
$line3 =~ s/-//g;   
@array=split(/\t/, $line);


if (  ( (length($line2)<2 && length($line3)>1 ) ||  $AMBcont==3 || $AMBcont==5 || $scorecont==3 || $scorecont==5  )  &&  ($line=~/AMB/)==0  && ($line=~/ERROR/)==0 )   {    #matches lines with aligned sequences, and identifies gap-only lines based on $AMBcont and $scorecont. 
$cont++;
#print "$cont and $line\n";
if ( $cont % 2 == 0) {  $seq1 = $seq1 . $line ; }          #cont increases when a seq is found. The two genomes alternate
else  {  $seq2 = $seq2 . $line ; }



                          }


    if (eof && $IDcont<2 )                                                              {          ##if only one alignment is found

	$line2=$seq1;
  	$line2 =~ s/-//g; 
  	$length1=length($line2);
  	$line2=$seq2;
  	$line2 =~ s/-//g; 
 	$length2=length($line2);

        $length11=length($seq1);
        $length21=length($seq2);
        
        $test2=($start2+$length2-1)<=$Larray2{$target};
	$test3=1;  

  if (  length($start1copy)==0) {
  $test1=-1;               }
  else			  {
  $test1=$start1copy+$length1-1<=$Larray1{$query}; #$start2<=$Larray1{$query};                              #start1 gives problems
  $test3=$start1copy>0; 
			  }
        #$test1=$start1+$length1-1<=$Larray1{$query}; #$start2<=$Larray1{$query};                              #start1 gives problems
        #$test3=$start1+$length1-1>0;
        #$test4=$Larray1{$query}>0;
print "heey $start1copy $query $start1 boh $test3\n";

	if   (   $test1==0 || $test2==0 || $start1<0 || $start2<0 )                                               {

  $query=~ s/ /_/g;
  $target=~ s/ /_/g;
  print $output3 "a score=$score\n";
  print $output3 "s $queryname.$query $start1 $length1 $strand1 $Larray1{$query} $seq1\n";
  print $output3 "s $targetname.$target $start2 $length2 $strand2 $Larray2{$target} $seq2\n\n";


   
											}
	elsif ( length($start1)==0 || length($start2)==0    ) 	              		{ #problematic alignment    

  #print "$test1 and $test2 and $test3 and $test4; $length11 and $length21, start1 is $start1, length1 is $length1, start2 is $start2, length2 is $length2, $Larray1{$query} and $Larray2{$target}\n";
  print $output2 "a score=$score\n";
  print $output2 "s $queryname.$query $start1 $length1 $strand1 $Larray1{$query} $seq1\n";
  print $output2 "s $targetname.$target $start2 $length2 $strand2 $Larray2{$target} $seq2\n\n";



											}


  	elsif (  $length11==$length21 && $seq1 =~ /\A[-ATGCNatgcn]+\Z/ && $seq2 =~ /\A[-ATGCNatgcn]+\Z/ && $start1+$length1-1<=$Larray1{$query} && $start2+$length2-1<=$Larray2{$target})                                             				{


       
	if (length($start1)>0 && length($start2)>0 )  {	

        $query=~ s/ /_/g;
        $target=~ s/ /_/g;
	print $output1 "a score=$score\n";
	print $output1 "s $queryname.$query $start1 $length1 $strand1 $Larray1{$query} $seq1\n";
        print $output1 "s $targetname.$target $start2 $length2 $strand2 $Larray2{$target} $seq2\n\n";	
                                 } 
											}

											}

	elsif (eof)  										{       #end of file


	$line2=$seq1;
  	$line2 =~ s/-//g; 
  	$length1=length($line2);
  	$line2=$seq2;
  	$line2 =~ s/-//g; 
 	$length2=length($line2);

        $length11=length($seq1);
        $length21=length($seq2);
	$test2=($start2+$length2-1)<=$Larray2{$target};
	$test3=1;  

  if (  length($start1copy)==0) {
  $test1=-1;               }
  else			  {
  $test1=$start1copy+$length1-1<=$Larray1{$query}; #$start2<=$Larray1{$query};                              #start1 gives problems
  $test3=$start1<0; #$start1copy>0;
  #print "test $test3\n"; 
			  }
	
print "heey $start1copy $query $start1 boh $test3\n";
	if   (  $test1==0 || $test2==0 ||  $start1<0 || $start2<0 )                                               {     #$test3==0

  $query=~ s/ /_/g;
  $target=~ s/ /_/g;
  print $output3 "a score=$score\n";
  print $output3 "s $queryname.$query $start1 $length1 $strand1 $Larray1{$query} $seq1\n";
  print $output3 "s $targetname.$target $start2 $length2 $strand2 $Larray2{$target} $seq2\n\n";


   
											}
	elsif ( length($start1)==0 || length($start2)==0     ) 	              		{ #problematic alignment    

  #print " $test1 and $test2 and $test3 and $test4; $length11 and $length21, start1 is $start1, length1 is $length1, start2 is $start2, length2 is $length2, $Larray1{$query} and $Larray2{$target}\n";
  print $output2 "a score=$score\n";
  print $output2 "s $queryname.$query $start1 $length1 $strand1 $Larray1{$query} $seq1\n";
  print $output2 "s $targetname.$target $start2 $length2 $strand2 $Larray2{$target} $seq2\n\n";



											}
 	elsif ( $length11==$length21 && $seq1 =~ /\A[\-ATGCNatgcn]+\Z/ && $seq2 =~ /\A[\-ATGCNatgcn]+\Z/  && $start1+$length1-1<=$Larray1{$query} && $start2+$length2-1<=$Larray2{$target}  )                                       					{
        if (length($start1)>0 && length($start2)>0 )  {

        $query=~ s/ /_/g;
        $target=~ s/ /_/g;
        print $output1 "a score=$score\n";
	print $output1 "s $queryname.$query $start1 $length1 $strand1 $Larray1{$query} $seq1\n";
        print $output1 "s $targetname.$target $start2 $length2 $strand2 $Larray2{$target} $seq2\n\n";
	
						       }
											}



												}

                                        }



		#print  "$chr\t$UTRstart\t$UTRend\tthree_prime_utr\t$query1\_$transcript\t$strand\n";
                #print  "$chr\t$UTRstart\t$UTRend\t$query1\tthree_prime_utr\t$strand\n";
