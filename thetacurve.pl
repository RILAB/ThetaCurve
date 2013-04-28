 ##########################################################################
#																									#	
#  Copyright (C) 2005-2009 Jeffrey Ross-Ibarra <rossibarra@gmail.com>  	  	#
#                                                                         	#
#  This program is free software: you can redistribute it and/or modify   	#
#  it under the terms of the GNU General Public License as published by   	#
#  the Free Software Foundation, either version 3 of the License, or      	#
#  (at your option) any later version.                                    	#
#                                                                         	#
#  This program is distributed in the hope that it will be useful,        	#
#  but WITHOUT ANY WARRANTY; without even the implied warranty of         	#
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          	#
#  GNU General Public License <http://www.gnu.org/licenses/>               #
# 	for more details.   	                                                   #
#																								  	#
 ##########################################################################

 ##########################################################################
#																									#	
#  Thanks are due to S. Wright															#
#																								  	#
 ##########################################################################


use strict;

my @p = ();
my @Q = ();
my $time = localtime(time); 
my %MLE = ();
my $bptotal;
my %multiMLE = ();

#setup file

my @lines = <>; die "no data in file\n" unless @lines;
chomp @lines;
my $loci = @lines/4;
open(OUTFILE,">outfile.txt");
print OUTFILE "Theta probability file run on $time\n\n";
my $runs = 5001;
my $multiMAXtheta = 0;

#probabilities

for my $l (1 .. $loci)
{
my $bigN = $lines[4*$l-3];
my $s = $lines[4*$l-2];
my $bp = $lines[4*$l-1]; $bptotal+=$bp;
print OUTFILE "Locus $lines[4*$l-4]\tn = $bigN\tS = $s\t$bp bp\n\nPER LOCUS THETA\tLIKELIHOOD\n\n";
my $MAXtheta = 0;

for my $t (0 .. $runs-1)
{
    my $theta = $t/100;
    
    #p2 array
    for my $f (0..$s)
    {
	$p[2]->[$f]=(1/(1+$theta))*($theta/($theta+1))**$f;
    }
    
    unless ($bigN<3)
    {
	for my $n (3..$bigN)
	{
	    for my $f (0..$s)
	    {
		$Q[$n]->[$f] = (($n-1)/($theta+$n-1)) * ($theta/($theta+$n-1))**$f;
		$p[$n]->[$f] = 0;
		for my $i (0..$f){$p[$n]->[$f] += $p[$n-1]->[$f-$i] * $Q[$n]->[$i];}
	    }
	}
    }

    print OUTFILE "$theta\t $p[$bigN]->[$s]\n";
    $MLE{$t} = $p[$bigN]->[$s];
    if($l==1){$multiMLE{$l}->{$t}=$MLE{$t};}
    else{$multiMLE{$l}->{$t}=$MLE{$t}*$multiMLE{$l-1}->{$t};}
    # sets max theta per locus
    if ($t>0 && $p[$bigN]->[$s]>$MLE{$t-1}){$MAXtheta=$theta;}
}

print OUTFILE "\nMLE of theta for $lines[4*$l-4] is $MAXtheta\n\n";
print OUTFILE '-' x 80;
print OUTFILE "\n\n";
}	

print OUTFILE "\n\n";
print OUTFILE '-' x 80;
print OUTFILE "\n\nMultilocus Likelihoods\nPER LOCUS THETA\tRELATIVE LOG LIKELIHOOD\n";
for my $a (0 .. $runs-1)
{
    my $o = $a/100;
    if ($a>0)
    {
	if ($multiMLE{$loci}->{$a}>$multiMLE{$loci}->{$a-1})
	{
	    $multiMAXtheta = $o;
	}
    }
}
my ($low, $high);
for my $j (1 .. $runs-1)
{
    my $b = $j/100; my $log;
    #print $multiMLE{$loci}->{$multiMAXtheta},"\n";
    if($multiMLE{$loci}->{$j}==0){$log=0;}
    else{
	$log = log(($multiMLE{$loci}->{$j})/($multiMLE{$loci}->{$multiMAXtheta*100}));
    }
    print OUTFILE "\n$b\t",$log;
    if ($b<=$multiMAXtheta){$low = $b unless $log>-2;}
    else{$high = ($j-1)/100 unless $log<-2;}
}
my $avbp = $bptotal/$loci;
my $multiMAX_bp=(int(0.5+10000*$multiMAXtheta/$avbp))/10000; my $lowbp = (int(0.5+10000*$low/$avbp))/10000; my $highbp = (int(0.5+10000*$high/$avbp))/10000;
print OUTFILE "\n\nMLE of the PER LOCUS theta for all loci is $multiMAXtheta (95% CI: $low - $high)\n";
print OUTFILE "MLE of the PER BP theta for all loci is $multiMAX_bp (95% CI: $lowbp - $highbp)\n";
close(OUTFILE);
