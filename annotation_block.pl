#!/usr/bin/perl -w
#===============================================================================================================
use warnings;
use strict;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"block:s","anof:s","out:s");
my $usage= <<"USAGE";
	Program: $0
	INPUT:
		-block         the blocks identified by ConsensusMatrix
		-anof          the annotation file
		-out           the output file name
USAGE
die $usage unless $opts{block};
die $usage unless $opts{anof};
die $usage unless $opts{out};
my $startTime=localtime();
#===============================================================================================================
my %bk          =();
my @idseq       =();
my $tag         =();
open I, "$opts{block}";
while (<I>)
{
	chomp;
	$_=~s/\r//g;
	if (!$_)
	{
		next;
	}
	if (/^Block\=/)
	{
		$tag=$_;
		next;
	}
	if (/^\d+\S+at/)
	{
		push(@idseq,$_);
		$bk{$_}=$tag;
	}
}
close I;
# Read gene annotation file
#===============================================================================================================
open I, "$opts{anof}";
my  %tl   =();
my  %anout=();
while (<I>)
{
	chomp;
	$_=~s/\r//g;
	if (!$_)
	{
		next;
	}
	if (/^\d+\S+at/)
	{
		my @ary=split/\t/,$_;
		$tl{$ary[0]}=$ary[21];
		if (!defined $anout{$ary[21]}) 
		{
			$anout{$ary[21]}{reprprobe} =$ary[0];
			$anout{$ary[21]}{qvalue}    =$ary[9];
			$anout{$ary[21]}{beta}      =$ary[11];
			$anout{$ary[21]}{phase}     =$ary[12];
			$anout{$ary[21]}{period}    =$ary[13];
			$anout{$ary[21]}{mean}      =$ary[10];
			$anout{$ary[21]}{geneID}    =$ary[23];
			$anout{$ary[21]}{sym}       =$ary[24];
			$anout{$ary[21]}{descrip}   =$ary[25];
			$anout{$ary[21]}{allprobe}  =$ary[0];
		}
		else
		{
			if ($anout{$ary[21]}{qvalue}>$ary[9])
			{
				$anout{$ary[21]}{allprobe}  =join ("\,",$ary[0],$anout{$ary[21]}{allprobe});
				$anout{$ary[21]}{reprprobe} =$ary[0];
				$anout{$ary[21]}{qvalue}    =$ary[9];
				$anout{$ary[21]}{beta}      =$ary[11];
				$anout{$ary[21]}{phase}     =$ary[12];
				$anout{$ary[21]}{period}    =$ary[13];
				$anout{$ary[21]}{mean}      =$ary[10];
				$anout{$ary[21]}{geneID}    =$ary[23];
				$anout{$ary[21]}{sym}       =$ary[24];
				$anout{$ary[21]}{descrip}   =$ary[25];
			}
			else
			{
				$anout{$ary[21]}{allprobe}  =join ("\,",($anout{$ary[21]}{allprobe},$ary[0]));
			}
		}
	}
}
close I;
# Output the anotation total information file
#===============================================================================================================
my @err=();
my %fag=();
open O, ">$opts{out}";
print O "TLid\tTypicalProbeID\tBlocks\tBlockPosition\tBeta\tPhase\tPeriod\tMeanExp\tq-value\tGeneID\tGeneSymbole\tGeneDescription\tAllPeriodProbes\n";
foreach my $pb (@idseq)
{
	my $id   =$tl{$pb};
	my $block=$bk{$pb};
	if (!defined $fag{$id}) 
	{
		$fag{$id}=1;
		print O join("\t",($id,$anout{$id}{reprprobe},$block,$anout{$id}{beta},$anout{$id}{phase},$anout{$id}{period},$anout{$id}{mean},$anout{$id}{qvalue},$anout{$id}{geneID},$anout{$id}{sym},$anout{$id}{descrip},$anout{$id}{allprobe})),"\n";
	}
}
if (@err)
{
	open E, ">$opts{out}\_err";
	print E "The following ProbeID has no annotation information:\n";
	my $errout=join "\t",@err;
	print E $errout,"\n";
}
close O;
#===============================================================================================================
my $options="-block  $opts{block}   -anof  $opts{anof}   -out  $opts{out}";
my $endTime=localtime();
open  LOG,">>ProgramRunning\_Log";
print LOG "From \<$startTime\> to \<$endTime\>\tperl $0 $options\n";
close LOG;
