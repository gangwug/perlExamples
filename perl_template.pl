#!/usr/bin/perl -w
#A:open file
#===============================================================================================================
use warnings;
use strict;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i:s","o:s");
my $usage= <<"USAGE";
	Program: $0
	INPUT:
		-i           inFile
		-o           outFile
	FILES REQUIRED:
	OUTPUT:
		outFile
USAGE
die $usage unless ($opts{i} && -e $opts{i});
$opts{o} = (defined $opts{o}) ? $opts{o} : "$opts{i}\_out";
my $startTime=localtime();
#===============================================================================================================
open I, "$opts{i}";
while (<I>)
{
	chomp;
	$_=~s/\r//g;
	next unless(/\S+/);
	my @t=split(/\t/,$_);
}
close I;
#===============================================================================================================
open O,">$opts{o}";
close O;
#===============================================================================================================
my $options="-i $opts{i}  -o $opts{o}";
my $endTime=localtime();
open  LOG,">>ProgramRunning\_Log";
print LOG "From \<$startTime\> to \<$endTime\>\tperl $0 $options\n";
close LOG;
#===============================================================================================================
################################################################################################################
#===============================================================================================================
#!/usr/bin/perl -w
#B:open dir
use warnings;
use strict;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"indir:s","outdir:s");
my $usage= <<"USAGE";
	Program: $0
	INPUT:
		-indir           the box number in every loci matrix
		-outdir           the output file dir
USAGE
die $usage unless $opts{indir};
die $usage unless $opts{outdir};
my $startTime=localtime();
# Read the dir including box_site files
#===============================================================================================================
opendir (IN,"$opts{indir}")||die "$!";
while (my $infile=readdir(IN))
{

	next if ($infile!~/\w/);
	open I, "$opts{indir}/$infile";
	my %allpos       =();
	my %segpos       =();
	my $tlname       =();
	my $outname  ="$infile"."_segmatrix";
	if ($infile=~/noperiod/)
	{	$tlname="noperiod";	}
	else
	{	$tlname="period";	}
	while (<I>)
	{
		chomp;
		$_=~s/\r//g;
		next if ( (!$_) || (/TLid/) )
	}
	close I;
	open OUT, ">$opts{outdir}/$outname";
	print OUT $tlname,$title,"\n";
	foreach my $segid (sort keys %segpos)
	{
		print OUT $segid;
		foreach my $spos ( sort {$a<=>$b} (keys %{$segpos{$segid}}) ) 
		{	print OUT "\t$segpos{$segid}{$spos}";	}
		print OUT "\n";
	}
	close OUT;
}
close IN;
#===============================================================================================================
my $options="-indir  $opts{indir}  -outdir $opts{outdir}";
my $endTime=localtime();
open  LOG,">>ProgramRunning\_Log";
print LOG "From \<$startTime\> to \<$endTime\>\tperl $0 $options\n";
close LOG;
#===============================================================================================================
################################################################################################################
#===============================================================================================================

