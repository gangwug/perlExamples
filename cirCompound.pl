#!/usr/bin/perl -w
#===============================================================================================================
use warnings;
use strict;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"inp:s","tis:s","tls:s","out:s");
my $usage= <<"USAGE";
	Program: $0
	INPUT:
		-inp           input gene file
		-tis           tissue list
		-tls           tlid gene symbol associated list
		-out           out file name
USAGE
die $usage unless ($opts{inp} && -e $opts{inp});
die $usage unless ($opts{tis} && -e $opts{tis});
die $usage unless ($opts{tls} && -e $opts{tls});
$opts{o} = (defined $opts{out}) ? $opts{out} : "$opts{i}\_out";
my $startTime=localtime();
#===============================================================================================================
my %tis=();
open I, "$opts{tis}";
while (<I>)
{
	chomp;
	$_=~s/\r//g;
	next unless(/\S+/);
	my @ary=split /\t/,$_;
	$tis{$ary[0]}=1;
}
close I;
#===============================================================================================================
my %tls=();
open I, "$opts{tls}";
while (<I>)
{
	chomp;
	$_=~s/\r//g;
	next if ( (!$_) || (/TLid/i) );
	my @ary=split /\t/,$_;
	$tls{$ary[4]}{$ary[0]}=1;
}
close I;
#===============================================================================================================
open I, "$opts{inp}";
open L, ">$opts{out}\_loci";
open N, ">$opts{out}\_noloci";
open S, ">$opts{out}\_sym";
my $tim="";
while (<I>)
{
	chomp;
	$_=~s/\r//g;
	next unless(/\S+/);
	$tim=&tisnm($_,\%tis,"\,") if (/symbol/);
	my @ary=split /\t/,$_;
	my $out="";
	foreach my $sn (split /\,/,$tim) 
	{	$out .=$ary[$sn]."\t" if ($sn);	}
	$out=~s/\t$//;
	print S "$ary[0]\t$out\n";
	if ($ary[0]=~/symbol/) 
	{	print L "tlID\t$out\n";	}
	if (defined $tls{$ary[0]}) 
	{
		foreach my $id (sort (keys %{$tls{$ary[0]}})) 
		{	print L $id,"\t",$out,"\n";	}
	}
	else
	{	print N "$ary[0]\t$out\n";	}
}
close I;
close L;
close N;
close S;
#===============================================================================================================
open I, "$opts{out}\_sym";
my %sym=();
while (<I>)
{
	chomp;
	$_=~s/\r//g;
	next if ( (!/$_/) || (/symbol/) );
	my @ary=split /\t/, $_;
	for (my $i=1;$i<@ary;$i++) 
	{
		$ary[0]=~s/\s+//;
		$ary[$i]=~s/\s+//;
		$sym{$ary[0]}{$i}=$ary[$i];
	}
}
close I;
#===============================================================================================================
open I, "$opts{out}\_loci";
my %loc=();
while (<I>)
{
	chomp;
	$_=~s/\r//g;
	next if ( (!/$_/) || (/tlID/) );
	my @ary=split /\t/,$_;
	for (my $i=1;$i<@ary ;$i++) 
	{
		$ary[0]=~s/\s+//;
		$ary[$i]=~s/\s+//;
		$loc{$ary[0]}{$i}=$ary[$i];
	}
}
close I;
#===============================================================================================================
open L, ">$opts{out}\_loci\_sumNumber";
open S, ">$opts{out}\_sym\_sumNumber";
my $tnm=keys %tis;
my $max=&maxcom($tnm);
our $data;
for (my $j=1;$j<=(keys %tis);$j++) 
{
	@$data="";
	print L $j;
	print S $j;
	my $tag=0;
	my @arr=(1..(keys %tis));
	compound(\@arr,$j-1,"\t",0);
	$$data[0]=~s/\t\,//;
	foreach my $cop (@$data) 
	{
		$cop=~s/^\t\,//;
		if ($cop)
		{
			my @tep=split /\,/,$cop;
			print L "\t",(&comgene(\%loc,\@tep));
			print S "\t",(&comgene(\%sym,\@tep));
			$tag++;
		}
	}
	print L "\tNA" x ($max-$tag) if($max-$tag);print L "\n";
	print S "\tNA" x ($max-$tag) if($max-$tag);print S "\n";
}
close L;
close S;
#===============================================================================================================
sub tisnm{
	my ($ttl,$tis,$sar)=@_;
	my @ary=split /\t/,$ttl;
	my @tep="0";
	for (my $i=0;$i<@ary;$i++) 
	{	push @tep,$i if (defined $$tis{$ary[$i]});	}
	return (join $sar, @tep);
}
#===============================================================================================================
#==========输出组合的种类==========
sub compound{
	my ($arr,$n,$sar,$i)=@_;
	my $k=$i;
	if ($n==0)
	{
		for(my $j=$k;$j<@$arr;$j++)
		{	push @$data, $sar.",".$$arr[$j];	}
		return;
	}
	for (;$k<@$arr-$n;$k++) 
	{	&compound($arr,$n-1,$sar.",".$$arr[$k],$i=$i+1);	}
}
#==========计算组合的最大值==========
sub maxcom{
	my ($n)=@_;
	my $sub=0;
	for (my $i=1;$i<=$n;$i++) 
	{
		my $u=1;
		my $d=1;
		for (my $j=0 ;$j<$i ;$j++) 
		{
			$u *=($n-$j);
			$d *=($j+1);
		}
		$sub=$u/$d if ($u/$d > $sub);
	}
	return $sub;
}
#==========计算某种组合的基因个数==========
sub comgene{
	my ($hash,$ary)=@_;
	my $sub=0;
	foreach my $ta (sort keys %$hash) 
	{
		next if (!$ta);
		my $flag=0;
		foreach my $tb (@$ary) 
		{
			next if (!$tb);
			$flag++ if ($$hash{$ta}{$tb} eq "NA"); 
		}
		$sub++ if($flag < (@$ary) );
	}
	return $sub;
}
#==================================
my $options="-inp $opts{inp}  -tis $opts{tis}  -tls $opts{tls}  -out $opts{out}";
my $endTime=localtime();
open  LOG,">>ProgramRunning\_Log";
print LOG "From \<$startTime\> to \<$endTime\>\tperl $0 $options\n";
close LOG;
#===============================================================================================================
