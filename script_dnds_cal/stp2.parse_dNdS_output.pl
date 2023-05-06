#!/usr/bin/env perl

my $file =$ARGV[0];
my ($id)=$file=~/(\S+)\.dNdS/;
open IN,"<","$file";
open OUT,">>","Summary_dnds.tbl";
my $head = <IN>;
while (<IN>)
{
	chomp;
	my @tmp = split (/\t/,$_);
	print OUT "$id\t$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$tmp[4]\t$tmp[5]\t";
	if ($tmp[5]>0)
	{
		$dnds = $tmp[4]/$tmp[5];
		print OUT "$dnds\n";
	}
	else
	{
		print OUT "NA\n";
	}
}
close IN;
close OUT;

