#!/usr/bin/perl

my $file=$ARGV[0];
open IN,"<","$file";
while (<IN>)
{
	chomp;
	my @tmp =split(/\t/,$_);
	my ($gene, $genome)=$tmp[0]=~/(\S+?_\d+?)_(Ne\d+?_\d+?_\d+)/;
	$change_r{$genome}{$tmp[1]}+=$tmp[4];
	$change_c{$genome}{$tmp[1]}+=$tmp[5];
	$count_r{$genome}{$tmp[1]}+=$tmp[2];
	$count_c{$genome}{$tmp[1]}+=$tmp[3];
	$count{$genome}{$tmp[1]}++;
	$genome_list{$genome}=1;
}
close IN;

foreach my $key (sort keys %genome_list)
{
	$drdc_charge=$drdc_MY=();
	if ($change_c{$key}{"MY"}==0)
	{
		print "$key\tcharge\tNA\n";
		print "$key\tMY\tNA\n";
	}
	else
	{
		$drdc_charge=($change_r{$key}{"charge"}/$count_r{$key}{"charge"})/($change_c{$key}{"charge"}/$count_c{$key}{"charge"});
		$drdc_MY=($change_r{$key}{"MY"}/$count_r{$key}{"MY"})/($change_c{$key}{"MY"}/$count_c{$key}{"MY"});
		print "$key\tcharge\t$drdc_charge\n";
		print "$key\tMY\t$drdc_MY\n";
	}
}	
