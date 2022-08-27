#!/usr/bin/env perl

$dir = $ARGV[0];
opendir DIR,"$dir";
my @list = readdir DIR;
my @list_grep = grep /Gene.*_out/,@list;
foreach my $input (@list_grep)
{
	my ($comp, $class)=$input=~/(\S+?)\.(\S+)_out/;
	open IN,"<","./$dir/$input";
	while (<IN>)
	{
		chomp;
		if (/^otu  1:\s+(\S+?)\s+(\S+)/)
		{
			$count_r = $1;
			$count_c = $2;
		}
		elsif (/^otu  2:\s+(\S+?)\s+(\S+)/)
		{
			$count_r2 = $1;
			$count_c2 = $2;
		}
		elsif (/Numbers of radical \(r, above diagonal/)
		{
			$count_r_avg = ($count_r+$count_r2)/2;
			$count_c_avg = ($count_c+$count_c2)/2;
			
			my $line_1 = <IN>;
			my @tmp = split(/\s+/,$line_1);
			my $line_2 = <IN>;
			my @tmp2 = split(/\s+/,$line_2);
			$change_r = $tmp[2];
			$change_c = $tmp2[1];
		}
	}
	print "$comp\t$class\t$count_r_avg\t$count_c_avg\t$change_r\t$change_c\n";
	close IN;
}
closedir DIR;

#print "$comp\t$class\t$count_r_avg\t$count_c_avg\t$change_r\t$change_c\n";
