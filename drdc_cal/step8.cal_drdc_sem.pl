#!/usr/bin/perl
#
my $file = $ARGV[0];
open IN,"<","$file";
while (<IN>)
{
	chomp;
	my @tmp = split(/\t/,$_);
	if ($tmp[1] eq "MY")
	{
		my ($aa)=$tmp[0]=~/(\S+?_\d+?)_\d+/;
		$name = "$aa"."_MY";
		push @{$name}, $tmp[2];
		$list{$aa}=1;
	}
	elsif ($tmp[1] eq "charge")
	{
		my ($aa)=$tmp[0]=~/(\S+?_\d+?)_\d+/;
		$name = "$aa"."_charge";
                push @{$name}, $tmp[2];
		$list{$aa}=1;
	}
}
close IN;

foreach my $key (sort keys %list)
{
	my $sum =$avg=$avg2=0;
	my $q1 = "$key"."_MY";
	my $q2 = "$key"."_charge";
	
	my $s1 = 0;
	grep {$s1 += $_}@{$q1};
	my $avg1 = $s1/($#{$q1}+1);
	my $d1 = 0;
	grep {$d1 += ($_-$avg1)**2;}@{$q1};
	my $tmp1 = $d1/($#{$q1}+1);
	my $sd1 = sqrt($tmp1);
	my $sem1 = $sd1/sqrt($#{$q1}+1);
	print "$key\tMY\t$avg1\t$sem1\n";
		
	my $s2=0;
	grep {$s2 += $_}@$q2;
	my $avg2 = $s2/($#{$q2}+1);
        my $d2 = 0;
        grep {$d2 += ($_-$avg2)**2;}@{$q2};
        my $tmp2 = $d2/($#{$q2}+1);
        my $sd2 = sqrt($tmp2);
        my $sem2 = $sd2/sqrt($#{$q2}+1);
        print "$key\tcharge\t$avg2\t$sem2\n";	
}

