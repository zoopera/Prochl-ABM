#!/usr/bin/perl

$file = $ARGV[0];
open IN,"<","$file";
while (<IN>)
{
        chomp;
        my @tmp= split (/\t/,$_);
        if ($tmp[1]=~/MED4/)
        {
				$count_n = $tmp[3]*$tmp[5];
				$count_s = $tmp[4]*$tmp[6];
				$site_n{$tmp[2]}+=$tmp[3];
				$site_s{$tmp[2]}+=$tmp[4];
				$sum_n{$tmp[2]}+=$count_n;
				$sum_s{$tmp[2]}+=$count_s;
				$list{$tmp[2]}++;
        }
	elsif ($tmp[2]=~/MED4/)
        {
		$count_n = $tmp[3]*$tmp[5];
                $count_s = $tmp[4]*$tmp[6];
                $site_n{$tmp[1]}+=$tmp[3];
                $site_s{$tmp[1]}+=$tmp[4];
                $sum_n{$tmp[1]}+=$count_n;
                $sum_s{$tmp[1]}+=$count_s;
		$list{$tmp[1]}++;
        }
}
close IN;

foreach my $key (sort keys %list)
{
	$dnds = ($sum_n{$key}/($site_n{$key}/$list{$key}))/($sum_s{$key}/($site_s{$key}/$list{$key}));
	
#	$dnds = ($dn{$key}/$count{$key})/($ds{$key}/$count{$key});
#       $avg = $sum_dnds{$key}/$count_dnds{$key};
        print "MED4_$key\t$dnds\n";
}
