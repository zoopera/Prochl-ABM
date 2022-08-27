#!/usr/bin/perl -w
###extract dN dS results from the YN00 method
  my $num_genes=$ARGV[0];
  opendir DH,"./core_genes" or die "can't opendir $!";
  my @list= grep{/fasta$/ && -f "./core_genes/$_"} readdir(DH);
  closedir DH;
 
 foreach my $file(@list)
 {
 	
	my ($famid)=$file=~/(\S+).fasta/;
	open DNDS,">./dNdS/$famid.dNdS" or die "can't open";
	print DNDS "gene1\tgene2\tN\tS\tdN\tdS\n";
	chomp $file;
	print $file."\n";
	open IN,"./core_genes/$file" or die "can't open $!";	
####transfer as phylip format
	open OUT,"> YN00/seqs.phylip";
	my $len=0;
	while(<IN>)
	{
	  chomp;
	 if(/>(\S+\|\S+)$/ && $len == 0)
	 {
	  my $name =$1;
	  my $line=<IN>;
	  $line=~s/\s//g;
	  $len=length($line);
	  print OUT "$num_genes $len\n";
	  print OUT "$name\n$line\n";
	 }	
	 elsif(/>(\S+\|\S+)$/ && $len != 0)
	 {
	  my $name=$1;
	  my $line=<IN>;
	  $line=~s/\s//g;
	  print OUT "$name\n$line\n";
	 }
	}
	
	close IN;
	close OUT;
#####run YN00 
	chdir("./YN00");
	system("/share/software/paml4.9e/src/yn00 yn00.ctl");	
	chdir("../");
#### extract YN00 output	
	$output = $famid.".out";
	`mv ./YN00/yn.out ./YN00/$output`;
	open YIN,"./YN00/$output" or die "can't open $!";
	my @dnds;
	my @names;
	while(<YIN>)
	{
	  if(/^\s+\d+\s+\d+\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+\S+\s+(\S+)\s+\+\-\s+\S+\s+(\S+)\s+\+\-/) 
	  {
		push @dnds,"$2\t$1\t$3\t$4";   ### keep the number of sits and dN, dS in the same order
	  }
	  if(/vs/)
	  {
		$_=~/\((\S+?)\|\d+\).*?\((\S+?)\|\d+\)/;  ###modify the regular expression based on you data
		push @names,"$1\t$2";
	  }
	}
	close YIN;
	for(my $i=0;$i<=$#names;$i++)
	{
		print DNDS $names[$i]."\t".$dnds[$i]."\n";
	}
 }
	close DNDS;
