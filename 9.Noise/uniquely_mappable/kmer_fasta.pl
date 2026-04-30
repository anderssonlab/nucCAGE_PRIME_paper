#!/usr/bin/perl

use strict;
use warnings;

my $usage = "Usage: $0 <infile.fa> <kmer length> <bin size> <out directory>\n";
my $infile = shift or die $usage;
my $read_length = shift or die $usage;
my $bin_size = shift or die $usage;
my $out = shift or die $usage;

my $chr = "";
my $seq = "";

open(IN,'<',$infile) || die "Could not open $infile: $\n";
my $iter = 0;
my $bin = 1;

open(OUT,">",$out."/bin.".$bin.".fa") || die "Could not open ".$out."/bin.".$bin.".fa for writing: $\n";

while (<IN>) {
	chomp;
	next if /^$/;
	if (/^>(\S*).*$/){
		if ($chr ne "")
		{
			my $length = length($seq);
			my $end = $length - $read_length;
			for (my $i = 0; $i <= $end; ++$i) {
				my $read = substr($seq,$i,$read_length);
				$read = uc($read);
				
				next if $read =~ /N/;
				die "$read\n" unless $read =~ /[ATGC]+/;

				my $chr_start = $i + 1;
				my $chr_end = $i + $read_length;
				my $id = $chr . ':' . $chr_start . '-' . $chr_end;

				print OUT ">$id\n$read\n";
				
				$iter++;
				if ($iter == $bin_size)
				{
					close(OUT);
					$iter = 0;
					$bin++;
					open(OUT,">",$out."/bin.".$bin.".fa") || die "Could not open ".$out."/bin.".$bin.".fa for writing: $\n";
				}
			}
		}
		$chr = $1;
		$seq = "";
		print STDERR "$chr\n";
	} else
	{
		$seq .= $_;
	}
}
close(IN);

my $length = length($seq);
my $end = $length - $read_length;
for (my $i = 0; $i <= $end; ++$i) {
	my $read = substr($seq,$i,$read_length);
	$read = uc($read);
	
	next if $read =~ /N/;
	die "$read\n" unless $read =~ /[ATGC]+/;
	
	my $chr_start = $i + 1;
	my $chr_end = $i + $read_length;
	my $id = $chr . ':' . $chr_start . '-' . $chr_end;
	
	print OUT ">$id\n$read\n";
	
	$iter++;
	if ($iter == $bin_size)
	{
		close(OUT);
		$iter = 0;
		$bin++;
		open(OUT,">",$out."/bin.".$bin.".fa") || die "Could not open ".$out."/bin.".$bin.".fa for writing: $\n";
	}
}

print STDERR "Done!\n";
close(OUT);

exit(0);
