#!/usr/bin/perl
use strict;
die "Usage:$0 <input list><out bam status>\n" unless @ARGV==2;
open (IN,"$ARGV[0]") or die "input file requirement!\n";
open (OUT,">$ARGV[1]") or die "out file requirement\n";

my @list;
while(<IN>){
	chomp;
	push @list,$_;
}
my %hash_reads;
my %hash_mapped;
my %hash_GC;
my %hash_coverage;
my %hash_X;
my %hash_Y;
my $sample;
my %hash_base;
my %hash_quality;
foreach my $index (0..$#list){
	open (IN1,"$list[$index]");
	while (<IN1>){
		chomp;
		if(/= (.*?).marked_fixed.bam/){
			$sample = $1;
		}elsif(/number of reads = (.*?)$/){
			$hash_reads{$sample} = $1;
		}elsif(/number of mapped reads = (.*?) \((.*?)\)$/){
			$hash_mapped{$sample} = $1."\t".$2;
		}elsif(/GC percentage = (.*?)$/){
			$hash_GC{$sample} = $1;
		}elsif(/mean coverageData = (.*?)$/){
			$hash_coverage{$sample} = $1;
		}elsif(/^\s+X\s+\d+\s+\d+\s+(.*?)\s+/){
                	$hash_X{$sample} = $1;
		}elsif(/^\s+Y\s+\d+\s+\d+\s+(.*?)\s+/){
			$hash_Y{$sample} = $1;
		}elsif(/^\s+number of sequenced bases = (.*?) bp/){
			$hash_base{$sample}= $1;
		}elsif(/^\s+mean mapping quality = (.*?)$/){
			$hash_quality{$sample}= $1;
		}
	}
}
foreach my $key (sort {$a cmp $b} keys %hash_reads){
	print OUT "$key\t$hash_reads{$key}\t$hash_mapped{$key}\t$hash_base{$key}\t$hash_GC{$key}\t$hash_coverage{$key}\t$hash_X{$key}\t$hash_Y{$key}\t$hash_quality{$key}\n";
}

close IN;
close IN1;
close OUT;
