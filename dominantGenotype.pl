#!/usr/bin/perl
use strict;
use warnings;

die "Usage: $0 <in gzvcf> <out gzvcf>" unless @ARGV == 2;

# Open the input gzipped VCF file for reading
open(IN, "gzip -dc $ARGV[0] |") or die "Cannot open gzvcf file";

# Open a pipe to `bcftools view` to compress output directly to a gzipped VCF
open(OUT, "| bcftools view -Oz -o $ARGV[1]") or die "Cannot open bcftools pipe";

while (<IN>) {
    chomp;
    if (/^#/) {
        print OUT "$_\n";
    } else {
        # Replace 1/1 with 0/0 and 1|1 with 0|0
        s/1\/1/0\/0/g;
        s/1\|1/0\|0/g;
        print OUT "$_\n";
    }
}

close IN;
close OUT;

# Index the output gzipped VCF file
system("bcftools index -t $ARGV[1]");
