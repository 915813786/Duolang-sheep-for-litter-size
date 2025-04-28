#!/usr/bin/perl
use strict;
use warnings;
die "Usage:$0 <input bam.stat> <input vcf.gz> <filtered vcf.gz>\n" unless @ARGV==3;

my $input_bamStat = $ARGV[0] or die "Usage: $0 bam.stat\n";
my $input_vcf = $ARGV[1] or die "Usage: $0 input.vcf.gz\n";
my $output_vcf = $ARGV[2] or die "Usage: $0 output.vcf.gz\n";

# Initialize counters
my ($total_variants, $filtered_genotypes, $output_variants) = (0, 0, 0);

# Open input (gzipped VCF) and output (temporary file)
open(my $depth_in, "$input_bamStat") or die "Cannot open input bam.stat: $!\n";
open(my $vcf_in, "bcftools view $input_vcf |") or die "Cannot open input VCF: $!\n";
open(my $vcf_out, ">", "filtered_samples.tmp.vcf") or die "Cannot write temporary file: $!\n";

my %hash;
while (<$depth_in>) {
    chomp;
    my @F = split/\t/,$_;
    my $depth = $F[6];
    $depth =~ s/X//;  # Remove 'X' from depth
    $hash{$F[0]} = $depth;
}

my %depth;
while (<$vcf_in>) {
    chomp;
    if (/^##/) {  # Keep header lines
        print $vcf_out "$_\n";
        next;
    }elsif(/^#C/){
        print $vcf_out "$_\n";
        my @id = split /\t/;
        foreach my $index (9..$#id){
            if (exists $hash{$id[$index]}) {
                $depth{$index} = $hash{$id[$index]};
            } else {
                print "Error: the $id[$index] sample has no depth input!\n";
            }
        }
        next;
    }

    $total_variants++;
    my @fields = split /\t/;
    my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format) = @fields[0..8];

    # Skip monomorphic sites (ALT=".") or reference-only sites
    if ($alt eq ".") {
        next;
    }

    # Get field indices for GQ/AD/DP
    my @fmt_fields = split /:/, $format;
    my ($gq_idx, $ad_idx, $dp_idx);
    for (my $i = 0; $i < @fmt_fields; $i++) {
        $gq_idx = $i if ($fmt_fields[$i] eq "GQ");
        $ad_idx = $i if ($fmt_fields[$i] eq "AD");
        $dp_idx = $i if ($fmt_fields[$i] eq "DP");
    }

    my $has_alt_allele = 0;  # Track if ALT allele exists after filtering

    # Process each sample's genotype
    for (my $col = 9; $col < @fields; $col++) {
        my @sample_data = split /:/, $fields[$col];
        my $gt = $sample_data[0];

        # Skip if already missing
        next if ($gt eq "./.");

        # Filter 1: GQ < 20
        if (defined $gq_idx && $sample_data[$gq_idx] ne "." && $sample_data[$gq_idx] < 20) {
            $fields[$col] = "./." . ( @fmt_fields > 1 ? ":" . join(":", @sample_data[1..$#sample_data]) : "" );
            $filtered_genotypes++;
            next;
        }

        # Filter 2: AD ratio < 0.2 (Handle different genotypes)
        if (defined $ad_idx && $sample_data[$ad_idx] ne ".") {
            my @ad = split /,/, $sample_data[$ad_idx];
            if (@ad == 2) {
                # Standard biallelic case (e.g., 0/1, 1/0)
                my $total_ad = $ad[0] + $ad[1];
                if ($total_ad > 0 && $ad[1]/$total_ad < 0.2 && ($gt =~ /0\/1|1\/0|0\|1|1\|0/)) {
                    $fields[$col] = "./." . ( @fmt_fields > 1 ? ":" . join(":", @sample_data[1..$#sample_data]) : "" );
                    $filtered_genotypes++;
                    next;
                }
            } else {
                # Handle cases with more than two alleles
                # You can add more complex logic here for multiallelic sites
            }
        }

        # Filter 3: DP < 1/3 depth or > 3 depth
        if (defined $dp_idx && $sample_data[$dp_idx] ne "." && ($sample_data[$dp_idx] < $depth{$col}/3 || $sample_data[$dp_idx] > $depth{$col}*3)) {
            $fields[$col] = "./." . ( @fmt_fields > 1 ? ":" . join(":", @sample_data[1..$#sample_data]) : "" );
            $filtered_genotypes++;
            next;
        }

        # Check for remaining ALT alleles
        $has_alt_allele = 1 if ($gt =~ /1/);
    }

    # Only keep variants with remaining ALT alleles
    if ($has_alt_allele) {
        print $vcf_out join("\t", @fields), "\n";
        $output_variants++;
    }
}

close $vcf_in;
close $vcf_out;

# Compress and index final output
system("bcftools view -Oz -o $output_vcf filtered_samples.tmp.vcf") == 0 
    or die "bcftools compression failed: $!\n";
system("bcftools index -t $output_vcf") == 0 
    or die "bcftools indexing failed: $!\n";

unlink("filtered_samples.tmp.vcf");

# Print summary statistics
print "\n=== FILTERING SUMMARY ===\n";
printf "Total variants processed:    %d\n", $total_variants;
printf "Genotypes filtered out:      %d\n", $filtered_genotypes;
printf "Biallelic variants retained: %d\n", $output_variants;
print "Output file: $output_vcf\n";
