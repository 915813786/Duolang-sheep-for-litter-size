#!/usr/bin/perl
use strict;
die "Usage:$0 <input list><out reads status>\n" unless @ARGV==2;
open (IN,"$ARGV[0]") or die "input file requirement!\n";
open (OUT,">$ARGV[1]") or die "out file requirement\n";

while(<IN>){
	chomp;
	my @F=split/\t/,$_;
	my $raw_base = 0;
	my $clean_base = 0;
	open (IN1,"$F[0]") or die "$F[0] cannot open!\n";
	my $i = 0;
	my $id = $F[0];
	$id =~ s/\/(.*)\///;
	$id =~ s/.stat//;
	while(<IN1>){
		chomp;
		$i++;
		s/,//g;
		my @fastq = split/\s+/,$_;
		if ($i == 2){
			$raw_base += $fastq[4];
		}elsif ($i == 4){
			$raw_base += $fastq[4];
		}elsif ($i == 6){
			$clean_base += $fastq[4];
		}elsif ($i == 8){
			$clean_base += $fastq[4];
		}
	}
	open (IN2,"$F[1]") or die "$F[1] cannot open!\n";
	my $i = 0;
	my $raw_gc;
	my $raw_q20;
	my $raw_q30;
	my $raw_seq;
	while(<IN2>){
		chomp;
		my @fastq = split/\s+/,$_;
		if (/^Total Sequences/){
			$raw_seq = $fastq[2];
		}
		if (/^%GC/){
			$raw_gc = $fastq[1];
		}
		if(/^#Quality/){
			$i++;
			next;
		}
		if ($i == 1){
			if ($fastq[0] > 20){
				$raw_q20 += $fastq[1];
			}
			if ($fastq[0] > 30){
				$raw_q30 += $fastq[1];
			}
		}
		if (/END_MODULE/){
			$i = 0;
		}
	}
        open (IN3,"$F[2]") or die "$F[2] cannot open!\n";
        my $i = 0;
        while(<IN3>){
                chomp;
                my @fastq = split/\s+/,$_;
                if (/^Total Sequences/){
                        $raw_seq += $fastq[2];
		}
                if (/^%GC/){
                        $raw_gc += $fastq[1];
                }
                if(/^#Quality/){
                        $i++;
                        next;
                }
                if ($i == 1){
                        if ($fastq[0] > 20){
                                $raw_q20 += $fastq[1];
                        }
                        if ($fastq[0] > 30){
                                $raw_q30 += $fastq[1];
                        }
                }
		if (/END_MODULE/){
			$i = 0;
		}
        }
	$raw_gc = $raw_gc/2;
	$raw_q20 = $raw_q20/$raw_seq;
	$raw_q30 = $raw_q30/$raw_seq;
	$raw_seq = $raw_seq/2;
	open (IN4,"$F[3]") or die "$F[3] cannot open!\n";
        my $i = 0;
        my $clean_gc;
        my $clean_q20;
        my $clean_q30;
        my $clean_seq;
	while(<IN4>){
		chomp;
		my @fastq = split/\s+/,$_;
                if (/^Total Sequences/){
                        $clean_seq = $fastq[2];
		}
                if (/^%GC/){
                        $clean_gc = $fastq[1];
                }
                if(/^#Quality/){
                        $i++;
                        next;
                }
                if ($i == 1){
                        if ($fastq[0] > 20){
                                $clean_q20 += $fastq[1];
                        }
                        if ($fastq[0] > 30){
                                $clean_q30 += $fastq[1];
                        }
                }
                if (/END_MODULE/){
                        $i = 0;
                }
        }
        open (IN5,"$F[4]") or die "$F[4] cannot open!\n";
        my $i = 0;
        while(<IN5>){
                chomp;
                my @fastq = split/\s+/,$_;
                if (/^Total Sequences/){
                        $clean_seq += $fastq[2];
		}
                if (/^%GC/){
                        $clean_gc += $fastq[1];
                }
                if(/^#Quality/){
                        $i++;
                        next;
                }
                if ($i == 1){
                        if ($fastq[0] > 20){
                                $clean_q20 += $fastq[1];
                        }
                        if ($fastq[0] > 30){
                                $clean_q30 += $fastq[1];
                        }
                }
                if (/END_MODULE/){
                        $i = 0;
                }
        }
	$clean_gc = $clean_gc/2;
        $clean_q20 = $clean_q20/$clean_seq;
        $clean_q30 = $clean_q30/$clean_seq;
	$clean_seq = $clean_seq/2;
	print OUT "$id\t$raw_base\t$clean_base\t$raw_seq\t$clean_seq\t$raw_gc\t$clean_gc\t$raw_q20\t$raw_q30\t$clean_q20\t$clean_q30\n";
}
close IN;
close IN1;
close IN2;
close IN3;
close IN4;
close IN5;
close OUT;
