#!/usr/bin/perl

##############################################################
#  Author     : Peng Hu
#  Email      ：penghu41@gmail.com
#  Last Edited: 08/26/2020
##############################################################

use strict;
use warnings;
use Getopt::Long;
my %opt;


GetOptions (\%opt,"read:s","tsv:s","qual:i","help");

my $help=<<USAGE;
This script could extract All mutation convertion tsv file which was generated by sam2tsv and then screen out the read name by defined quality score and gene annotation.

Usage: perl $0 -read read.txt -tsv sample.tsv -qual 27

USAGE

if ($opt{help} or keys %opt < 1) {
	print "$help\n";
	exit();
}

my $in = $opt{tsv};
my $qual = $opt{qual};

my $T_to_A = 0;
my $T_to_C = 0;
my $T_to_G = 0;
my $T_to_T = 0;
my $total_T = 0;

my $A_to_A = 0;
my $A_to_C = 0;
my $A_to_G = 0;
my $A_to_T = 0;
my $total_A = 0;

my $G_to_A = 0;
my $G_to_C = 0;
my $G_to_G = 0;
my $G_to_T = 0;
my $total_G = 0;

my $C_to_A = 0;
my $C_to_C = 0;
my $C_to_G = 0;
my $C_to_T = 0;
my $total_C = 0;


my $in2=$opt{read};
my %hash;
open (IN2,$in2) or die $!;
while (<IN2>) {
   chomp;
   $hash{$_} = "exist";
}
close IN2;

open (IN,$in) or die $!;
open (OUT,">${in}_q${qual}_gene_anno_stat.txt") or die $!;
#my $line= <IN>;
#print $line;
while (my$line=<IN>) {
	chomp $line;
	my @array = split "\t",$line;
	if ($hash{$array[0]}) {
        if (ord($array[6])-33 > $qual) { # modify the assay number for NOTE-seq use.
		    if ($array[1] & 16) { # negative strand
			    if ($array[8] =~ /A/){
				    $total_T++;
				    $T_to_T++ if ($array[5] =~ /A/);
                    $T_to_G++ if ($array[5] =~ /C/);
                    $T_to_C++ if ($array[5] =~ /G/);
                    $T_to_A++ if ($array[5] =~ /T/);
				} elsif ($array[8] =~ /G/){
				    $total_C++;
				    $C_to_T++ if ($array[5] =~ /A/);
                    $C_to_G++ if ($array[5] =~ /C/);
                    $C_to_C++ if ($array[5] =~ /G/);
                    $C_to_A++ if ($array[5] =~ /T/);
				} elsif ($array[8] =~ /C/){
				    $total_G++;
				    $G_to_T++ if ($array[5] =~ /A/);
                    $G_to_G++ if ($array[5] =~ /C/);
                    $G_to_C++ if ($array[5] =~ /G/);
                    $G_to_A++ if ($array[5] =~ /T/);
				} elsif ($array[8] =~ /T/){
				    $total_A++;
				    $A_to_T++ if ($array[5] =~ /A/);
                    $A_to_G++ if ($array[5] =~ /C/);
                    $A_to_C++ if ($array[5] =~ /G/);
                    $A_to_A++ if ($array[5] =~ /T/);
				} 
			} else { # positive strand
			    if ($array[8] =~ /A/){
				    $total_A++;
				    $A_to_T++ if ($array[5] =~ /T/);
                    $A_to_G++ if ($array[5] =~ /G/);
                    $A_to_C++ if ($array[5] =~ /C/);
                    $A_to_A++ if ($array[5] =~ /A/);
				} elsif ($array[8] =~ /G/){
				    $total_G++;
				    $G_to_T++ if ($array[5] =~ /T/);
                    $G_to_G++ if ($array[5] =~ /G/);
                    $G_to_C++ if ($array[5] =~ /C/);
                    $G_to_A++ if ($array[5] =~ /A/);
				} elsif ($array[8] =~ /C/){
				    $total_C++;
				    $C_to_T++ if ($array[5] =~ /T/);
                    $C_to_G++ if ($array[5] =~ /G/);
                    $C_to_C++ if ($array[5] =~ /C/);
                    $C_to_A++ if ($array[5] =~ /A/);
				} elsif ($array[8] =~ /T/){
				    $total_T++;
				    $T_to_T++ if ($array[5] =~ /T/);
                    $T_to_G++ if ($array[5] =~ /G/);
                    $T_to_C++ if ($array[5] =~ /C/);
                    $T_to_A++ if ($array[5] =~ /A/);
				}
			}
		}
	}
}
close IN;


print OUT "total_T\t${total_T}\n";
print OUT "T_to_A\t${T_to_A}\n";
print OUT "T_to_G\t${T_to_G}\n";
print OUT "T_to_C\t${T_to_C}\n";
print OUT "T_to_T\t${T_to_T}\n";

print OUT "total_A\t${total_A}\n";
print OUT "A_to_A\t${A_to_A}\n";
print OUT "A_to_G\t${A_to_G}\n";
print OUT "A_to_C\t${A_to_C}\n";
print OUT "A_to_T\t${A_to_T}\n";

print OUT "total_G\t${total_G}\n";
print OUT "G_to_A\t${G_to_A}\n";
print OUT "G_to_G\t${G_to_G}\n";
print OUT "G_to_C\t${G_to_C}\n";
print OUT "G_to_T\t${G_to_T}\n";

print OUT "total_C\t${total_C}\n";
print OUT "C_to_A\t${C_to_A}\n";
print OUT "C_to_G\t${C_to_G}\n";
print OUT "C_to_C\t${C_to_C}\n";
print OUT "C_to_T\t${C_to_T}\n";

close OUT;
