#!/usr/bin/perl

######################################################
# This script is for converting amino acid 
# three-letter code one-letter code for all
# amino acids in a given input file.
# @author Bian Li
# Report any bugs to Bian Li via bian.li@vanderbilt.edu
######################################################

use strict;
use warnings;

unless($#ARGV == 1)
{
  die "Usage: $0 <input_file> <output_file> ...";
}

# create streams for reading and writing
open(IN, "< $ARGV[0]") or die "Could not open file $ARGV[0] for reading: $!";
my @lines = <IN>;
close(IN);
open(OUT,"> $ARGV[1]") or die "Could not open file $ARGV[1] for writing: $!";

# amino acid conversion hash, only natural amino acids considered
my %aa_hash=(
  Ala=>'A',
  Arg=>'R',
  Asn=>'N',
  Asp=>'D',
  Cys=>'C',
  Glu=>'E',
  Gln=>'Q',
  Gly=>'G',
  His=>'H',
  Ile=>'I',
  Leu=>'L',
  Lys=>'K',
  Met=>'M',
  Phe=>'F',
  Pro=>'P',
  Ser=>'S',
  Thr=>'T',
  Trp=>'W',
  Tyr=>'Y',
  Val=>'V'
);

my @three_letter_codes = ("Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val");

# interate over all line in the input file and do the conversion
foreach my $line (@lines)
{
  chomp($line);
  foreach my $three_letter_code (@three_letter_codes)
  {
    $line =~ s/$three_letter_code/$aa_hash{$three_letter_code}/g;
  }
  print OUT "$line\n";
}

close(OUT)
