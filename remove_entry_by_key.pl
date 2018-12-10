#!/usr/bin/perl

###############################################################
# This script takes as input a list of records keyed by keys
# and a list of keys. It removes records that begin with keys.
# @author Bian Li
# Report any bugs to Bian Li via bian.li@vanderbilt.edu
###############################################################


use strict;
use warnings;
use Getopt::Long;

# make sure the number of command line arguments is correct
if ($#ARGV != 7)
{
	die help();
}

###########################
sub help
###########################
{
  print "./remove_entry_by_key.pl\n".
        "--input(-i) <filename>\n".
        "--output(-o) <filename>\n".
        "--keys(-k) <filename>\n".
        "--removed(-r) <filename>\n";
}

# get command line flags and options
my $input = "";
my $output = "";
my $removed = "";
my $keys = "";

GetOptions
(
  "input|i=s" => \$input,
  "output|o=s" => \$output,
  "removed|r=s" => \$removed,
  "keys|k=s" => \$keys
) or die "Error in command line arguments\n";

# open input file
open(IN, "< $input") or die "Could not open $input for reading: $!";
my @entries = <IN>;
close(IN);
open(KEYS, "< $keys") or die "Could not open $keys for reading: $!";
my @keys = <KEYS>;
close(KEYS);

# open file for writing
open(OUT, "> $output") or die "Could not open $output for writing: $!";
open(RM, "> $removed") or die "Could not open $removed for writing: $!";

# remove undesired entries
my @undesired_entries;
foreach my $key (@keys)
{
	chomp($key);
	print "Removing entries that begin with $key\n";
	my @matched_entries = grep { $_ =~ /$key/} @entries;
	push @undesired_entries, @matched_entries; 
	@entries = grep { $_ !~ /$key/} @entries;
}

# write output
foreach my $entry (@entries)
{
	print OUT "$entry";
}
close(OUT);
foreach my $entry (@undesired_entries)
{
	print RM "$entry";
}
close(RM);
