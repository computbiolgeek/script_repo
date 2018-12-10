#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use List::MoreUtils qw(uniq);

unless($#ARGV==5)
{
  die &Usage();
}

sub Usage()
{
  print "./delete_command_lines\n".
	"--input(-i) <filename>\n".
	"--common(-c) <filename>\n".
	"--output(-o) <filename>\n";
}

my $input = "";
my $common = "";
my $output = "";
# my @outputs = ();

GetOptions
(
  "input|i=s" => \$input,
  "common|c=s" => \$common,
  "output|o=s" => \$output
) or die "Error in command line arguments\n";

# open input file
open(IN, "< $input") or die "Could not open $input for reading: $!";
my @inputs = <IN>;
close(IN);

# open file containing common lines
open(CM, "< $common") or die "Could not open $common for reading: $!";
my @commons = <CM>;
close(CM);

# open file for writing
open(OUT, "> $output") or die "Could not open $common for writing: $!";

# remove newline charater from each entry
foreach my $line (@inputs)
{
  chomp($line);
}

# remove common lines from inputs
foreach my $common_line (@commons)
{
  chomp($common_line);
  my $index = 0;
  foreach my $line (@inputs)
  {
    if(index($inputs[$index], $common_line) == -1)
    {
      $index++;
    }
  }
  splice(@inputs, $index, 1);
}

# my @uniq_lines = uniq @outputs;

foreach my $line (@inputs)
{
  print OUT "$line\n";
}

close(OUT);
