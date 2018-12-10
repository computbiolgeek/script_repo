#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my $pdb_filename;
my $old_id          ||= "A";
my $new_id          ||= "B";
my $output_filename ||= "output.pdb";

# process command line options
GetOptions(
	"pdb|p=s"        => \$pdb_filename,
	"current_id|c=s" => \$old_id,
	"new_id|n=s"     => \$new_id,
	"output|o=s"     => \$output_filename,
) or die "Failed to process command line options: $!";

# check command line options
unless ( defined $pdb_filename ) {
	die "Usage: $0 -p <input.pdb> -c <chain_id> "
	  . "-n <chain_id> -o <output.pdb>" . "\n";
}

# open pdb file
open( my $FH, $pdb_filename )
  or die "Cannot open pdb file $pdb_filename for reading: $!";
my @pdb_records = <$FH>;
close($FH);

# do the modification and write to a new pdb file
open( $FH, ">", $output_filename )
  or die "Cannot open file $output_filename for writting: $!";
foreach my $record (@pdb_records) {

	# modify the chain id
	if ( $record =~ /^ATOM/ && substr( $record, 21, 1 ) eq $old_id ) {
		substr( $record, 21, 1 ) = $new_id;
	}

	# remove anything in fields from 67 to 77
	if ( $record =~ /^ATOM/ ) {
		substr( $record, 66, 10 ) =~ s/[A-Z0-9]+/ /g;
	}
	print $FH $record;
}

exit;

=head1 NAME

rename_pdb_chain_id.pl

=head1 SYNOPSIS

rename_pdb_chain_id.pl -p <input.pdb> -c <chain_id> -n <chain_id> -o <output.pdb>

=head1 DESCRIPTION

This script was written for renaming chain IDs in a PDB file.

=head1 AUTHOR

Bian Li, bian.li@vanderbilt.edu

You can use this script as well as modify it to suit your own need, or let me know if you want me to modify it to suit your need. 
Also, you can sue me if you find any annoying bugs. But before you do that please do give me a chance to fix the bug.
