#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;
use Getopt::Long qw(GetOptions);
use Carp qw(carp croak);

my $dbpath ||= "/dors/meilerlab/apps/PISCES/opm/2014";
my $pdbid_file;

GetOptions(
	"dbpath=s" => \$dbpath,
	"pdbs=s"   => \$pdbid_file,
) or die "Failed to parse options: $!";

unless ( defined $pdbid_file ) {
	print
"Usage: $0 --pdbs <file with a list of 5-letter pdbs> --dbpath <path to database>\n";
	die "A file with a list of 5 letter pdb ids must be supplied";
}

open( my $fh, "<", $pdbid_file )
  or die "Could not open $pdbid_file for reading: $!";
chomp( my @pdb_ids = <$fh> );
close($fh) or die "Could not close $pdbid_file: $!";

foreach my $pdb_id (@pdb_ids) {
	my $midid     = substr( $pdb_id, 1, 2 );
	my $fourlt_id = substr( $pdb_id, 0, 4 );
	my $chain_id  = substr( $pdb_id, 4, 1 );
	my $pdbdir    = $dbpath . "/" . $midid;

	# cd to the directory
	chdir($pdbdir);
	my $protomer         = $pdb_id . ".pdb";
	my $oligomer         = $fourlt_id . ".bcl.pdb";
	my $naccess_protomer = $pdb_id . ".rsa";
	my $naccess_oligomer = $fourlt_id . ".rsa";

	# run naccess on protomer and oligomer
	if ( -s $naccess_protomer ) {
		print "$naccess_protomer already existed and is also valid\n";
	} else {
		system("naccess $protomer -p 1.4") == 0
                  or die "naccess failed to run for $protomer: $!";
	}
	if ( -s $naccess_oligomer ) {
		print "$naccess_oligomer already existed and is also valid\n";
	} else {
		system("naccess $oligomer -p 1.4") == 0
                  or die "naccess failed to run for $oligomer: $!";
	}

	# extract ASA and RSA from naccess output files and write the extracted
	# data to separate files
	unless ( -s $naccess_protomer && -s $naccess_oligomer ) {
		carp "$naccess_protomer or $naccess_oligomer or both are either not existing or empty";
	}

	# extracting
	my $protomer_extract = $pdb_id . ".proto_naccess";
	my $oligomer_extract = $pdb_id . ".oligo_naccess";
	if ( -s $protomer_extract && -s $oligomer_extract ) {
		print "Both $protomer_extract and $oligomer_extract exist and are valid, skip\n";
		next;
	}
	extract( $naccess_protomer, $protomer_extract, $chain_id );
	extract( $naccess_oligomer, $oligomer_extract, $chain_id );
}

sub extract {
	my ( $ifile, $ofile, $chain_id ) = @_;
	open $fh, "<", $ifile or die "Could not open $ifile for reading: $!";
	chomp( my @records = <$fh> );
	close($fh);
	open $fh, ">", $ofile or die "Could not open $ofile for writing: $!";
	foreach my $record (@records) {
		if ( $record !~ /^RES/ ) {
			next;
		}
		my @entries = split( /\s+/, $record );
		if ( $entries[2] ne $chain_id ) {
			next;
		}
		print $fh
		  sprintf( "%5d,%7.2f,%7.2f\n", $entries[3], $entries[4], $entries[5] );
	}
	close($fh);
}
