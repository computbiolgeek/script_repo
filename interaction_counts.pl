#!/usr/bin/perl

###################################################################
# @brief: this script is to be used for enumerating frequency of
# occurrence of interactions in a docked ensemble. It takes
# a list containing the pdb files of the ensemble, sequence
# ids of residues in the docking partners, and outputs a heatmap
# @version: 1.1
# @author: Bian Li
# @last modified: 02/29/2016
# @contact: bian.li@vanderbilt.edu
###################################################################

use strict;
use warnings;
use Getopt::Long;
use List::Util qw(sum);
use Data::Dumper;

my $pdb_list_file;
my $a_id_file;
my $b_id_file;
my $outfile;

# get command line arguments
GetOptions (
	"list=s" => \$pdb_list_file,
	"a_seq_ids=s" => \$a_id_file,
	"b_seq_ids=s" => \$b_id_file,
	"output=s" => \$outfile,
) or die "Could not parse command line options: $!";

# read in sequence ids of interest
open my $fh_a, "<", $a_id_file or die "Could not open $a_id_file for reading: $!";
chomp(my @a_seq_ids = <$fh_a>);
close $fh_a;
open my $fh_b, "<", $b_id_file or die "Could not open $b_id_file for reading: $!";
chomp(my @b_seq_ids = <$fh_b>);
close $fh_b;

# read in paths to pdb files
open my $fh_pdb, "<", $pdb_list_file or die "Could not open $pdb_list_file for reading: $!";
chomp(my @pdbs = <$fh_pdb>);
close $fh_pdb;

my $interaction_counts = "";

# check whether each residue pair is interacting
my $a_chainid = shift @a_seq_ids;
my $b_chainid = shift @b_seq_ids;
foreach my $name_id_pair_a (@a_seq_ids) 
{
	my ($resname_a, $seq_id_a) = split(/\s+/, $name_id_pair_a);
	$interaction_counts .= "$resname_a"."$seq_id_a";
	foreach my $name_id_pair_b (@b_seq_ids) 
	{
		my ($resname_b, $seq_id_b) = split(/\s+/, $name_id_pair_b);
		my $count = 0;
		# iterate over all complexes
		foreach my $pdb (@pdbs) 
		{
			# read in the current pdb
			open my $fh_pdb, "<", $pdb or die "Could not open $pdb for reading: $!";
			chomp(my @pdb_lines = <$fh_pdb>);
			# find the residues
			my @residue_records_a = get_residue_records($a_chainid, $seq_id_a, @pdb_lines);
			my @residue_records_b = get_residue_records($b_chainid, $seq_id_b, @pdb_lines);
			# check whether they are interacting
			if(is_interacting(\@residue_records_a, \@residue_records_b))
			{
				$count++;
			}
			
		}
		$interaction_counts .= ",$count";	
	}
	$interaction_counts .= "\n";
}

# write results in to file
open my $fh_out, ">", $outfile;
my $header = "Residues";
foreach my $name_id_pair_b (@b_seq_ids)
{
	$name_id_pair_b =~ s/\s+//;
	$header .= ",$name_id_pair_b"
}
$header .= "\n";
print $fh_out $header;
print $fh_out $interaction_counts;
close $fh_out;


############################################################
# Determine whether two residues are interacting. They
# are interacting if the separation between their side
# chain centroids is less than 10 angstroms and they are
# pointing toward each other.
# Parameters: pdbs records of the residues in consideration
# Return: a Boolean indicating whether they are interacting
############################################################
sub is_interacting
{
	my ($resA, $resB) = @_;
	my @resA = @{$resA};
	my @resB = @{$resB};
	my @centroidA = get_sidechain_centroid(@resA);
	my @centroidB = get_sidechain_centroid(@resB);
	my $centroid_separation = distance(\@centroidA, \@centroidB);
	my @caA = ca_coordinates(@resA);
	my @caB = ca_coordinates(@resB);
	my $ca_separation = distance(\@caA, \@caB);
	return $centroid_separation <= 10 && $ca_separation > $centroid_separation;	
}

#######################################################
# Extracts the pdb records of the residue of interest
# Parameters: chain id, sequence id, all records of the
#             pdb file
# Return: the pdb records of the residue of interest
#######################################################
sub get_residue_records
{
	my ($chainid, $seqid, @pdb_lines) = @_;
	my @residue_records;
	foreach my $line (@pdb_lines) 
	{
		# skip non-atom records
		if($line !~ /^ATOM/)
		{
			next;
		}
		if
		(
			substr($line, 21, 1) eq $chainid 
			&& substr($line, 22, 4) == $seqid 
			&& substr($line, 12, 4) !~ /H/ # remove hydrogen atoms
		) 
		{
			push(@residue_records, $line);
		}
	}
	return @residue_records;
}

#################################################
# Extract the three-letter code from pdb record
#################################################
sub get_aa_3lcode
{
	my @residue_records = @_;
	return substr($residue_records[0], 17, 3);
}


##########################################################
# Given the three-letter code, return the one-letter code
##########################################################
sub get_aa_1lcode
{
	# amino acid conversion hash, only natural amino acids considered
	my %aa_hash=(
	  Ala => 'A',
	  Arg => 'R',
	  Asn => 'N',
	  Asp => 'D',
	  Cys => 'C',
	  Glu => 'E',
	  Gln => 'Q',
	  Gly => 'G',
	  His => 'H',
	  Ile => 'I',
	  Leu => 'L',
	  Lys => 'K',
	  Met => 'M',
	  Phe => 'F',
	  Pro => 'P',
	  Ser => 'S',
	  Thr => 'T',
	  Trp => 'W',
	  Tyr => 'Y',
	  Val => 'V',
	);
	my $three_letter_code = ucfirst @_;
	return $aa_hash{$three_letter_code};
}

#####################################################
# Extract the sequence id of the residue of interest
#####################################################
sub get_residue_seqid
{
	my @residue_records = shift @_;
	return substr($residue_records[0], 22, 4);
}


###########################################################
# Computes the centroid of side chain atoms
# Parameters: pdb records for the residue of interest
#             stored in an array
# Return: an array containing the Cartesian coordiantes
#         of the coordinates of the side chain centroid
###########################################################
sub get_sidechain_centroid
{
	my @residue_records = @_;
	my ($x, $y, $z) = (0, 0, 0);
	my $num_atoms = 0;
	foreach my $record (@residue_records)
	{
		my $atom = substr($record, 12, 4);
		if
		(
			$atom eq " N  "
			|| $atom eq " CA "
			|| $atom eq " C  "
			|| $atom eq " O  "
		)
		{
			next;
		}
		$x += substr($record, 30, 8);
		$y += substr($record, 38, 8);
		$z += substr($record, 46, 8);
		$num_atoms++;
	}
	my @centroid = ($x / $num_atoms, $y / $num_atoms, $z / $num_atoms);
	return @centroid;
}

######################################################
# Computes the Euclidean distance between two points
# Parameters: array references to the coordinates of
#             the two points
# Return: the Euclidean distance
######################################################
sub distance
{
	my ($pointA, $pointB) = @_;
	my @pointA = @{$pointA};
	my @pointB = @{$pointB};
	return sqrt(sum(map {($pointA[$_] - $pointB[$_])**2} (0 .. 2)));
}

#######################################################
# Extracts the coordinates of the CA atom of the
# the residue of interest
# Parameter: pdb records for the residue of interest
#             stored in an array
# Return: the coordinates of the CA atom of the
#         the residue of interest
#######################################################
sub ca_coordinates
{
	my @residue_records = @_;
	my @ca_coordinates;
	foreach my $record (@residue_records)
	{
		my $atom = substr($record, 12, 4);
		if($atom eq " CA ")
		{
			@ca_coordinates = (substr($record, 30, 8), 
							   substr($record, 38, 8), 
							   substr($record, 46, 8));
			last;
		}
	}
	return @ca_coordinates;
}

######################################################
# Compute the vector from CA to side chain centroid
# Parameter: pdb records of the residue of interest
# Return: the vector stored in an array
######################################################
sub ca_centroid_vector
{
	my @residue_records = @_;
	my @sidechain_centroid = get_sidechain_centroid(@residue_records);
	my @ca_coordinates = ca_coordinates(@residue_records);
	my @ca_centroid_vec = map {$sidechain_centroid[$_] - $ca_coordinates[$_]} (0 .. 2);
	return @ca_centroid_vec;
}

############################################
# Compute the dot product of two vectors
# Parameters: two vectors represented as arrays
# Return: the dot product
############################################
sub dot_product
{
	my ($vecA, $vecB) = @_;
	my @vecA = @{$vecA};
	my @vecB = @{$vecB};
	return sum(map {$vecA[$_] * $vecB[$_]} (0 .. 2));
}