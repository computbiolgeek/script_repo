#!/usr/bin/perl

###########################################################
# Warning: please do not make any attempt to try to
# understand this code. It's Perl and it also has a lot
# of hackings!
###########################################################

use strict;
use warnings;
use Carp qw(carp croak);
use Scalar::Util qw(looks_like_number);
use List::Util qw(sum);
use List::MoreUtils qw(uniq);
use Getopt::Long qw(GetOptions);
use Data::Dumper;
use LWP::UserAgent qw(get);

my $path2db = "/dors/meilerlab/apps/PISCES/opm/2014";
my $pdb_list;
my $outfile;

GetOptions(
	"list=s" => \$pdb_list,
	"out=s"  => \$outfile,
) or die usage();

open my $fh_in, "<", $pdb_list;
chomp(my @pdbids = <$fh_in>);
close($fh_in);

# connects the output file to the file handle for writing
open my $fh_out, ">", $outfile;
# print the header
my $header = sprintf("%-8s%8s%14s%14s%14s%14s%14s%14s\n", 
	"pdbid", "total", "tmh_res", "helices", "tm_helices", "bent_helices", "subunits", "family_id"
);
print $fh_out $header;

# iterate over all pdbids
foreach my $pdbid (@pdbids) {
	print "Processing $pdbid\n";
	my $four_ltr_id = substr($pdbid,0,4);
	my $mid_id = substr($pdbid, 1, 2);
	my $chainid = substr($pdbid, -1, 1);
	# parse the pdb file for interesting records
	my $pdb_file_opm = $path2db."/".$mid_id."/".$four_ltr_id.".pdb"; # this pdb file is only used for getting half-thickness
	my $pdb_file_bcl = $path2db."/".$mid_id."/".$four_ltr_id.".bcl.pdb";
	my %pdb_hash_chain = parse_pdb($pdb_file_bcl, $chainid);
	my %pdb_hash_all = parse_pdb($pdb_file_bcl);
	# parse the MAHSSMI file
	my $mahssmi_file = $path2db."/".$mid_id."/".$four_ltr_id.".mahssmi";
	my @mahssmi_lines = parse_mahssmi($mahssmi_file, $chainid);
	# extract membrane half-thickness from the first line of the pdb file
	open my $opm_fh, "<", $pdb_file_opm or die "Could not open $pdb_file_opm for reading: $!";
	my $first_line = <$opm_fh>;
	close ($opm_fh);
	chomp ($first_line);
	my @fl_tokens = split (/\s+/, $first_line);
	my $half_thickness = $fl_tokens[-1];
	# warn the user when the half-thickness is not a number
	unless(looks_like_number($half_thickness)) {
		carp "$pdb_file_opm does not seem to contain a valid half-thickness, please check";
		next;
	}
	print "The membrane half_thickness for $four_ltr_id is $half_thickness\n";
	# get the total number of residues
	my $total = get_size($pdbid, $pdb_hash_chain{SEQRES});
	# get the number of helices
	my @helices = get_helices_mahssmi(\@mahssmi_lines);
	my $helices = scalar @helices;
	print "Found the following helices in the MAHSSMI file for $four_ltr_id\n";
	print Data::Dumper->Dump([\@helices], [qw(*Helices)]);
	# get the number of transmembrane helical residues
	my $tmh_res = get_number_tmh_residues($pdbid, $half_thickness, \@helices, $pdb_hash_chain{ATOM});
	# get the number of transmemrbane helices
	my $tm_helices = get_number_tmh($pdbid, $half_thickness, \@helices, $pdb_hash_chain{ATOM});
	# get the number of bent helices
	my $bent_helices = get_number_bent_helix($pdbid, \@helices, $pdb_hash_chain{ATOM});
	# get the number of subunits
	my $subunits = get_number_subunits($four_ltr_id, $pdb_hash_all{SEQRES});
	# get the superfamily id
	my $family_id = get_family_id($four_ltr_id);
	my $out_line = sprintf("%-8s%8d%14d%14d%14d%14d%14d%14d\n", 
		$pdbid, $total, $tmh_res, $helices, $tm_helices, $bent_helices, $subunits, $family_id);
	print $fh_out $out_line;
}

close $fh_out;

########################################
# Tells the user how to use this script
########################################
sub usage {
	my $help_string = "";
	$help_string .= "./mp_summary.pl --list <pdb_list> --out output_filename\n";
	print $help_string;
}

sub parse_pdb {
	my ($pdb_file, $chainid) = @_;
	open my $fh_pdb, "<", $pdb_file;
	chomp(my @pdb_lines = <$fh_pdb>);
	my %pdb_hash = (
		FIRSTLINE => $pdb_lines[0],
		SEQRES    => [],
		HELIX     => [],
		ATOM      => [],
	);
	foreach my $line (@pdb_lines) {
		if(defined $chainid) {
			if($line =~ /^SEQRES/ && substr($line, 11, 1) eq $chainid) {
				push(@{$pdb_hash{SEQRES}}, $line);
			}
			if($line =~ /^HELIX/ && substr($line, 19, 1) eq $chainid) {
				push(@{$pdb_hash{HELIX}}, $line);
			}
			if($line =~ /^ATOM/ && substr($line, 21, 1) eq $chainid) {
				push(@{$pdb_hash{ATOM}}, $line);
			}
		}
		else {
			if($line =~ /^SEQRES/) {
				push(@{$pdb_hash{SEQRES}}, $line);
			}
			if($line =~ /^HELIX/) {
				push(@{$pdb_hash{HELIX}}, $line);
			}
			if($line =~ /^ATOM/) {
				push(@{$pdb_hash{ATOM}}, $line);
			}			
		}
	}
	return %pdb_hash;
}

#######################################################################################
# Parse the MAHSSMI file to extract the records for the interested chain into an array
#######################################################################################
sub parse_mahssmi {
	my ($mahssmi_file, $chainid) = @_;
	print "Parsing mahssmi file $mahssmi_file\n";
	open my $fh_mahssmi, "<", $mahssmi_file 
		or die "Could not open $mahssmi_file for reading, make it exists: $!";
	chomp(my @mahssmi_lines = <$fh_mahssmi>);
	my @chain_lines;
	for my $line (@mahssmi_lines) {
		if($line =~ /^#/) { # skip comments
			next;
		}
		my @line_tokens = split(/\s+/, $line);
		if($line_tokens[2] ne $chainid) {
			next;
		}
		push(@chain_lines, $line);
	}
	return @chain_lines;
}

sub get_size {
	my ($pdbid, $seqres_records) = @_;
	unless(@{$seqres_records}) {
		croak "No SEQRES records found for $pdbid";
	}
	my @seqres_elements = split(/\s+/, $seqres_records->[0]);
	return $seqres_elements[3];
}

######################################################################################
# Count the number of HELIX records for a particular chain. Note that this does
# not necessarily equal to the actual number of helices. For example, a helix broken 
# by a proline might split up into two records in some pdb files
######################################################################################
sub get_number_helices {
	my ($pdbid, $helix_records) = @_;
	unless(@{$helix_records}) {
		carp "No HELIX records found for $pdbid";
		return 0;
	}
	return scalar @{$helix_records};
}

######################################################################################
# Count the number of helices based on the more accurate MAHSSMI file
######################################################################################
sub get_helices_mahssmi {
	my $mahssmi_lines = shift @_;
	my @helices;
	my $first_encounter = 1;
	my $helix_start;
	my $helix_end;
	my $inserted = 0;
	my $first_line = 1;
	foreach my $line (@{$mahssmi_lines}) {
		my @line_tokens = split(/\s+/, $line);
		if ($line_tokens[3] ne "H" && $line_tokens[4] ne "M" && !$first_encounter) {
			if ($inserted) {
				next;
			}
			$helix_end = $line_tokens[0] - 1;
			push(@helices, $helix_start." ".$helix_end);
			$inserted = 1;
			$first_encounter = 1;
		}
		else {
			if (!$first_encounter) {
				next;
			}
			if ($line_tokens[3] eq "H") {
				$helix_start = int($line_tokens[0]);
				$inserted = 0;
				$first_encounter = 0;
			}
		}
		$first_line = 0;
	}
	return @helices;
}

######################################################################################
# Count the number of transmembrane helices based on the more accurate MAHSSMI file
######################################################################################
sub get_number_tmh_mahssmi {
	my $mahssmi_lines = shift @_;
	my $count = 0;
	my $first_encounter = 1;
	foreach my $line (@{$mahssmi_lines}) {
		my @line_tokens = split(/\s+/, $line);
		if ($line_tokens[4] ne "M") {
			$first_encounter = 1;
		}
		else {
			if (!$first_encounter) {
				next;
			}
			$count++;
			$first_encounter = 0;
		}		
	}
	return $count;
}

sub get_number_tmh_residues {
	my ($pdbid, $half_thickness, $helices, $atom_records) = @_;
	unless (@{$helices}) {
		carp "No HELIX records found for $pdbid";
		return 0;
	}
	my $count = 0;
	foreach my $helix (@{$helices}) {
		if (!is_tmh_pdb($half_thickness, $helix, $atom_records)) {
			next;
		}
		my ($start, $end) = split(/\s+/, $helix);
		$count += $end - $start + 1;
	}
	return $count;
}

sub get_number_tmh {
	my ($pdbid, $half_thickness, $helices, $atom_records) = @_;
	my $count = 0;
	foreach my $helix (@{$helices}) {
		if (is_tmh_pdb($half_thickness, $helix, $atom_records)) {
			$count++;
		}
	}
	unless ($count) {
		carp "No transmembrane helices found in $pdbid, please verify manually";
	}
	return $count;
}

########################################################################################
# A helix is considered a transmembrane helix if its center locates inside the memrbane
########################################################################################
sub is_tmh_pdb {
	my ($half_thickness, $helix, $atom_records) = @_;
	my ($start, $end) = split(/\s+/, $helix);
	my @ca_start = get_ca_coordinates($start, $atom_records);
	my @ca_end = get_ca_coordinates($end, $atom_records);
	my $middle_z = ($ca_start[2] + $ca_end[2]) / 2;
	if ($middle_z <= $half_thickness && $middle_z >= -$half_thickness && ($end - $start) >= 12)
	{
		return 1;
	}
	return 0;
}

sub get_ca_coordinates {
	my ($seqid, $atom_records) = @_;
	my @ca_coordinates;
	my @residue_records = get_residue_records($seqid, $atom_records);
	unless (@residue_records) {
		croak "Missing coordinates for residue $seqid in the pdb file, please check";
	}
	foreach my $record (@residue_records) {
		my $atom = substr($record, 12, 4);
		if ($atom eq " CA ") {
			@ca_coordinates = (substr($record, 30, 8), 
							   substr($record, 38, 8), 
							   substr($record, 46, 8));
			last;
		}
	}
	return @ca_coordinates;
}

sub get_residue_records {
	my ($seqid, $atom_records) = @_;
	my @residue_records;
	foreach my $line (@{$atom_records}) {
		if (
			substr($line, 22, 4) == $seqid 
			&& substr($line, 12, 4) !~ /H/ # remove hydrogen atoms
		) {
			push(@residue_records, $line);
		}
	}
	return @residue_records;
}

sub get_cb_coordinates {
	my $residue_records = shift @_;
	my @cb_coordinates;
	foreach my $record (@{$residue_records}) {
		my $atom = substr($record, 12, 4);
		if ($atom eq " CB ") {
			@cb_coordinates = (substr($record, 30, 8), 
							   substr($record, 38, 8), 
							   substr($record, 46, 8));
			last;
		}
	}
	return @cb_coordinates;	
}

sub get_sidechain_centroid(\@) {
	my $residue_records = shift @_;
	my ($x, $y, $z) = (0, 0, 0);
	my $num_atoms = 0;
	foreach my $record (@{$residue_records}) {
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

###################################################################
#
###################################################################
sub get_number_subunits {
	my ($pdbid, $seqres_records) = @_;
	unless(@{$seqres_records}) {
		croak "No SEQRES records found for $pdbid";
	}
	my @chain_ids = ();
	foreach my $record (@{$seqres_records}) {
		my @seqres_elements = split(/\s+/, $record);
		push(@chain_ids, $seqres_elements[2]);
	}
	my @uniq_chain_ids = uniq(@chain_ids);
	return scalar @uniq_chain_ids;
}

sub distance {
	my ($pointA, $pointB) = @_;
	return sqrt (sum(map {($pointA->[$_] - $pointB->[$_])**2} (0 .. 2)));
}

#########################################################################
# First compute the length of the line connecting the two end points
# of the helix. Then compute the actual length of the helix:
# rise_per_residue * number_of_residues. If actual_length - ideal_length
# is greater than a certain threshold, then the helix is bent.
#########################################################################
sub get_number_bent_helix {
	my ($pdbid, $helices, $atom_records) = @_;
	my $count = 0;
	foreach my $helix (@{$helices}) {
		my ($start, $end) = split (/\s+/, $helix);
		my @ca_start = get_ca_coordinates($start, $atom_records);
		my @ca_end = get_ca_coordinates($end, $atom_records);
		my $ideal_length = distance(\@ca_start, \@ca_end);
		my $rise_per_residue = 1.5;
		my $actual_length = $rise_per_residue * ($end - $start + 1);
		if ($ideal_length == 0) {
			next;
		}
		my $is_bent_helix = ($actual_length - $ideal_length) / $ideal_length >= 0.1 ? 1 : 0;
		if ($is_bent_helix) {
			$count++;
		}
	}
	unless ($count) {
		carp "No bent helices found in $pdbid";
	}
	return $count;
}

####################################################################
# Fetch the html page of the given pdb, then extract the family id.
# If the html page of the given pdb does not contain the family id,
# fetch the html page of the representative pdb, then extract the
# family id.
####################################################################
sub get_family_id {
	my $pdbid = shift @_;
	my $url = "http://opm.phar.umich.edu/protein.php?search=".$pdbid;
	my $user_agent = LWP::UserAgent->new();
	$user_agent->timeout(30);
	my $response = $user_agent->get($url);
	unless ($response->is_success) {
		die $response->status_line;
	}
	my $html = $response->decoded_content;
	my $family_id;
	if ($html =~ /families.php\?family=([0-9]+)/) {
		$family_id = $1;
	}
	elsif ($html =~ /protein.php\?pdbid=([0-9a-z]+)/) { # get representative pdb
		$url = "http://opm.phar.umich.edu/protein.php?search=".$1;
		$response = $user_agent->get($url);
		$html = $response->decoded_content;
		if ($html =~ /families.php\?family=([0-9]+)/) {
			$family_id = $1;
		}
	}
	return $family_id;
}
