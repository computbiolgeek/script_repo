#!/usr/bin/perl

###############################################################################
# This script batch renames files with a <old> pattern in the filename to files
# with <old> replaced by <new> in the filename.
# @author Bian Li
# Please report bugs to Bian Li via bian.li@vanderbilt.edu
###############################################################################

use strict;
use warnings;

# make sure there are two command line arguments
unless($#ARGV == 2) {
	die "Usage: $0 <old> <new> <pattern between quotes>";
}

# rename every file with <old> in the filename to file with <old> replaced by <new>
foreach my $file (glob "$ARGV[2]") {
	my $new_file = $file;
	$new_file =~ s/$ARGV[0]/$ARGV[1]/;
	if(-e $new_file) {
		warn "can't rename $file to $new_file: $new_file already exists $!\n";
	}
	else {
		print "renaming $file to $new_file\n";
		rename $file, $new_file or warn "rename $file to $new_file failed $!\n";	
	}
}
