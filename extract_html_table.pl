#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long qw(GetOptions);
use HTML::TableExtract;
use LWP::Simple qw(get);
use Data::Dumper;

my $url;
my $headers;
my $output;

GetOptions(
	"url=s"     => \$url,
	"headers=s" => \$headers,
	"output=s"  => \$output,
) or die usage();

# check whether required parameters are passed in
unless ( defined $url && defined $headers && defined $output ) {
	die usage();
}

# print messages
print "The URL you passed in is: $url\n";

# split input headers
my @headers = split( /,/, $headers );
my @cleaned_headers = map { $_ =~ s/^\s+// } @headers;
print "The table headers you requested are: \n";
print Dumper( \@headers );

my $table_extract = HTML::TableExtract->new( headers => \@headers );
my $html = get($url);
$table_extract->parse($html);

# extract the first table
my ($table) = $table_extract->tables();
unless ( defined $table ) {
	die "No tables found in $url: $!";
}

# write extracted columns
open my $output_fh, ">", $output;
print $output_fh join( ",", @headers ), "\n";
foreach my $row ( $table->rows() ) {
	foreach ( @{$row} ) {
		s/,//g;
	}
	print $output_fh join( ",", @{$row} ), "\n";
}
close($output_fh);

#########################
# print out help info
#########################
sub usage {
	print
	  "Usage: $0 --url <URL> --headers <string of headers separated my commas> "
	  . "--output <output filename> \n";
}

# The following is POD documentation, use "perldoc extract_html_table.pl" from command line to read i

=head1 NAME

extract_html_table

=head1 SYNOPOSIS

extract_html_table.pl --url <URL> --headers <string of headers separated my commas>
--output <output filename>

=head1 DESCRIPTION

This script is intented for extracting information from a table in a HTML-based
web page. It assumes that the first table on the page is the one of interest.

=head1 AUTHOR

Bian Li, bian.li@vanderbilt.edu

You can use this script as well as modify it to suit your own need, or let me know
if you want me to modify it to suit your need. Also, you can sue me if you find any
annoying bugs. But before you do that please do give me a chance to fix the bug.
