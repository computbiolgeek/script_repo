#!/usr/bin/perl

$opm_pdb_dir="/dors/meilerlab/apps/PISCES/opm/2014/all/";
$local_pdb_dir="/dors/meilerlab/apps/PISCES/opm/2014/";


if($#ARGV != 0) {
  die "usage: $0 <pdb_list> ..."
}

my @pdb_list;
if( -f $ARGV[0]) {
  open(IN, "<$ARGV[0]") or die "could not open file $ARGV[0]";
  # push the entries into an array
  @pdb_list=<IN>;
  close IN;
} else {
  # get pdb ids from command line
  @file=@ARGV;
}

my $multimer_info="./oligomeric_state_bcl_pdb.txt";
open(OUT,">$multimer_info");
print OUT "PDB_ID\tTotal_Chains\tTM_Chains\n";

# loop over all entries
foreach my $entry (@pdb_list) {
  # check if the current entry is 4 letter code for 5 letter code
  chomp($entry);
  if(!(length($entry)==4 || length($entry)==5)) {
    print("The entry $entry should be either 4-letter pdb id or 5-letter chain code. Skip ...\n");
    next;
  }

  # get the 4-letter pdb id
  my $pdb_id=lc(substr($entry,0,4));
  my $mid_id=lc(substr($entry,1,2));

  my $pdb_filename=$opm_pdb_dir.$pdb_id.".pdb";
  my $local_pdb_filename=$local_pdb_dir.$mid_id."/".$pdb_id.".bcl.pdb";

  if(! -f $pdb_filename) {
    print "$pdb_filename does not exist, now trying to download it from the opm web server ...";
    system("curl -m 60 -s http://opm.phar.umich.edu/pdb/".$pdb_id.".pdb -o ".$pdb_filename);
    open(IN, "<$pdb_filename");
    my $first_line=<IN>;
    close IN;
    if( $first_line=~/^REMARK/) {
      print "Download succeeded ...";
    } elsif($first_line=~/^html/) {
      print "Download failed, go to http://opm.phar.umich.edu/pdb/ to manually download $pdb_id ... ";
      unlink($pdb_filename);
      next;
    }
  }

  # get oligomeric information
  my $get_chains_opm = "grep '^ATOM ' $local_pdb_filename | awk '{if( ".'$3 == "N"){print substr($0,22,1)}}\' | sort | uniq | tr -d \'\\n\'';
  my $chains_opm = qx( $get_chains_opm);
  my $number_chains=length($chains_opm);
  my $get_tm_chains="awk '/^ATOM/{if("."\$9 < 1 && \$9 > -1){print \$5}}' ".$local_pdb_filename." | uniq | tr -d '\n'";
  my $tm_chains=qx($get_tm_chains);
  my $number_tm_chains=length($tm_chains);
  if($number_tm_chains>$number_chains) {
    $tm_chains=$chains_opm;
    $number_tm_chains=$number_chains;
  }
  print "$pdb_id has chains $chains_opm and ";
  print "$tm_chains are transmembrane domains.\n";
  print OUT "$pdb_id\t$number_chains\t$number_tm_chains\n";
}

close OUT;
exit
