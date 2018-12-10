#!/usr/bin/env python3

"""
    @summary: This script takes in a list of PDB IDs and retrieves the PDB coordinate file for each ID.
    @author: Bian Li
    @contact: bian.li@vanderbilt.edu
    @change: 1/17/2017
"""

from argparse import ArgumentParser
import os, sys
from Bio.PDB import PDBList, PDBParser, PDBIO

def main():
    """
    """
    # parse command line options and arguments
    parser = ArgumentParser( description = "Hi there, I download PDB file for you in batch!" )
    parser.add_argument( "-l", "--list", dest = "list",
                         help = "a file containing a list of PDB IDs one per line" )
    parser.add_argument( "-c", "--chains", dest = "chains", action = "store_true",
                         help = "whether the PDB IDs contains chain IDs" )
    parser.add_argument( "-d", "--dir", dest = "dir",
                         help = "path to the directory where to store the PDB files" )
    parser.add_argument( "-b", "--database", dest = "db",
                         help = "name of the protein structure database from which to download files" )
    parser.add_argument( "-v", "--verbose", dest = "verbose", action = "store_true",
                         help = "print details while downloading" )
    args = parser.parse_args()

    # print fetched arguments
    if args.verbose:
        print( "list: " + args.list )
        print( "dir:  " + args.dir )

    # read the PDB IDs into a python list
    with open( args.list, "rt" ) as f:
        pdb_ids = f.read().splitlines()

    # retrieve PDB files and extract essential headers
    pdbl = PDBList()
    for pdb_id in pdb_ids:
        four_letter = pdb_id[:4].lower()
        mid_letters = pdb_id[1:3].lower()
        
        # create dir if not exists
        if not os.path.exists( args.dir + "/" + mid_letters ):
            os.mkdir( args.dir + "/" + mid_letters )
        
        # retrieve current PDB file and store it under a directory tree
        pdb_file = args.dir + "/" + mid_letters + "/" + four_letter + ".pdb"
        if not os.path.exists( pdb_file ) or os.stat( pdb_file ).st_size == 0:
            print(args.db)
            if args.db == "OPM":
                curl_pdb = "curl http://opm.phar.umich.edu/pdb/" + four_letter + ".pdb -o " + pdb_file 
            else: 
                curl_pdb = "curl https://files.rcsb.org/view/" + four_letter + ".pdb -o " + pdb_file 
            os.system( curl_pdb )
        else:
            print( pdb_file, "already exists and is not empty, skip ..." )

        # if given PDB IDs are actually chains, then get the chain
        if args.chains:
            chain_id = pdb_id[-1]
            chain_file = args.dir + four_letter + chain_id + ".pdb"
            if os.path.exists( chain_file ):
                print( chain_file + " already exists, skip ..." )
            else:
                pdb_parser = PDBParser( PERMISSIVE = 1 )
                structure = pdb_parser.get_structure( four_letter, pdb_file )
                try:
                    chain = structure[0][chain_id]
                except KeyError:
                    print("No chain " + chain_id + " was found in " + pdb_file)
                    sys.exit(1) 
                pdb_writer = PDBIO()
                pdb_writer.set_structure( chain )
                pdb_writer.save( chain_file )


if __name__ == "__main__":
    main()
