#!/usr/bin/env python3

"""
This script reads in a fasta file containing two records and does pairwise
alignment on them using the identidy matrix. -10 points are deducted from
opening a gap.
"""

# imports from the standard library
from argparse import ArgumentParser
# imports from external libraries
from Bio import pairwise2, Seq, SeqIO, SeqRecord, Align, AlignIO
from Bio.Alphabet import IUPAC


def count_mismatches(seq_a, seq_b):
    """
    Hamming distance between seq2 and seq2 on aligned regions.

    Parameters
    ----------
    seq_a : str
        The first sequence as a str.
    seq_b : str
        The second sequence as a str.

    Returns
    -------
    int
        The Hamming distance between the given two sequences.
    """
    count = 0
    for a, b in zip(seq_a, seq_b):
        if a != '-' and b != '-' and a != b:
            count += 1
    return count


def parse_cmd_arguments():
    # setting up
    parser = ArgumentParser()
    parser.add_argument('-a', '--sequence-a', dest='seq_a',
                        help='a fasta file containing the first sequence')
    parser.add_argument('-b', '--sequence-b', dest='seq_b',
                        help='a fasta file containing the second sequence')
    parser.add_argument('-o', '--output', dest='output',
                        help='a name for the output fasta file')
    # do any checking on command-line arguments here, if necessary
    return parser.parse_args()


def main():
    # parse command-line arguments
    args = parse_cmd_arguments()

    # parse amino acid sequences
    record_a = next(SeqIO.parse(args.seq_a, 'fasta', alphabet=IUPAC.protein))
    record_b = next(SeqIO.parse(args.seq_b, 'fasta', alphabet=IUPAC.protein))

    # print the given two sequences
    print(record_a, '\n')
    print(record_b, '\n')

    # now align the two given sequences
    print('Now aligning the two given sequences using global alignment and '
          'the identity matrix with -10 gap openning penalty ...')
    alignments = pairwise2.align.globalms(record_a.seq, record_b.seq, 1, -0.5,
                                          -10, 0)
    multiple_alignment = Align.MultipleSeqAlignment([
        SeqRecord.SeqRecord(seq=Seq.Seq(alignments[0][0]), id=record_a.id,
                            description=record_a.description),
        SeqRecord.SeqRecord(seq=Seq.Seq(alignments[0][1]), id=record_b.id,
                            description=record_b.description)
    ])
    AlignIO.write(multiple_alignment, args.output, 'fasta')
    print('The following alignment has been written to ' + args.output, '\n')

    # flag suspicious alignment
    mismatches = count_mismatches(alignments[0][0], alignments[0][1])
    if mismatches >= 0.1 * len(alignments[0][1]):
        print('Suspicious alignment between UniProt sequence and PDB sequence,'
              ' please check', args.output)

    # print the alignment
    print(pairwise2.format_alignment(*alignments[0]))


if __name__ == '__main__':
    main()
