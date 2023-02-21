#!/usr/bin/env python3

from argparse import ArgumentParser
from Bio.PDB import PDBParser, PPBuilder
from Bio.Alphabet import ProteinAlphabet
from Bio import Seq, SeqRecord, Align, AlignIO, SeqIO
from Bio import pairwise2


def get_chain_seq(pdb_file, chain_id):
    """

    Parameters
    ----------
    pdb_file : str

    chain_id : str

    Returns
    -------
    Bio.SeqRecord

    """
    pdb_parser = PDBParser(PERMISSIVE=1)
    structure = pdb_parser.get_structure(id='tmp', file=pdb_file)
    try:
        chain = structure[0][chain_id]
    except KeyError:
        print('No chain ' + chain_id + ' was found in ' + pdb_file)
        return None

    # build a polypeptide object from given chain
    ppb = PPBuilder()
    chain_seq = Seq.Seq('', ProteinAlphabet())
    for pp in ppb.build_peptides(chain):
        chain_seq += pp.get_sequence()
    return SeqRecord.SeqRecord(chain_seq, id=chain_id)


def count_mismatches(seq_a, seq_b):
    """
    Hamming distance between seq2 and seq2 on aligned regions.

    Parameters
    ----------
    seq_a
    seq_b

    Returns
    -------

    """
    count = 0
    for a, b in zip(seq_a, seq_b):
        if a != '-' and b != '-' and a != b:
            count += 1
    return count


def pairwise_align(seq_a, seq_b):
    """

    Parameters
    ----------
    seq_a
    seq_b

    Returns
    -------

    """
    # now align the two given sequences
    alignments = pairwise2.align.globalms(
        seq_a, seq_b, 1, -0.5, -10, 0
    )
    multiple_alignment = Align.MultipleSeqAlignment(
        [
            SeqRecord.SeqRecord(seq=Seq.Seq(alignments[0][0]), id='PDB seq'),
            SeqRecord.SeqRecord(seq=Seq.Seq(alignments[0][1]), id='Input seq')
        ]
    )

    mismatches = count_mismatches(alignments[0][0], alignments[0][1])
    if mismatches >= 0.1 * len(alignments[0][1]):
        print('Suspicious alignment between UniProt sequence and PDB sequence, '
              'please check')

    # print the alignment
    print(pairwise2.format_alignment(*alignments[0]))

    # return the alignment
    return multiple_alignment


def create_id_mapping(aligned_seq_a, aligned_seq_b, start_id=1):
    """

    Parameters
    ----------
    aligned_seq_a
    aligned_seq_b
    start_id

    Returns
    -------

    """
    id_mapping = {}
    seq_a_id = start_id - 1
    seq_b_id = 0
    for x, y in zip(aligned_seq_a, aligned_seq_b):
        if x != '-':
            seq_a_id += 1
        if y != '-':
            seq_b_id += 1
        if x != '-':
            if y == '-':
                id_mapping[seq_a_id] = None
            else:
                id_mapping[seq_a_id] = seq_b_id
    return id_mapping


def parse_cmd_args():
    """

    Returns
    -------

    """
    # parse command-line arguments
    parser = ArgumentParser(description='')
    parser.add_argument('-m', '--mapping', dest='mapping', required=False,
                        type=str, help='Residue ID mapping from PDB IDs to '
                                       'target residue IDs.')
    parser.add_argument('-i', '--input', dest='input', required=True, type=str,
                        help='Input PDB file.')
    parser.add_argument('-a', '--alignment', dest='alignment', required=False,
                        type=str, help='Pairwise alignment of the PDB sequence'
                                       'with the target sequence.')
    parser.add_argument('-s', '--sequence', dest='sequence', required=False,
                        type=str, help='Pairwise alignment of the PDB sequence'
                                       'with the target sequence.')
    parser.add_argument('--start-seqres', dest='start_seqres', required=False,
                        type=int, help='Pairwise alignment of the PDB sequence'
                                       'with the target sequence.')
    parser.add_argument('-c', '--chain-id', dest='pdb_chain', required=True,
                        type=str, help='ID of the chain to be renumbered.')
    parser.add_argument('-o', '--output', dest='output', required=True,
                        type=str, help='Output PDB file.')
    return parser.parse_args()


def main():
    """

    Returns
    -------

    """
    # parse command-line arguments
    args = parse_cmd_args()

    # parse the input ID mapping into a dict
    if args.mapping is not None:
        mapping = {}
        with open(args.mapping, 'rt') as ipf:
            for l in ipf:
                old_id, new_id = l.strip().split()
                try:
                    mapping[int(old_id)] = int(new_id)
                except ValueError:
                    print(
                        old_id, 'or', new_id, 'is not a valid '
                        'residue sequence number, skipped.'
                    )
                    continue
    else:  # no ID mapping given, create one by doing pairwise alignment
        if args.alignment is not None:
            print('Using alignment', args.alignment, 'for renumbering.')
            alignment = AlignIO.read(args.alignment, format='fasta')
        else:
            seq_a = get_chain_seq(args.input, args.pdb_chain)
            seq_b = SeqIO.read(args.sequence, format='fasta')

            # now align the two given sequences
            print('Now aligning the two given sequences using global alignment '
                  'and the identity matrix with -10 gap opening penalty...')
            alignment = pairwise_align(seq_a.seq, seq_b.seq)

            # store the alignment to disk file
            AlignIO.write(alignment, './alignment.fasta', 'fasta')

        if args.start_seqres is not None:
            start_seqres = args.start_seqres
        else:
            start_seqres = 1
        mapping = create_id_mapping(
            alignment[0].seq, alignment[1].seq, start_seqres
        )

    # renumber ATOM records of the requested chain
    print('Renumbering records according to alignment:', './alignment.fasta')
    new_records = []
    with open(args.input, 'rt') as ipf:
        for l in ipf:
            if l.startswith('ATOM') and l[21] == args.pdb_chain:
                # renumber the record
                old_id = int(l[22:26])
                new_id = mapping[old_id]
                new_record = l[:22] + '%4s' % new_id + l[26:]
                new_records.append(new_record)
            else:
                new_records.append(l)

    # write renumbered record to output file
    with open(args.output, 'wt') as opf:
        opf.writelines(new_records)

    # print final status
    print('Renumbered PDB records written to', args.output)
    print('Done!')


if __name__ == '__main__':
    main()
