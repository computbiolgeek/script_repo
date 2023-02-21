#!/usr/bin/env python3

from argparse import ArgumentParser


def main():
    """
    Convert a multiple sequence alignment in Clustal format into Grishin format.
    """
    # parse command-line arguments
    parser = ArgumentParser()
    parser.add_argument('-i', '--input', dest='input', type=str, required=True,
            help='Input sequence alignment in Clustal format.')
    parser.add_argument('-n', '--n-seqs', dest='nseqs', type=int, required=True,
            help='Number of sequences aligned.')
    args = parser.parse_args()

    # parse Clustal Omega alignment
    seq_ids = []
    nseqs = args.nseqs
    seqs = [''] * args.nseqs
    with open(args.input, 'rt') as ipf:
        for i, l in enumerate(ipf):
            # each chunk of alignment has two extra lines
            n = i % (nseqs + 2)
            if n < nseqs:
                try:
                    seq_id, seq, _ = l.split()
                except ValueError:
                    seq_id, seq = l.split()
                seqs[n] += seq
                if len(seq_ids) < nseqs:
                    seq_ids.append(seq_id)

    # write pairwise alignment in Grishin format
    for i in range(1, nseqs):
        out_file = seq_ids[0] + '_' + seq_ids[i] + '.grishin'
        with open(out_file, 'wt') as opf:
            opf.write('## ' + seq_ids[0] + ' ' + seq_ids[i] + '.pdb\n')
            opf.write('#\n')
            opf.write('score from program: 0\n')
            opf.write('0 ' + seqs[0] + '\n')
            opf.write('0 ' + seqs[i])
        print('Wrote a pairwise alignment in Grishin format written to', out_file)


if __name__ == '__main__':
    main()
