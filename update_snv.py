#!/usr/bin/env python3

from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Alphabet import IUPAC


class SNV:
    """
    A class that abstracts the concept of single nucleotide variants.
    """
    def __init__(self, uniprot=None, position=None, wild_type=None, variant=None):
        """

        Parameters
        ----------
        wild_type : str
            The wild-type residue as a single-letter string.
        variant : str
            The variant residue as a single-letter string.
        position : int
            Position where the mutation occurred.
        """
        self._uniprot = uniprot
        self._position = position
        self._wild_type = wild_type
        self._variant = variant

    @classmethod
    def create_snv(cls, snv_str):
        fields = snv_str.split()
        if len(fields) != 4:
            raise ValueError('Invalid number of fields in ' + fields + '. '
                             'Expected 4 fields.')
        try:
            pos = int(fields[1])
            return cls(fields[0], pos, fields[2], fields[3])
        except ValueError:
            pos = int(fields[2])
            return cls(fields[0], int(fields[2]), fields[1], fields[3])

    @property
    def uniprot(self):
        return self._uniprot

    @uniprot.setter
    def uniprot(self, up):
        if not isinstance(up, str):
            raise TypeError('Invalid argument type! Must be a str.')
        self._uniprot = up

    @property
    def wild_type(self):
        return self._wild_type

    @wild_type.setter
    def wild_type(self, wt):
        if len(wt) != 1:
            raise ValueError('Invalid argument value: ' + wt)
        if wt not in IUPAC.protein.letters:
            raise ValueError('Invalid letter: ' + wt)
        self._wild_type = wt

    @property
    def variant(self):
        return self._variant

    @variant.setter
    def variant(self, vt):
        if not isinstance(vt, str):
            raise TypeError('Argument must be a str! ' + vt + ' is a ' +
                            type(vt))
        if len(vt) != 1:
            raise ValueError('Invalid argument length: ' + len(vt) + '. '
                             'Expected single letter.')
        if vt not in IUPAC.protein.letters:
            raise ValueError('Invalid argument value: ' + vt)
        self._variant = vt

    @property
    def position(self):
        return self._position

    @position.setter
    def position(self, pos):
        if not isinstance(pos, int):
            raise TypeError('Invalid argument type: ' + type(pos) + '. '
                            'Expected an int')
        if pos < 1:
            raise ValueError('Invalid argument value: ' +  pos + '. Expected '
                             'a positive int')
        self._position = pos

    def __str__(self):
        return '%s %d %s %s' % (
            self._uniprot, self._position, self._wild_type, self._variant
        )

    def __repr__(self):
        return 'SNV(%r, %d, %r, %r)' % (
            self._uniprot, self._position, self._wild_type, self._variant
        )


def parse_cmd_arguments():
    # setting up
    parser = ArgumentParser()
    parser.add_argument('-i', '--input', dest='input', type=str,
                        required=True, help='A file that contains the list '
                        'of SNVs for which the numbering needs to be updated.')
    parser.add_argument('-o', '--output', dest='output', type=str,
                        required=True, help='A file in which to write the '
                        'updated SNVs.')
    parser.add_argument('-a', '--alignment', dest='alignment', type=str,
                        required=True, help='Sequence alignment in FASTA '
                        'format. The old sequence must be the first '
                        'sequence record in the FASTA file.')
    # do any checking on command-line arguments here, if necessary
    return parser.parse_args()


def main():
    # parse command-line arguments
    args = parse_cmd_arguments()

    # parse the given FASTA file into Biopython SeqRecords
    records = list(SeqIO.parse(args.alignment, 'fasta', alphabet=IUPAC.protein))

    # there must be two and only two sequence records in the FASTA file
    if len(records) != 2:
        raise SystemExit('There must be two and only two sequence records in '
                         'the FASTA file: ' + args.input)

    # create a list of SNVs from the given file
    with open(args.input, 'rt') as infile:
        snvs = [SNV.create_snv(l.replace(',', ' ').replace(':', ' ')) for l
                in infile if not l.strip().startswith('#')]

    # get all the positions where gaps were inserted
    non_gap_indices = [i for i, a in enumerate(records[0].seq) if a != '-']

    # update SNV numbering
    updated_snvs = list()
    for snv in snvs:
        # add SNVs that do not need numbering updated
        if snv.wild_type == records[1].seq[snv.position - 1]:
            updated_snvs.append(snv)
            continue
        # update numbering
        offset = non_gap_indices[0]
        new_position = snv.position + offset
        # make sure that the amino acid at the new position in the other
        # sequence matches the wild-type amino acid
        if snv.wild_type != records[1].seq[new_position - 1]:
            print('In %s: wild-type residue in the SNV does not match that '
                  'in the new sequence: %s vs %s at %d, skip to the next SNV.' 
                  % (args.input, snv.wild_type, records[1].seq[new_position - 1], new_position))
            continue
        snv.position = new_position
        updated_snvs.append(snv)

    if not updated_snvs:
        raise SystemExit('No valid variant was found in %s.' % args.input)
    else:
        with open(args.output, 'wt') as outfile:
            outfile.write('\n'.join([str(v) for v in updated_snvs]))
            # append a newline character
            outfile.write('\n')


if __name__ == '__main__':
    main()
