#!/usr/bin/env python3

import os
from argparse import ArgumentParser


def main():
    """
    Create a membrane span file from PPM predictions for Rosetta.
    """
    parser = ArgumentParser()
    parser.add_argument('-i', '--input', dest='input', type=str, required=True,
            help='Transmembrane segnments predicted by the PPM Server.')
    parser.add_argument('-n', '--num-res', dest='len', type=str, required=True,
            help='Number of residues of the protein sequence.')
    parser.add_argument('-o', '--output', dest='output', type=str, required=True,
            help='Output file to which to write the span file for Rosetta.')
    args = parser.parse_args()

    # parse PPM predictions into a list
    segments = []
    with open(args.input, 'rt') as ipf:
        for l in ipf:
            segments += [s.replace('(', ':').replace(')', '').split(':')[1] \
                    for s in l.strip().split(',')]
    
    print('PPM Server detected ', len(segments), ' transmembrane segments.')
    print('Now writing span file to ' +  args.output)

    # write segments into a Rosetta span file
    basename = os.path.basename(args.output).split('.')[0]
    with open(args.output, 'wt') as opf:
        opf.write('TM region definitions for ' + basename + ' using PPM\n') 
        opf.write(str(len(segments)) + ' ' + args.len + '\n')
        opf.write('antiparallel\n')
        opf.write('n2c\n')
        formatted_segments = []
        for s in segments:
            start, end = [x.strip() for x in s.strip().split('-')]
            formatted_segments.append(start + '\t' + end + '\t' + start + '\t' + end)
        opf.write('\n'.join(formatted_segments))

    print('span file written successfully to ' + args.output)


if __name__ == '__main__':
    main()
