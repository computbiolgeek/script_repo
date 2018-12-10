#!/usr/bin/env python3

"""
This script takes two lists of records and removes records from the first
list that are also in the second list.
"""

import sys


def main():
    # get command-line arguments
    args = sys.argv

    # if -h is given on the command line, print help information
    if '-h' in args:
        print('usage: rm_comm.py <first list> <second list>')
        sys.exit(0)

    # read in the first list
    with open(args[1], 'rt') as f:
        first = [r.strip() for r in f.readlines()]

    # read in the second list
    with open(args[2], 'rt') as f:
        second = [r.strip() for r in f.readlines()]

    # remove elements in first_set that are also in the second_set from
    # first_set
    first_new = [e for e in first if e not in second]

    # print remaining elements
    for e in first_new:
        print(e)


if __name__ == "__main__":
    main()
