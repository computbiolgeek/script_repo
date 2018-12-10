#!/usr/bin/env python3

import sys
from argparse import ArgumentParser
from Bio.PDB import PDBParser, is_aa
import numpy as np
import pandas as pd
import matplotlib as plt
import pylab
from itertools import combinations_with_replacement


def compute_sc_centroid(residue):
    """
        C = (v1 + ... + vn) / n
    """
    # return CA coordinates if residue is GLY
    if residue.get_resname() == "GLY":
        return residue['CA'].get_coord()
    # for other residues, return sidechain centroid
    backbone_atoms = {'N', 'CA', 'C', 'O'}
    centroid = np.array( [0.0, 0.0, 0.0] )
    number_sidechain_atoms = 0
    for atom in residue.get_list():
        if atom.get_name() not in backbone_atoms:
            centroid += atom.get_coord()
            number_sidechain_atoms += 1
    return centroid / number_sidechain_atoms


def compute_distance(res_a, res_b):
    '''
    The Euclidean distance between the centroid of residue A and the centroid of residue B.
    '''
    centroid_a = compute_sc_centroid(res_a)
    centroid_b = compute_sc_centroid(res_b)
    return np.linalg.norm(centroid_a - centroid_b)


def compute_shortest_distance(res_a, res_b):
    '''
    Returns the smallest of the distances between the heavy atoms of residue a and those of residue b.
    '''
    distances = [np.linalg.norm(atom_a - atom_b) for atom_a in res_a.get_list() for atom_b in res_b.get_list()]
    return min(distances)

def px2pt(p):
    '''
    '''
    return p * 72. / 96


def init_spines(hidden=[]):
    ax = pylab.gca()

    all_spines = ['bottom', 'left', 'top', 'right', 'polar']

    for spine in all_spines:
        if spine in hidden:
            ax.spines[spine].set_visible(False)
        else:
            try:
                ax.spines[spine].set_visible(True)
                ax.spines[spine].set_linewidth(px2pt(0.75))
            except KeyError:
                pass


def init_pylab(font_kwargs={}):
    plt.rc('lines', linewidth=px2pt(2))
    plt.rc('xtick', **{'direction': 'out'})
    plt.rc('ytick', **{'direction': 'out'})
    plt.rc('legend', frameon=False, fontsize=font_kwargs['size'], numpoints=1)
    plt.rc('font', **font_kwargs)

    pylab.tick_params(axis='x', which='both', top='off')
    pylab.tick_params(axis='y', which='both', right='off')

    init_spines()


def compute_dist_matrix(residues, symmetric=True):
    '''
    '''
    dist_mat = np.zeros((len(residues), len(residues)), dtype='float64')

    #
    dist_mat[:] = np.NaN

    # residue pair indices
    pair_indeces = combinations_with_replacement(range(len(residues)), 2)

    for i, j in pair_indeces:
        res_i = residues[i]
        res_j = residues[j]
        dist_ij = compute_shortest_distance(res_i, res_j)
        dist_mat[i, j] = dist_ij 
        
        if symmetric:
            dist_mat[j, i] = dist_ij
    return dist_mat


def main():
    # command line argument parser
    parser = ArgumentParser(description='Compute a distance matrix given two list of residues.')
    parser.add_argument('-p', '--pdb', dest='pdb', help='input file in PDB format')
    parser.add_argument('-o', '--output', dest='output', help='name to the file to which to write the distance matrix')

    # parse command line arguments
    args = parser.parse_args()

    # parse the given PDB file
    pdb_parser = PDBParser(PERMISSIVE=1)
    structure = pdb_parser.get_structure(id='pdb', file=args.pdb)
    model = structure[0]

    # get all the residues
    residues = [r for r in model.get_residues() if is_aa(r)]

    # compute pair distances between residues in the first chain and residues in the second chain
    dist_mat = compute_dist_matrix(residues)

    # make a contact map
    init_pylab({'family': 'sans-serif', 'size': 24})
    init_spines(hidden=['bottom', 'left', 'top', 'right'])
    pylab.gcf().set_figwidth(10)
    pylab.gcf().set_figheight(10)
    pylab.xlim(1, len(residues) + 1)
    pylab.ylim(1, len(residues) + 1)
    pylab.xlabel('Residue index')
    pylab.ylabel('Residue index')
    ax, fig = pylab.gca(), pylab.gcf()
    cmap = plt.cm.jet
    # map_object = pylab.pcolormesh(dist_mat, shading='flat', edgecolors=None, cmap=cmap)
    contour_map = plt.pyplot.contour(dist_mat, levels=list(range(0, 55, 6)), cmap=cmap)
    map_object = ax.imshow(dist_mat, interpolation='spline16', cmap=cmap)
    box = ax.get_position()
    pad, width = 0.02, 0.02
    cax = fig.add_axes([box.xmax + pad, box.ymin, width, box.height])
    color_bar = pylab.colorbar(map_object, drawedges=False, cax=cax)
    color_bar.set_ticks(ticks = list(range(0, 55, 6)), update_ticks=True)
    color_bar.outline.set_visible(False)
    pylab.ylabel('Distance (' + r'$\AA$' + ')')
    # ax.set_title(args.pdb.split('.')[0] + 'contact map', fontweight='bold')

    # save the map
    pylab.savefig(args.output, bbox_inches='tight', dpi=300)



if __name__ == '__main__':
    main()

