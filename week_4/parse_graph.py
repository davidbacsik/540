'''
Program to convert FASTA sequences to edit graphs. Constrained to only 3 input sequences.

To see usage:
    python process_graph.py --help

By David Bacsik
'''

import os
import time
import argparse
import glob
import numpy as np

def parseArgs():
    '''Takes in command line flags and passes them to variables.'''

    parser = argparse.ArgumentParser(
            description = 'Find the highest weight path in a graph. See '
            'template.csv for graph file format.')
    parser.add_argument('--f', help='FASTA file containing exactly 3 sequences',
            required=True)
    parser.add_argument('--out', help='Prefix for output file (optional)',
            required=False, default='')

    return parser

# Global args
reads = list()
all_verts = None
all_edges = None

class read(object):
    '''
    An individual read from a FASTA file

    Args:
        name (str): name
        seq (str): sequence of string
        length (int): length
    '''

    def __init__(self):
        '''
        Initilize relevant variables.
        '''

        # Basic from FASTA
        self.name = str()
        self.seq = str()
        self.length = int()

def parseFASTA(fasta_file):
    '''
    Takes in FASTA file and returns a parsed sequence, name, and length.

    Args:
        fasta_file (str): name of FASTA file containing exactly 3 sequences.

    Returns:
        read.name (str): short name of read
        read.seq (str): sequence of read
        read.length (int): length of sequence
    '''

    # Global vars
    global reads

    with open(fasta_file, 'r') as f:
        current_read = read()
        for line in f:
            line = line.strip('\n')
            if line != '':
                if line[0] == '>':

                    # Parse Short Name
                    line = line.split('|')
                    name = line[2]
                    name = name.split()

                    # Set read name
                    current_read.name = name[0]
                else:
                    current_read.seq = line
                    current_read.length = len(line)
                    reads.append(current_read)
                    current_read = read()

def seqToGraph(read):
    '''
    Converts a FASTA sequence to an edit graph. Edges are residues.

    Args:
        read.name (str): name of sequence
        read.seq (str): sequence

    Returns:
        read.verts (np.array): array of vertices
        read.edges (np.array): array of edges
    '''

    # Make vertices
    read.verts = np.empty(read.length, dtype=int)

    # Make edges
    read.edges = np.empty(read.length, dtype='object')

    # Iterate and fill verts and edges.
    for i in range(read.length):
        read.verts[i] = i

    for i in range(read.length):
        read.edges[i] = read.seq[i]

def mergeArrays():
    '''
    Merges each read's graph into a 3-D graph.

    Args:
        reads (list): List of read info, including graphs of vertices and edges.

    Returns:
        all_verts (np.array): Array of all vertices, with dimensions of sequence length for each read.
        all_edges (dict): List of all vertices, including vertices that
    '''

    # Global vars
    global reads
    global all_verts
    global all_edges

    # Read in individual graphs
    x = reads[0]
    y = reads[1]
    z = reads[2]

    # Merge vertices
    all_verts = np.empty(shape=(x.length+1, y.length+1, z.length+1), dtype='object')

    # Merge edges
    all_edges = dict()

    for i in range(0, x.length):
        all_edges['edge_x_{0}'.format(i)] = [(i,0,0), (i+1,0,0), x.seq[i]]
    for i in range(1, y.length):
        all_edges['edge_y_{0}'.format(i)] = [(0,i,0), (0,i+1,0), y.seq[i]]
    for i in range(1, z.length):
        all_edges['edge_z_{0}'.format(i)] = [(0,0,i), (0,0,i+1), z.seq[i]]

    # Add all pairwise combos:
    # XY
    for i in range(0, max(x.length, y.length)):
        if i >= y.length:
            all_edges['edge_xy_{0}'.format(i)] = [(i,y.length,0), (i+1,y.length,0), '{0}--'.format(x.seq[i])]
        elif i >= x.length:
            all_edges['edge_xy_{0}'.format(i)] = [(x.length,i,0), (x.length,i+1,0), '-{0}-'.format(y.seq[i])]
        else:
            all_edges['edge_xy_{0}'.format(i)] = [(i,i,0), (i+1,i+1,0), '{0}{1}-'.format(x.seq[i],y.seq[i])]

    # XZ
    for i in range(0, max(x.length, z.length)):
        if i >= z.length:
            all_edges['edge_xz_{0}'.format(i)] = [(i,0,z.length), (i+1,0,z.length), '{0}--'.format(x.seq[i])]
        elif i >= x.length:
            all_edges['edge_xz_{0}'.format(i)] = [(x.length,0,i), (x.length,0,i+1), '--{0}'.format(z.seq[i])]
        else:
            all_edges['edge_xz_{0}'.format(i)] = [(i,0,i), (i+1,0,i+1), '{0}-{1}'.format(x.seq[i],z.seq[i])]

    # YZ
    for i in range(0, max(y.length, z.length)):
        if i >= z.length:
            all_edges['edge_yz_{0}'.format(i)] = [(0,i,z.length), (0,i+1,z.length), '-{0}-'.format(y.seq[i])]
        elif i >= y.length:
            all_edges['edge_yz_{0}'.format(i)] = [(0,y.length,i), (0,y.length,i+1), '--{0}'.format(z.seq[i])]
        else:
            all_edges['edge_yz_{0}'.format(i)] = [(0,i,i), (0,i+1,i+1), '-{0}{1}'.format(y.seq[i],z.seq[i])]

    # Add triplet combos:
    # XYZ
    for i in range(0, max(x.length, y.length, z.length)):
        x_end = True
        y_end = True
        z_end = True
        x_char = '-'
        y_char = '-'
        z_char = '-'
        if i < x.length:
            x_end = False
            x_char = x.seq[i]
        if i < y.length:
            y_end = False
            y_char = y.seq[i]
        if i < z.length:
            z_end = False
            z_char = z.seq[i]

        all_edges['edge_xyz_{0}'.format(i)] = [([i,x.length][x_end],[i,y.length][y_end],[i,z.length][z_end]),
                ([i+1,x.length][x_end],[i+1,y.length][y_end],[i+1,z.length][z_end]),
                '{0}{1}{2}'.format(x_char,y_char,z_char)]

def writeOutput():
    '''
    Writes graph to text file.

    Args:
        all_verts (np.array): Empty array of all vertices. Dimensions of 3 sequence lengths.
        all_edges (np.array): Sparse matrix of all edges. Edges are stored as single character.
    '''

# MAIN
def main():
    '''Main body of script.'''

    # Global vars
    global reads

    # Start timer
    start_time = time.time()

    # Parse command line arguments
    parser = parseArgs()
    print("\nExecuting {0} at {1}".format(parser.prog, time.asctime()))
    args = parser.parse_args()
    print("\nParsed the following arguments:\n\t{0}".format(
            '\n\t'.join(['{0} = {1}'.format(arg, val) for (arg, val) in
            vars(args).items()])))

    # Open FASTA and collect reads.
    parseFASTA(args.f)

    # Convert reads to graph format
    for x in reads:
        seqToGraph(x)

    # Merge arrays
    mergeArrays()

    # End timer
    end_time = time.time()
    runtime = end_time-start_time

    # Write output
    #writeOutput(args.out, runtime, path_score, first_vertex, last_vertex, edge_path)

if __name__ == '__main__':
    main()
