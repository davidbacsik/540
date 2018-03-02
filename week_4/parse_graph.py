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
    parser.add_argument('--b', help='2-D Matrix of amino acid substitution scores.',
            required=False, default='blosum62.txt')
    parser.add_argument('--out', help='Prefix for output file (optional)',
            required=False, default='')

    return parser

# Global args
reads = list()
subs = list()
gap_penalty = int()

class sub(object):
    '''
    An individual amino acid substitution and its score.

    Args:
        initial (str): Character representing original amino acid
        final (str): Character representing derived amino acid
        score (int): Value for this substitution score
    '''

    def __init__(self):
        '''
        Initilize relevant variables.
        '''

        self.initial = str()
        self.final = str()
        self.score = int()

def parseBlosum(blosum_file):
    '''
    Reads in blosum substitution file and organizes scores for calcuations.

    Args:
        blosum_file (str): name of blosum substitution file

    Returns:
        sub.initial (str): Character representing original amino acid
        sub.final (str): Character representing derived amino acid
        sub.score (int): Value for this substitution score
    '''

    # Global vars
    global subs
    global gap_penalty

    with open(blosum_file, 'r') as f:
        for line in f:
            if 'Gap' in line:
                line = line.split()
                for item in line:
                    if type(item) == int:
                        print(item)


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
    read.verts = np.empty(shape=(1,read.length), dtype=int)

    # Make edges
    read.edges = np.empty(shape=(1,read.length), dtype='object')

    # Iterate and fill verts and edges.
    for i in range(read.length):
        read.verts[0,i] = i
        read.edges[0,i] = read.seq[i]

# MAIN
def main():
    '''Main body of script.'''

    # Global vars
    global reads
    global subs

    # Start timer
    start_time = time.time()

    # Parse command line arguments
    parser = parseArgs()
    print("\nExecuting {0} at {1}".format(parser.prog, time.asctime()))
    args = parser.parse_args()
    print("\nParsed the following arguments:\n\t{0}".format(
            '\n\t'.join(['{0} = {1}'.format(arg, val) for (arg, val) in
            vars(args).items()])))

    # Open Blosum Matrix and parse
    parseBlosum(args.b)

    # Open graph and iterate through reads
    parseFASTA(args.f)

    # Convert reads to graph format
    for x in reads:
        seqToGraph(x)


    # End timer
    end_time = time.time()
    runtime = end_time-start_time

    # Write output
    #writeOutput(args.out, runtime, path_score, first_vertex, last_vertex, edge_path)

if __name__ == '__main__':
    main()
