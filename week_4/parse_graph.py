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
    parser.add_argument('--g', help='Gap Penalty. Should be negative number',
            required=False, default='-6')
    parser.add_argument('--out', help='Prefix for output file (optional)',
            required=False, default='')

    return parser

# Global args
reads = list()
subs = list()
gap_penalty = int()
blosum_lookup = dict()
blosom_array = None

def parseBlosum(blosum_file, gap_value):
    '''
    Reads in blosum substitution file and organizes scores for calcuations.

    Args:
        blosum_file (str): name of blosum substitution file
        gap_value (int): Gap penalty. Should be negative number

    Returns:
        gap_penalty (int): Gap penalty. Should be negative number
        blosum_array (np.array): Substitution Matrix
        blosum_lookup (dict): Lookup table converting amino acid symbol to index value for blosum array

    '''

    # Global vars
    global gap_penalty
    global blosum_lookup
    global blosum_array

    # Parse Gap penalty
    gap_penalty = gap_value

    # Parse blosum matrix
    with open(blosum_file, 'r') as f:
        count = 0
        for line in f:
            blosum_lookup[line[0]] = count - 1
            count += 1

    with open(blosum_file, 'r') as f:
        blosum_array = np.loadtxt(f, dtype=object, skiprows=1, usecols=range(1,len(blosum_lookup)+1))


def calculateScore(aa_init, aa_final):
    '''
    Takes initial amino acid, subsistuted amino acid, and looks up substitution score.

    Args:
        aa_init (str): Initial amino acid
        aa_final (str): Substituted amino acid

    Returns:
        sub_score (int): Score looked up from blosum_array
    '''

    sub_score = int()

    aa_init_index = blosum_lookup[aa_init]
    aa_final_index = blosum_lookup[aa_final]

    sub_score = blosum_array[aa_init_index,aa_final_index]

    return sub_score

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
    parseBlosum(args.b, args.g)

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
