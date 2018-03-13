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
import operator
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
# Reads
reads = list()

# Graph stuff
all_verts = None
all_edges = None
edge_histo = dict()

# Blosum stuff
gap_penalty = int()
blosum_lookup = dict()
blosom_array = None

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
    read.edges = np.empty(read.length+1, dtype='object')

    # Iterate and fill verts and edges.
    for i in range(read.length):
        read.verts[i] = i

    for i in range(read.length):
        read.edges[i+1] = read.seq[i]

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

class vertex(object):
    '''Vertex in 3D graph

    Args:
        score (int) = Best score at vertex
        parent (tuple) = Coordinates of parent in max path
    '''

    def __init__(self):
        '''
        Initilize relevant variables.
        '''

        self.score = 0
        self.parent = None
        self.alignment = None

class edge(object):
    ''' Edge in 3D graph

    Args:
        start (tuple) = Coordinates of starting vertex
        end (tuple) = Coordinates of ending vertex
        parent_score (int) = Score of vertex the edge originates from
        alignment (str) = Alignment along this edge
        score (int) = total score for Vertex
        xy_score (int) = score of xy substitution
        xz_score (int) = score of xz substittuion
        yz_score (int) = score of yz substitution
        edge_score (int) = score the subsitution adds to the path
    '''

    def __init__(self):
        '''
        Initilize relevant variables.
        '''

        self.start = None
        self.parent_score = None
        self.end = None
        self.alignment = None
        self.score = None
        self.xy_score = None
        self.xz_score = None
        self.yz_score = None
        self.edge_score = None

def makeGraph():
    '''
    Makes 3D graph from 3 input sequences.

    Args:
        reads (list): List of read info, including graphs of vertices and edges.

    Returns:
        all_verts (np.array): Array of all vertices, with dimensions of sequence length+1 for each read.
        all_edges (dict): List of all vertices, including vertices that. Total length = all seq length
        edge_histo (dict): Count of how many times each edge comes up.
    '''

    # Global vars
    global reads
    global all_verts
    global all_edges
    global edge_histo
    global gap_penalty

    # Read in individual graphs
    x = reads[0]
    y = reads[1]
    z = reads[2]

    # Make Vertex Graph:
    all_verts = np.empty((x.length+1,y.length+1,z.length+1), dtype=object)

    # Generate edges
    all_edges = np.empty(x.length*y.length*z.length*7, dtype=object)
    count = 0

    for i in range(x.length+1):
        for j in range(y.length+1):
            for k in range(z.length+1):
                current_vertex = vertex()
                all_verts[i,j,k] = current_vertex
                if i != 0 and j != 0 and k != 0:
                    current_edges = dict()
                    ce_scores = dict()
                    for i2 in [i,i-1]:
                        for j2 in [j,j-1]:
                            for k2 in [k,k-1]:
                                if (i2 != i or j2 != j or k2 != k):
                                    if i == i2:
                                        x_char = '-'
                                    else:
                                        x_char = x.seq[i-1]
                                    if j == j2:
                                        y_char = '-'
                                    else:
                                        y_char = y.seq[j-1]
                                    if k == k2:
                                        z_char = '-'
                                    else:
                                        z_char = z.seq[k-1]

                                    # Edge
                                    current_edge = edge()
                                    current_edge.start = (i2,j2,k2)
                                    current_edge.end = (i,j,k)

                                    if x_char != '-' and y_char != '-':
                                        current_edge.xy_score = int(calculateScore(x_char,y_char))
                                    else:
                                        current_edge.xy_score = int(gap_penalty)
                                    if x_char != '-' and z_char != '-':
                                        current_edge.xz_score = int(calculateScore(x_char,z_char))
                                    else:
                                        current_edge.xz_score = int(gap_penalty)
                                    if y_char != '-' and z_char != '-':
                                        current_edge.yz_score = int(calculateScore(y_char,z_char))
                                    else:
                                        current_edge.yz_score = int(gap_penalty)

                                    current_edge.edge_score = current_edge.xy_score + current_edge.xz_score + current_edge.yz_score
                                    current_edge.parent_score = all_verts[i2,j2,k2].score
                                    current_edge.score = current_edge.parent_score + current_edge.edge_score

                                    current_edge.alignment = '{0}{1}{2}'.format(x_char,y_char,z_char)

                                    # Add to All Edges
                                    all_edges[count] = current_edge
                                    count += 1

                                    # Add to current_edges
                                    current_edges[current_edge.alignment] = current_edge
                                    ce_scores[current_edge.alignment] = current_edge.score

                                    # Add to Histo
                                    if current_edge.alignment in edge_histo.iterkeys():
                                        edge_histo[current_edge.alignment] += 1
                                    else:
                                        edge_histo[current_edge.alignment] = 1

                    # Pick Parent
                    best_edge = max(ce_scores.iterkeys(), key=(lambda key: ce_scores[key]))
                    current_vertex.parent = current_edges[best_edge].start
                    current_vertex.score = current_edges[best_edge].score
                    current_vertex.alignment = current_edges[best_edge].alignment
                    all_verts[i,j,k] = current_vertex

def writeOutput(runtime):
    '''
    Writes graph to text file.

    Args:
        all_verts (np.array): Empty array of all vertices. Dimensions of 3 sequence lengths.
        all_edges (np.array): Sparse matrix of all edges. Edges are stored as single character.
        runtime (float): Runtime in seconds
    '''

    with open('results.txt', 'w') as out_file:
        out_file.write('Run time: {0}\n\n'.format(runtime))
        for i in range(reads[0].length+1):
            for j in range(reads[1].length+1):
                for k in range(reads[2].length+1):
                    out_file.write('Vertex: {0}, Parent: {1}, Score: {2}\n'.format((i,j,k), all_verts[i,j,k].parent, all_verts[i,j,k].score))

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

    # Open blosum file and load scores
    parseBlosum(args.b, args.g)

    # Open FASTA and collect reads.
    parseFASTA(args.f)

    # Convert reads to graph format
    for x in reads:
        seqToGraph(x)

    # Make graph from all 3 sequences.
    makeGraph()

    # End timer
    end_time = time.time()
    runtime = end_time-start_time

    # Write output
    writeOutput(runtime)

    print('Done.')

if __name__ == '__main__':
    main()
