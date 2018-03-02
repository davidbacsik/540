'''Script to find the highest weight segment in a sequence.

To see usage:
    python wll.py --help

by David Bacsik'''

import os
import time
import argparse
import glob

def parseArgs():
    '''Takes in command line flags and passes them to variables.'''

    parser = argparse.ArgumentParser(
            description = 'Find the highest weight path in a graph. See '
            'template.csv for graph file format.')
    parser.add_argument('--fasta', help='FASTA file',
            required=True)
    parser.add_argument('--score', help='Score scheme file', required=True)
    parser.add_argument('--out', help='Prefix for output file (optional)',
            required=False, default='')

    return parser

def writeOutput(out_prefix, run_time, path_score, path, first_vertex, last_vertex, seq_content):
    '''
    Writes output text file.

    Args:
        out_prefix (str): Name of output file, from command line
        run_time (float): Total number of seconds the script took to run
        path_score (float): The total weight of the segment
        path (str): Sequence of highest weight segment
        first_vertex (int): Index of the first character in the path
        last_vertex (int): Index of the last character in the path

    Returns:
        out_prefix_output.txt (txt): Text file containing output

    '''

    # Header
    if out_prefix != '':
        out_prefix = out_prefix+'_'
    with open ('{0}output.txt'.format(out_prefix), 'w') as out_file:
        # Header
        out_file.write('Assignment: GS 540 HW2\n')
        out_file.write('Name: David Bacsik\n')
        out_file.write('Email: dbacsik@uw.edu\n')
        out_file.write('Language: Python\n')
        out_file.write('Runtime: '+str(run_time)+' sec\n\n')

        # Seq Summary
        out_file.write('Part X\n')
        out_file.write(seq_content['carrot']+'\n')
        out_file.write('A='+str(seq_content['a'])+'\n')
        out_file.write('C='+str(seq_content['c'])+'\n')
        out_file.write('T='+str(seq_content['t'])+'\n')
        out_file.write('G='+str(seq_content['g'])+'\n')
        out_file.write('N='+str(seq_content['n'])+'\n')
        out_file.write('*='+str(seq_content['length'])+'\n\n')

        # Path Summary
        out_file.write('Score: '+str(path_score)+'\n')
        out_file.write('Begin: '+str(first_vertex)+'\n')
        out_file.write('End: '+str(last_vertex)+'\n')
        out_file.write('Path: '+path+'\n')
        out_file.write('Description: X\n')

def readFasta(fasta_file, score_file):
    '''
    Reads FASTA file and calculates base content.
    Also reads score scheme file and passes it in.

    Args:
        fasta_file (str): Name of FASTA file containing sequence.
        score_file (str): Name of CSV containing scores for each base in FASTA.

    Returns:
        seq (str): Sequence from FASTA file
        seq_content (dict): Dictionary listing name and base content of sequence
        base_scores (dict): Dictionary of edge weights for each base
    '''
# Open FASTA file
    genome = open(fasta_file, 'r')

    seq_list = []

    carrot_line = None

    # Parse out name line
    for line in genome:
        line = line.strip('\n')
        if line[0] == ">":
            carrot_line = line
        else:
            seq_list.append(line)

    genome.close()

    # Concatenate sequence back
    seq = ''.join(map(str, seq_list))
    seq = seq.upper()

    # Calculate genome contents
    # Make content dicitonary
    seq_content = {
    'carrot': carrot_line,
    'a': seq.count('A'),
    'c': seq.count('C'),
    'g': seq.count('G'),
    't': seq.count('T'),
    'n': seq.count('N'),
    'length': len(seq)
    }

    score_scheme = dict()
    # Open score score_scheme
    with open(score_file) as score_file:
        for line in score_file:
            line = line.split(',')
            score_scheme[str(line[0])] = float(line[1])

    return seq, seq_content, score_scheme

def maxSegment(seq, score_scheme):
    '''
    For a sequence, finds and returns maximum scoring segment.

    Args:
        seq (str): Sequence from FASTA file
        score_scheme (dict): Dictionary of edge weights for each base

    Returns:
        score (float): Score of segment
        best_start (int): Starting index of highest scoring segment
        best_end (int): Ending index of highest scoring segment
    '''

    cumul = 0
    max_score = 0
    start = 0
    end = None
    best_start = None
    best_end = None

    # Iterate through each base
    for i in range(len(seq)):
        cumul += score_scheme[seq[i]]
        if cumul <= 0:
            cumul = 0
            start = i + 1
        elif cumul >= max_score:
            max_score = cumul
            best_end = i + 1
            best_start = start

    return max_score, best_start, best_end

def main():
    '''Main body of script.'''
    # Start timer
    start_time = time.time()

    # IO
    # Parse command line arguments
    parser = parseArgs()
    print("\nExecuting {0} at {1}".format(parser.prog, time.asctime()))
    args = parser.parse_args()
    print("\nParsed the following arguments:\n\t{0}".format(
            '\n\t'.join(['{0} = {1}'.format(arg, val) for (arg, val) in
            vars(args).items()])))

    fasta_file = args.fasta
    score_file = args.score

    # Read fasta
    seq, seq_content, score_scheme = readFasta(fasta_file, score_file)

    print(seq_content['carrot'])
    print('A='+str(seq_content['a']))
    print('C='+str(seq_content['c']))
    print('T='+str(seq_content['t']))
    print('G='+str(seq_content['g']))
    print('N='+str(seq_content['n']))
    print('*='+str(seq_content['length']))

    # find highest weight maxSegment
    max_score, best_start, best_end = maxSegment(seq, score_scheme)
    print('\nScore: '+str(max_score))
    print('Begin: '+str(best_start))
    print('End: '+str(best_end))
    print('Path: '+seq[best_start:best_end])

    # End timer
    end_time = time.time()
    runtime = end_time-start_time

    # Write output
    writeOutput(args.out, runtime, max_score, seq[best_start:best_end], best_start, best_end, seq_content)

if __name__ == '__main__':
    main()
