'''Script to finds D-segments in a sequence which pass a certain threshold score.

To see usage:
    python wll.py --help

by David Bacsik'''

import os
import time
import argparse
import glob
import operator
import gzip
import numpy as np

def parseArgs():
    '''Takes in command line flags and passes them to variables.'''

    parser = argparse.ArgumentParser(
            description = 'Find the highest weight path in a graph. See '
            'template.TSV for graph file format.')
    parser.add_argument('--counts', help='Formatted file of read count starts',
            required=True)
    parser.add_argument('--score', help='Score scheme file', required=True)
    parser.add_argument('--out', help='Prefix for output file (optional)',
            required=False, default='')

    return parser

def writeOutput(out_prefix, run_time, elevated_segs, elevated_stats, normal_stats):
    '''
    Writes output text file.

    Args:
        out_prefix (str): Name of output file, from command line
        run_time (float): Total number of seconds the script took to run
        elevated_segs (list): List of dictionaries containing all the elevated segments.
                                Format is: start, end, score
        elevated_stats (dict): Dictionary containing summary of stats for elevated segments:
                                total = number of segments, 0 = number of positions with no read starts,
                                1 = number of positions with 1 read start, 2 = number of positions with 2 read starts,
                                3 = number of positions with 3 or more read starts
        normal_stats (dict): Dictionary containing summary of stats for normal segments:
                                total = number of segments, 0 = number of positions with no read starts,
                                1 = number of positions with 1 read start, 2 = number of positions with 2 read starts,
                                3 = number of positions with 3 or more read starts

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

        # Segment Histogram
        out_file.write('Segment Histogram:\n')
        out_file.write('Non-Elevated CN Segments={0}\n'.format(normal_stats['total']))
        out_file.write('Elevated CN Segments={0}\n\n'.format(elevated_stats['total']))

        # Elevated CN segment List:
        out_file.write('Elevated CN Segment List:\n')
        for segment in elevated_segs:
            out_file.write('{0} {1} {2}\n'.format(segment['start'],segment['end'],segment['score']))
        out_file.write('\n')

        # Annotations
        out_file.write('Annotations:\n\n')

        for i in range(0,3):
            out_file.write('Start: {0}\n'.format(elevated_segs[i]['start']))
            out_file.write('End: {0}\n'.format(elevated_segs[i]['end']))
            out_file.write('Description: X\n\n')

        # Read Start Histograms
        out_file.write('Read Start Histogram for Non-Elevated CN Segments:\n')
        out_file.write('0={0}\n'.format(normal_stats[0]))
        out_file.write('1={0}\n'.format(normal_stats[1]))
        out_file.write('2={0}\n'.format(normal_stats[2]))
        out_file.write('>=3={0}\n\n'.format(normal_stats[3]))

        out_file.write('Read Start Histogram for Elevated CN Segments:\n')
        out_file.write('0={0}\n'.format(elevated_stats[0]))
        out_file.write('1={0}\n'.format(elevated_stats[1]))
        out_file.write('2={0}\n'.format(elevated_stats[2]))
        out_file.write('>=3={0}\n'.format(elevated_stats[3]))

class highSegments(object):
    '''
    All high segments in a given set of Copy Number Variants.
    '''

    def __init__(self):
        '''
        Initilize relevant variables.
        '''

        self.counts_file = None
        self.score_file = None
        self.counts = list()
        self.score_dict = dict()
        self.D = None
        self.S = None

        self.elevated_segs = list()
        self.normal_segs = list()
        self.elevated_stats = {'total': 0,
                                0: 0,
                                1: 0,
                                2: 0,
                                3: 0}
        self.normal_stats = {'total': 0,
                                0: 0,
                                1: 0,
                                2: 0,
                                3: 0}

    def setInFiles(self, counts_in_file, score_in_file):
        '''
        Sets highSegments file variables properly.

        Args:
            counts_in_file (str): Name of file containing chromosome position and counts.
            score_in_file (str): Name of TSV containing scores for each count scheme.

        Returns:
            counts_file (str): Name of file containing chromosome position and counts.
            score_file (str): Name of TSV containing scores for each count scheme.
        '''

        assert(counts_in_file != None),'No counts file.'
        print('Reading in counts file:')
        self.counts_file = counts_in_file
        print(self.counts_file)
        print('\n')

        assert(score_in_file != None),'No score file.'
        print('Reading in score file:')
        self.score_file = score_in_file
        print(self.score_file)
        print('\n')


    def readCounts(self):
        '''
        Reads formatted TSV file with counts of how many reads start at that
            chromosomal position.
        Also reads score scheme file and passes it in.

        Args:
            counts_file (str): Name of file containing chromosome position and counts.

        Returns:
            counts (list): List of dicts; chromosome, position, and number of counts starting at that position.
        '''

        if '.gz' in self.counts_file:
            with gzip.open(self.counts_file) as cf:
                for line in cf:
                    line = line.strip('\n')
                    line = line.split('\t')
                    line_dict = {'chr':line[0],
                                'pos':line[1],
                                'count':line[2]}
                    self.counts.append(line_dict)
        else:
            with open(self.counts_file) as cf:
                for line in cf:
                    line = line.strip('\n')
                    line = line.split('\t')
                    line_dict = {'chr':line[0],
                                'pos':line[1],
                                'count':line[2]}
                    self.counts.append(line_dict)


    def readScores(self):
        '''
        Reads formatted score file of how to analayze counts file.

        Args:
            score_file (str): Name of file containing chromosome position and counts.

        Returns:
            score_dict (dict): Dictionary of reference scores. Describes how to handle each number.
                        Format is:
                        {(Count: (score, str('=', '>', '<', '>=', '<='))}
        '''

        with open(self.score_file) as sf:
            for line in sf:
                line = line.strip('\n')
                line = line.split(',')
                if line[0] == 'score':
                    self.score_dict[int(line[1])] = (str(line[3]), np.float64(line[2]))
                elif line[0] == 'D':
                    self.D = np.float64(line[1])
                elif line[0] == 'S':
                    self.S = np.float64(line[1])
                else:
                    pass

        print('D: ',self.D)
        print('S: ',self.S)
        print('Score Scheme: ')
        print(self.score_dict)
        print('\n')


    def scorePos(self, position):
        '''
        Takes a position's copy number and returns its score.

        Args:
            counts (list): List of dicts; chromosome, position, and number of counts starting at that position.
            score_dict (dict): Dictionary of reference scores. Describes how to handle each number.
                        Format is:
                        {(Count: (score, str('=', '>', '<', '>=', '<='))}
            position (int): Specific position to look at in counts

        Returns:
            pos_score (float): Score for that position.
        '''

        ops = {'=': operator.eq,
                '<': operator.lt,
                '>': operator.gt,
                '<=': operator.le,
                '>=': operator.ge}

        pos_score = np.float64()

        read_starts = int(self.counts[position]['count'])

        for key, value in self.score_dict.iteritems():
            if ops[value[0]](read_starts, key):
                pos_score = np.float64(value[1])
                break

        return pos_score, read_starts

    def makeNegatives(self):
        '''
        Finds negative segment in front of positive segment.
        '''

        # Populate ist of starts and stops for negative segments
        boundary_list = list()

        for segment in self.elevated_segs:
            boundary_list.append(segment['start']-1)
            boundary_list.append(segment['end'])

        # Include beginning and end
        if 0 not in boundary_list:
            boundary_list.append(0)

        if len(self.counts) not in boundary_list:
            boundary_list.append(len(self.counts))

        # Sort
        boundary_list.sort()

        for x in range(0,len(boundary_list),2):
            self.normal_segs.append({'start': boundary_list[x],'end': boundary_list[x+1]})

        # Calculate counts distribution
        current_seg_counts = {0:0, 1:0, 2:0, 3:0}

        for seg in self.normal_segs:
            for i in range(seg['start'],seg['end']):
                current_count = int(self.counts[i]['count'])
                if current_count >= 3:
                    current_count = 3
                current_seg_counts[current_count] += 1

        self.normal_stats['total'] = len(self.normal_segs)

        # Update stats for normal segments
        for key, value in current_seg_counts.iteritems():
            self.normal_stats[key] += value

    def findSegs(self):
        '''
        Finds D segments in a series of positions with scores.

        Args:
            counts (list): List of dicts; chromosome, position, and number of counts starting at that position.
            score_dict (dict): Dictionary of reference scores. Describes how to handle each number.
                        Format is:
                        {(Count: (score, str('=', '>', '<', '>=', '<='))}
            D (float): Dropoff threshold
            S (float): Score threshod

        Returns:
            segment (dict): Dictionary containing segment information.
        '''

        cumul = 0
        max_score = 0
        start = 0
        end = 0
        current_count_list = dict()
        current_seg_counts = {0:0, 1:0, 2:0, 3:0}

        # Iterate through each base
        for i in range(len(self.counts)):
            # Get score and read_starts at that site
            # Update cumulative score
            score, count = self.scorePos(i)
            cumul += score

            # Update Max, if needed
            if cumul >= max_score:
                max_score = cumul
                end = i

            # Update count list for this sigmenet
            current_count_list[i] = count

            # End segment
            if cumul <= 0 or cumul <= (max_score + self.D) or i == len(self.counts):
                current_seg = {'start': start + 1,
                                'end': end + 1,
                                'score': "%.2f" % np.float64(max_score)}

                # Calculate counts distribution
                for x in range(start,end+1):
                    current_count = current_count_list[x]
                    if current_count >= 3:
                        current_count = 3
                    current_seg_counts[current_count] += 1

                if max_score >= self.S:
                    self.elevated_segs.append(current_seg)
                    self.elevated_stats['total'] += 1

                    # Update stats for elevated highSegments
                    for key, value in current_seg_counts.iteritems():
                        self.elevated_stats[key] += value

                cumul = 0
                max_score = 0
                start = i + 1
                end = i + 1
                current_count_list = dict()
                current_seg_counts = {0:0, 1:0, 2:0, 3:0}

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

    # Read in counts and score scheme.
    solution = highSegments()
    solution.setInFiles(args.counts,args.score)
    solution.readCounts()
    solution.readScores()

    # Find all D-segments
    solution.findSegs()
    print('Elevated Segs:')
    for x in solution.elevated_segs:
        print(x['start'], x['end'], x['score'])

    solution.makeNegatives()

    print('Normal Stats: ', solution.normal_stats)
    print('Elevated Stats: ', solution.elevated_stats)


    # End timer
    end_time = time.time()
    runtime = end_time-start_time

    # Write output
    writeOutput(args.out, runtime, solution.elevated_segs, solution.elevated_stats, solution.normal_stats)

if __name__ == '__main__':
    main()
