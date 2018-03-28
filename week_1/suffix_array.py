# ## 540-Week 1
#
# This script aims to read in two fasta files, and find the longest string (or strings) that match.
#
# They can be in either forward or reverse orientation.

# Imports

from pympler import asizeof
import argparse
import time
import numpy as np

# File I/O

# Function to load sequence into

def parseArgs():
    """Returns `argparse.ArgumentParser` for script."""

    parser = argparse.ArgumentParser(
            description='Process FASTQ files for homology search'
            )
    parser.add_argument('--f1', help="FASTQ file with genome a",
            required=True)
    parser.add_argument('--f2', help="FASTQ file with genome b",
            required=True)
    parser.add_argument('--out', help="Prefix for output file",
            required=False)
    return parser

# Global args
# Reads
genomes = list()
suffix_array = None

class genome(object):
    '''
    Contains a whole genome input as f1 or f2.
    '''

    def __init__(self):
        '''
        Initilize relevant variables.
        '''

        # Basic from FASTA
        self.name = str()
        self.long_name = str()
        self.seq = str()
        self.length = int()
        self.content = dict()

def getSequence(in_file):
    current_genome = genome()
    seq_list = list()

    # Open input file
    with open(in_file, 'r') as fasta:
        for line in fasta:
            line = line.strip('\n')
            if line[0] == ">":
                current_genome.long_name = line
                line = line.strip('>')
                line = line.split()
                current_genome.name = line[0]
            else:
                seq_list.append(line.upper())

    current_genome.seq = ''.join(map(str, seq_list))

    return current_genome

def countChars(current_genome):
    '''
    Takes a genome as input, and counts the various characters.
    '''

    current_genome.length = len(current_genome.seq)
    current_genome.content['A'] = current_genome.seq.count('A')
    current_genome.content['C'] = current_genome.seq.count('C')
    current_genome.content['T'] = current_genome.seq.count('T')
    current_genome.content['G'] = current_genome.seq.count('G')
    current_genome.content['N'] = current_genome.seq.count('N')

def revComp(input_genome):
    '''
    Takes in a genome, opens its sequence, rev comps it, and returns a new genome.
    '''
    new_genome = genome()
    input_sequence = input_genome.seq

    intermediate_seq = input_sequence.replace('A','W').replace('T','X').replace('C','Y').replace('G','Z')
    new_seq = intermediate_seq.replace('W','T').replace('X','A').replace('Y','G').replace('Z','C')
    new_seq = new_seq[::-1]

    new_genome.seq = new_seq
    new_genome.name = input_genome.name+'_rev_comp'

    return new_genome

class suffix(object):
    '''
    A suffix of a sequence. Attributes include sequence, length, forward or reverse,
    '''

    def __init__(self):
        '''
        Initilize relevant variables.
        '''

        # Basic from FASTA
        self.seq = str()
        self.length = int()
        self.source = str()
        self.orientation = str() # True if forward, False if reverse.
        self.position = int()


def populateSuffix():
    '''
    Makes suffixes for each base in each sequence. Adds them to suffix_array variable.
    '''
    global genomes
    global suffix_array

    # Make array
    total_length = 0
    for g in genomes:
        total_length += g.length

    #print('Size of suffix array: {0}'.format(asizeof.asizeof(suffix_array)))
    suffix_array = np.empty(total_length, dtype=object)
    #print('Size of suffix array: {0}'.format(asizeof.asizeof(suffix_array)))

    # Make suffices and put into array
    total_count = 0
    for g in genomes:
        for i in range(g.length):
            current_suffix = suffix()
            current_suffix.length = len(g.seq[i:])
            current_suffix.source = g.name
            current_suffix.position = i
            if 'rev_comp' in g.name:
                current_suffix.orientation = 'reverse'
            else:
                current_suffix.orientation = 'forward'
            suffix_array[total_count] = current_suffix
            #print('Current suffix ({0}) size: {1}; array size: {2}'.format(total_count, asizeof.asizeof(current_suffix), asizeof.asizeof(suffix_array)))
            total_count += 1

def getSuf(suffix):
    '''
    Compares two suffix sequences
    '''

    global genomes

    seq = None

    for g in genomes:
        if suffix.source == g.name:
            seq = g.seq[suffix.position:]

    return seq

def cmpSufs(suf_a, suf_b):
    #print(suf_a.source, suf_a.position, suf_b.source, suf_b.position)

    a = getSuf(suf_a)
    b = getSuf(suf_b)

    return cmp(a, b)

def sortSuffix():
    '''
    Sorts Suffix using custom compare function
    '''
    global suffix_array

    suffix_array = sorted(suffix_array, cmp = cmpSufs)


def findLongest():
    '''
    Finds longest matching strings in suffix_array
    '''
    global suffix_array
    global max_counter
    global max_suffices

    max_counter = None
    max_suffices = list()

    for i in range(len(suffix_array)-1):
        suffix_a = suffix_array[i]
        suffix_b = suffix_array[i+1]

        #print(suffix_a.source, suffix_a.position, suffix_a.length, getSuf(suffix_a))
        #print(suffix_b.source, suffix_b.position, suffix_b.length, getSuf(suffix_b))

        check_length = min(suffix_a.length, suffix_b.length)
        counter = 0
        for x in range(check_length):
            #print(getSuf(suffix_a)[x], getSuf(suffix_b)[x])
            if getSuf(suffix_a)[x] == getSuf(suffix_b)[x]:
                counter += 1
            else:
                if counter > max_counter:
                    max_counter = counter
                    max_suffices = [(suffix_a,suffix_b,counter)]
                elif counter == max_counter:
                    max_suffices.append((suffix_a,suffix_b,counter))
                break

def writeOutput(args, runtime):
    '''
    Writes output in a specific format to text file.
    '''

    global genomes
    global suffix_array
    global max_counter
    global max_suffices

    with open('{0}_results.txt'.format(args.out), 'w') as out_file:
        # Header
        out_file.write('Assignment: GS540 HW1\n')
        out_file.write('Name: David Bacsik\n')
        out_file.write('Email: dbacsik@uw.edu\n')
        out_file.write('Language: Python\n')
        out_file.write('Run time: {0}\n\n'.format(runtime))

        #QC
        out_file.write('Fasta 1: {0}\n'.format(args.f1))
        g = genomes[0]
        out_file.write(g.long_name+'\n')
        for key, val in g.content.iteritems():
            out_file.write('{0}={1}\n'.format(key,val))
        out_file.write('*={0}\n\n'.format(g.length))


        out_file.write('Fasta 2: {0}\n'.format(args.f2))
        g = genomes[1]
        out_file.write(g.long_name+'\n')
        for key, val in g.content.iteritems():
            out_file.write('{0}={1}\n'.format(key,val))
        out_file.write('*={0}\n\n'.format(g.length))

        # results
        out_file.write('Match length: {0}\n'.format(max_counter))
        out_file.write('Number of match strings: {0}\n\n'.format(len(max_suffices)))

        for match in max_suffices:
            out_file.write('Match string: {0}\n'.format(getSuf(match[0])[:match[2]]))
            out_file.write('Description:\n\n')

            out_file.write('Fasta: {0}\n'.format(match[0].source.strip('_rev_comp')))
            out_file.write('Position: {0}\n'.format(match[0].position))
            out_file.write('Strand: {0}\n\n'.format(match[0].orientation))

            out_file.write('Fasta: {0}\n'.format(match[1].source.strip('_rev_comp')))
            out_file.write('Position: {0}\n'.format(match[1].position))
            out_file.write('Strand: {0}\n\n'.format(match[1].orientation))


# MAIN
def main():
    '''Main body of script.'''

    # Global vars
    global genomes

    # Start timer
    start_time = time.time()

    # Parse command line arguments
    parser = parseArgs()
    print("\nExecuting {0} at {1}".format(parser.prog, time.asctime()))
    args = parser.parse_args()
    print("\nParsed the following arguments:\n\t{0}".format(
            '\n\t'.join(['{0} = {1}'.format(arg, val) for (arg, val) in
            vars(args).items()])))

    # Load genomes
    print('Loading genomes.\n')
    for f in [args.f1, args.f2]:
        print(f)
        current_genome = getSequence(f)
        genomes.append(current_genome)

    # Make Reverse complements
    print('Making reverse compliments.\n')
    for g in genomes[:]:
        g1 = revComp(g)
        genomes.append(g1)

    print('Calculating sequence makeup.\n')
    # Calculate sequence makeup
    for g in genomes[:]:
        countChars(g)

    print('Populating suffix array with start positions.\n')
    # Make Suffix Array
    populateSuffix()

    print('Sorting.\n')
    # Sort Suffix Array
    sortSuffix()

    print('Finding longest match.\n')
    # Find longest matching string
    findLongest()

    # End timer
    end_time = time.time()
    runtime = end_time-start_time

    print('Writing output.\n')
    # Write output
    writeOutput(args, runtime)

    print('Done.')

if __name__ == '__main__':
    main()
