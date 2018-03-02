
# coding: utf-8

# ## 540-Week 1
#
# This script aims to read in two fasta files, and find the longest string (or strings) that match.
#
# They can be in either forward or reverse orientation.
#
# Done:
#
#     -Load in fasta files
#     -Write Reverse Compliment Function
#     -Build suffix array
#     -Sort array
#     -Implement numpy arrays
#     -Check that a suffix array works the way you think it does
#     -Write output file
#     -Match substrings from suffix table
#     -Figure out a scheme to track what fasta file a suffix comes from.
#     -Write code to count bases for each sequence and report
#     -Store results in dictionary that can handle more than one result
#
# To Do:
#
#     -Format Output
#

# Imports

import argparse
import time
from numba import jit

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
    parser.add_argument('--out_prefix', help="Prefix for output file",
            required=False)
    return parser


def get_sequence(in_file,file_num,out_prefix):
    # Open input file
    genome = open(in_file, 'r')

    seq_list = []

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
    a_content = seq.count('A')
    c_content = seq.count('C')
    g_content = seq.count('G')
    t_content = seq.count('T')
    n_content = seq.count('N')
    all_content = len(seq)

    # Write genome contents
    with open('{0}_output.txt'.format(out_prefix),'a+') as out:
        out.write('FASTA {0}: {1}\n'.format(file_num, in_file))
        out.write(carrot_line+'\n')
        out.write('A='+str(a_content)+'\n')
        out.write('C='+str(c_content)+'\n')
        out.write('G='+str(g_content)+'\n')
        out.write('T='+str(t_content)+'\n')
        out.write('N='+str(n_content)+'\n')
        out.write('*='+str(all_content)+'\n\n')

    return seq

# Get reverse compliment for each sequence

def rev_comp(input_sequence):
    intermediate_seq = input_sequence.replace('A','W').replace('T','X').replace('C','Y').replace('G','Z')
    new_seq = intermediate_seq.replace('W','T').replace('X','A').replace('Y','G').replace('Z','C')
    new_seq = new_seq[::-1]
    return new_seq

# Build suffix array function
suffix_array = [()]
def populate_suffix(input_sequence, seq_name):
    global suffix_array
    total_length = len(input_sequence)
    for x in range(total_length):
        suffix_array.append((seq_name, input_sequence[x:], x))
    suffix_array.pop(0)

# Sort the suffix array.

def sort_suffix(input_array):
    global sorted_array
    sorted_array = sorted(input_array, key=lambda tup: tup[1])

# Write output

def write_out(out_prefix, run_time):
    with open('{0}_output.txt'.format(out_prefix),'a+') as out:
        out.write('Runtime: {0}\n'.format(run_time))
        out.write('\n')

        out.write('Match length: ' + str(max(max_counter['max_length'])) + '\n')
        out.write('Number of match strings: ' + str(len(max_counter['max_length'])) + '\n')

        for x in range(len(max_counter['max_length'])):
            out.write('\nMatch string: ' + str(max_counter['string_1'][x]) + '\n')
            out.write('\n')
            out.write('Fasta: ' + max_counter['string_1_source'][x].strip('_rev_comp') + '\n')
            out.write('Position: ' + str(max_counter['string_1_position'][x]) + '\n')
            if '_rev_comp' in max_counter['string_1_source'][x]:
                out.write('Strand: reverse\n')
            else:
                out.write('Strand: forward\n')
            out.write('\n')
            out.write('Fasta: ' + max_counter['string_2_source'][x].strip('_rev_comp') + '\n')
            out.write('Position: ' + str(max_counter['string_2_position'][x]) + '\n')
            if '_rev_comp' in max_counter['string_2_source'][x]:
                out.write('Strand: reverse\n')
            else:
                out.write('Strand: forward\n')


# Find the longest string

counter = int()
max_counter = {'max_length': [0], 'string_1': [0], 'string_2': [0], 'string_1_index': [0], 'string_2_index': [0], 'string_1_source': [0], 'string_2_source': [0], 'string_1_position': [0], 'string_2_position': []}

def find_longest_string(string_1, string_2, string_1_index, string_2_index, string_1_source, string_2_source, string_1_position, string_2_position):
    global counter
    global max_counter
    for x in range(min(len(string_1),len(string_2))):
        if string_1[x] == string_2[x]:
            counter = x + 1
        else:
            break
    if counter == int(max(max_counter['max_length'])):
        max_counter['max_length'].append(counter)
        max_counter['string_1'].append(string_1[0:counter])
        max_counter['string_2'].append(string_2[0:counter])
        max_counter['string_1_index'].append(string_1_index)
        max_counter['string_2_index'].append(string_2_index)
        max_counter['string_1_source'].append(string_1_source)
        max_counter['string_2_source'].append(string_2_source)
        max_counter['string_1_position'].append(string_1_position)
        max_counter['string_2_position'].append(string_2_position)
    elif counter > int(max(max_counter['max_length'])):
        max_counter['max_length'] = []
        max_counter['string_1'] = []
        max_counter['string_2'] = []
        max_counter['string_1_index'] = []
        max_counter['string_2_index'] = []
        max_counter['string_1_source'] = []
        max_counter['string_2_source'] = []
        max_counter['string_1_position'] = []
        max_counter['string_2_position'] = []
        max_counter['max_length'].append(counter)
        max_counter['string_1'].append(string_1[0:counter])
        max_counter['string_2'].append(string_2[0:counter])
        max_counter['string_1_index'].append(string_1_index)
        max_counter['string_2_index'].append(string_2_index)
        max_counter['string_1_source'].append(string_1_source)
        max_counter['string_2_source'].append(string_2_source)
        max_counter['string_1_position'].append(string_1_position)
        max_counter['string_2_position'].append(string_2_position)
    else:
        pass

def main():
    """Main body of script."""
    start_time = time.time()

    parser = parseArgs()
    args = parser.parse_args()
    print("\nParsed the following arguments:\n\t{0}".format(
            '\n\t'.join(['{0} = {1}'.format(arg, val) for (arg, val) in
            vars(args).items()])))

    # global vars
    global suffix_array

    # Load sequences
    genome_a_file = args.f1
    genome_b_file = args.f2

    sequences = {}
    genome_file_counter = 0
    for genome_file in [genome_a_file, genome_b_file]:
        genome_file_counter += 1
        sequences[genome_file]=(get_sequence(genome_file, genome_file_counter, args.out_prefix))

    # rev compliment
    for seq_name, seq in sequences.items():
        sequences[seq_name+'_rev_comp'] = rev_comp(seq)

    # Populate the suffix array with suffices from all sequences
    for seq_name, seq in sequences.items():
        populate_suffix(seq, seq_name)

    #suffix_array = np.array(suffix_array)

    # Sort the suffix array
    sort_suffix(suffix_array)

    # find longest strings
    for x in range(len(sorted_array)):
        string_1_index = x
        if string_1_index == len(sorted_array)-1:
            string_2_index = 0
        else:
            string_2_index = x+1
        if sorted_array[string_1_index][0] in sorted_array[string_2_index][0] or sorted_array[string_2_index][0] in sorted_array[string_1_index][0]:
            pass
        else:
            find_longest_string(sorted_array[string_1_index][1],sorted_array[string_2_index][1], string_1_index, string_2_index, sorted_array[string_1_index][0], sorted_array[string_2_index][0], sorted_array[string_1_index][2], sorted_array[string_2_index][2])

    #Calculate run time
    run_time = time.time() - start_time

    # write an output files
    write_out(args.out_prefix, run_time)
if __name__== "__main__":
  main()
