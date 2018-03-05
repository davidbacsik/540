    parser.add_argument('--b', help='2-D Matrix of amino acid substitution scores.',
            required=False, default='blosum62.txt')
    parser.add_argument('--g', help='Gap Penalty. Should be negative number',
            required=False, default='-6')

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
