#!/usr/bin/python

"""
Template for dynamic programming assignment.
The code in the template is compatible with both Python 2 and Python 3
When you finish this code, it should be at least compatible with Python 3.
"""

# Packages for commandline options:
import argparse
import sys
import pickle as pk

# Built-in exchange matrices.
with open('exchange_matrices/identity.pkl', 'rb') as f:
    identity = pk.load(f)

with open('exchange_matrices/pam250.pkl', 'rb') as f:
    pam250 = pk.load(f)

with open('exchange_matrices/blosum62.pkl', 'rb') as f:
    blosum62 = pk.load(f)


def get_args():
    """Collect the inputs."""
    parser = argparse.ArgumentParser(
        prog='PROG',
        usage='%(prog)s [options]',
        description='Aligning two sequences',
        epilog='The code was co-opted from Anton Feenstra\'s and'
        'modified by Cico Zhang'
    )
    parser.add_argument('-f', '--fasta', dest='fasta', metavar='FILE',
                        required=True, help='input alignment file (fasta)')
    parser.add_argument('-e,', '--exchange_matrix', dest='exchange_matrix',
                        metavar='EXCHANGE MATRIX NAME', help='Exchange '
                        'matrix: pam250, blosum62 or identity',
                        default='pam250')
    parser.add_argument('-l', '--local', dest='align_local',
                        action='store_true', help='Local alignment',
                        default=False)
    parser.add_argument('-g', '--global', dest='align_global',
                        action='store_true', help='Global alignment',
                        default=False)
    parser.add_argument('-s', '--semi_global', dest='align_semiglobal',
                        action='store_true', help='Semi-global alignment',
                        default=False)
    parser.add_argument('-p', '--penalty', dest='gap_penalty', type=int,
                        help='Gap penalty', default=2)
    parser.add_argument('-o', '--output', dest='alignment', required=True,
                        metavar='FILE', default='output.align',
                        help='The file to store the alignment')
    parser.add_argument('-m', '--score_matrix', dest='score_matrix',
                        required=True, metavar='FILE', default='output.align',
                        help='The file to store the score matrix')
    parser.add_argument('-v', dest='print_on_screen', action='store_true',
                        help='Print the output (alignment(s) and score '
                        'matrix) on the screen', default=False)

    args = parser.parse_args()

    if args.fasta is None:
        sys.exit('Error: no input file (fasta)')

    if not (args.align_local or args.align_global or args.align_semiglobal):
        sys.exit('Error: No alignment strategy is given: global, local or '
                 'semi-global')
    if args.align_local + args.align_global + args.align_semiglobal > 1:
        sys.exit('Error: More than one alignment strategy is given.')

    if args.exchange_matrix not in ['pam250', 'blosum62', 'identity']:
        sys.exit('Unknown exchange matrix ' + args.exchange_matrix)

    return args


class Sequence:
    """Stores a sequence object."""

    def __init__(self, Label="", Sequence=""):
        """Initialize a new Sequence object.

        Label -- identifier of sequence (text)
        Sequence -- sequence string in single-letter alphabet
        """
        self.Label = Label
        self.Sequence = Sequence

    # this makes that you can do 'print sequence' and get nice output:
    def __str__(self):
        """Return string representation of a Sequence object."""
        # newline-delimited values of all the attributes
        return ">%s\n%s" % (self.Label, self.Sequence)


def readSequences(lines):
    """Return Sequences object.

    lines -- list of lines or any object that behaves like it

    This routine parses a fasta file and returns a list of Sequence objects
    containing the sequences with label and sequence data set
    """
    seqs = []
    label = None
    seq_lines = []
    for line in lines:
        line = line.strip()      # strip off white space
        if not line:             # skip empty lines
            continue
        if line.startswith(';'):  # ignore comment lines
            continue
        # check for start of next sequence:
        if line.startswith('>'):  # label line
            # first, store the previous sequence if we had one:
            if seq_lines:
                seqs.append(Sequence(label, ''.join(seq_lines)))
                seq_lines = []
            # get the label (name) for the next sequence
            label = line[1:].strip()
        else:
            # collect all lines with sequence information for this sequence:
            seq_lines.append(line)
    # take care of the last sequence in the file
    seqs.append(Sequence(label, ''.join(seq_lines)))
    return seqs


def do_global_alignment(sequences, matrix, penalty):
    """Do pairwise global alignment using DP."""
    sequence_x = sequences[0].Sequence
    sequence_y = sequences[1].Sequence
    dimension_of_column = len(sequence_x) + 2
    dimension_of_row = len(sequence_y) + 2

    # Initialize the alignment score matrix with zeros.
    score_matrix = [[0] * (dimension_of_column) for i in range(dimension_of_row)]
    score_matrix[0][0], score_matrix[0][1], score_matrix[1][0] = ' ', '-', '-'
    for i in range(2, dimension_of_column):
        score_matrix[0][i] = sequence_x[i - 2]
    for j in range(2, dimension_of_row):
        score_matrix[j][0] = sequence_y[j - 2]

    # Fill in the score matrix
    for i in range(1, dimension_of_row):
        score_matrix[i][1] = -penalty * (i - 1)
    for j in range(1, dimension_of_column):
        score_matrix[1][j] = -penalty * (j - 1)
    for i in range(2, dimension_of_row):
        for j in range(2, dimension_of_column):
            gap_horizontal = score_matrix[i][j - 1] - penalty
            gap_vertical = score_matrix[i - 1][j] - penalty
            diagonal = score_matrix[i - 1][j - 1] + matrix[ord(sequence_y[i - 2]) - ord('A')][
                ord(sequence_x[j - 2]) - ord('A')]
            score_matrix[i][j] = max(gap_horizontal, gap_vertical, diagonal)

    # Traceback
    seq1, seq2, match_line = [], [], []
    i, j = dimension_of_row - 1, dimension_of_column - 1
    score = score_matrix[i][j]

    while i > 1 and j > 1:
        score_current = score_matrix[i][j]
        diagonal = score_matrix[i - 1][j - 1] + matrix[ord(sequence_y[i - 2]) - ord('A')][
                ord(sequence_x[j - 2]) - ord('A')]
        gap_vertical = score_matrix[i - 1][j] - penalty
        gap_horizontal = score_matrix[i][j - 1] - penalty

        if score_current == gap_vertical:
            seq1.append(sequence_y[i - 2])
            seq2.append('-')
            match_line.append(' ')
            i -= 1
        elif score_current == diagonal:
            seq1.append(sequence_y[i - 2])
            seq2.append(sequence_x[j - 2])
            if sequence_y[i - 2] == sequence_x[j - 2]:
                match_line.append('|')
            else:
                match_line.append(' ')
            i, j = i - 1, j - 1
        elif score_current == gap_horizontal:
            seq1.append('-')
            seq2.append(sequence_x[j - 2])
            match_line.append(' ')
            j -= 1

    end_gap_in_x, end_gap_in_y, match_line1 = [], [], []
    if i == 1:
        for k in range(0, j - 1):
            end_gap_in_y.append('-')
            end_gap_in_x.append(sequence_x[k])
            match_line1.append(' ')
    elif j == 1:
        for k in range(0, i - 1):
            end_gap_in_x.append('-')
            end_gap_in_y.append(sequence_y[k])
            match_line1.append(' ')

    seq1 = end_gap_in_y + seq1[::-1]
    seq2 = end_gap_in_x + seq2[::-1]
    match_line = match_line1 + match_line[::-1]
    alignment_matrix = [seq2, match_line, seq1]

    print( "score = %d" % score )
    return alignment_matrix, score_matrix

def do_local_alignment(sequences, matrix, penalty):
    """Do pairwise local alignment using DP."""
    sequence_x = sequences[0].Sequence
    sequence_y = sequences[1].Sequence
    dimension_of_column = len( sequence_x ) + 2
    dimension_of_row = len( sequence_y ) + 2
    best = 0
    optloc = (0, 0)

    # Initialize the alignment score matrix with zeros.
    score_matrix = [[0] * (dimension_of_column) for i in range( dimension_of_row )]
    score_matrix[0][0], score_matrix[0][1], score_matrix[1][0] = ' ', '-', '-'
    for i in range( 2, dimension_of_column ):
        score_matrix[0][i] = sequence_x[i - 2]
    for j in range( 2, dimension_of_row ):
        score_matrix[j][0] = sequence_y[j - 2]

    # Fill in the score matrix with the right order
    for i in range(2, dimension_of_row):
        for j in range(2, dimension_of_column):
            gap_horizontal = score_matrix[i][j - 1] - penalty
            gap_vertical = score_matrix[i - 1][j] - penalty
            diagonal = score_matrix[i - 1][j - 1] + matrix[ord(sequence_y[i - 2]) - ord('A')][
                ord(sequence_x[j - 2]) - ord('A')]
            score_matrix[i][j] = max( gap_horizontal, gap_vertical, diagonal, 0 )

            # Trace the cell with the largest score
            if score_matrix[i][j] >= best:
                best = score_matrix[i][j]
                optloc = (i, j)

    # Traceback
    i, j = optloc[0], optloc[1]
    match_line = []
    seq1, seq2 = [], []
    while i > 1 and j > 1:
        score_current = score_matrix[i][j]
        diagonal = score_matrix[i - 1][j - 1] + matrix[ord(sequence_y[i - 2]) - ord('A')][
                ord(sequence_x[j - 2]) - ord('A')]
        gap_vertical = score_matrix[i - 1][j] - penalty
        gap_horizontal = score_matrix[i][j - 1] - penalty

        if score_current == gap_vertical:
            seq1.append(sequence_y[i - 2])
            seq2.append('-')
            match_line.append(' ')
            i -= 1
        elif score_current == diagonal:
            seq1.append(sequence_y[i - 2])
            seq2.append(sequence_x[j - 2])
            if sequence_y[i - 2] == sequence_x[j - 2]:
                match_line.append('|')
            else:
                match_line.append(' ')
            i, j = i - 1, j - 1
        elif score_current == gap_horizontal:
            seq1.append('-')
            seq2.append(sequence_x[j - 2])
            match_line.append(' ')
            j -= 1
        elif score_current == 0:
            break

    seq1 = seq1[::-1]
    seq2 = seq2[::-1]
    match_line = match_line[::-1]
    alignment_matrix = [seq2, match_line, seq1]

    print("The best score is: %d" % best)
    return alignment_matrix, score_matrix

def do_semiglobal_alignment(sequences, matrix, penalty):
    """Do pairwise semi-global alignment using DP."""
    sequence_x = sequences[0].Sequence
    sequence_y = sequences[1].Sequence
    dimension_of_column = len(sequence_x) + 2
    dimension_of_row = len(sequence_y) + 2
    best = 0
    optloc = (0, 0)

    # Initialize the alignment score matrix with zeros.
    score_matrix = [[0] * (dimension_of_column) for i in range( dimension_of_row )]
    score_matrix[0][0], score_matrix[0][1], score_matrix[1][0] = ' ', '-', '-'
    for i in range( 2, dimension_of_column ):
        score_matrix[0][i] = sequence_x[i - 2]
    for j in range( 2, dimension_of_row ):
        score_matrix[j][0] = sequence_y[j - 2]

    # Fill in the score matrix
    for i in range(2, dimension_of_row):
        for j in range(2, dimension_of_column):
            gap_horizontal = score_matrix[i][j - 1] - penalty
            gap_vertical = score_matrix[i - 1][j] - penalty
            diagonal = score_matrix[i - 1][j - 1] + matrix[ord(sequence_y[i - 2]) - ord('A')][
                ord(sequence_x[j - 2]) - ord('A')]
            score_matrix[i][j] = max(gap_horizontal, gap_vertical, diagonal)

    # Trace the cell with the best score
    i = dimension_of_row - 1
    j = dimension_of_column - 1
    while i > 1:
        if score_matrix[i][dimension_of_column - 1] >= best:
            best = score_matrix[i][dimension_of_column - 1]
            optloc = (i, dimension_of_column - 1)
        i -= 1

    while j > 1:
        if score_matrix[dimension_of_row - 1][j] > best:
            best = score_matrix[dimension_of_row - 1][j]
            optloc = (dimension_of_row - 1, j)
        elif score_matrix[dimension_of_row - 1][j] == best:
            if dimension_of_column - 1 - j < dimension_of_row - 1 - optloc[0]:
                best = score_matrix[dimension_of_row - 1][j]
                optloc = (dimension_of_row - 1, j)
        j -= 1
 
    i, j = optloc[0], optloc[1]
    match_line = []
    seq1, seq2 = [], []

    while i > 1 and j > 1:
        score_current = score_matrix[i][j]

        if score_current == score_matrix[i - 1][j] - penalty:
            seq1.append(sequence_y[i - 2])
            seq2.append('-')
            match_line.append(' ')
            i -= 1
        elif score_current == score_matrix[i - 1][j - 1] + matrix[ord(sequence_y[i - 2]) - ord('A')][
            ord(sequence_x[j - 2]) - ord('A')]:
            seq1.append( sequence_y[i - 2] )
            seq2.append( sequence_x[j - 2] )
            if sequence_y[i - 2] == sequence_x[j - 2]:
                match_line.append('|')
            else:
                match_line.append(' ')
            i, j = i - 1, j - 1
        elif score_current == score_matrix[i][j - 1] - penalty:
            seq1.append('-')
            seq2.append( sequence_x[j - 2] )
            match_line.append(' ')
            j -= 1

    # Gaps in N or 5' terminal
    end_gap_in_x, end_gap_in_y, match_line1 = [], [], []
    if i == 1:
        for k in range(0, j-1):
            end_gap_in_y.append('-')
            end_gap_in_x.append(sequence_x[k])
            match_line1.append(' ')
    elif j == 1:
        for k in range(0, i-1):
            end_gap_in_x.append('-')
            end_gap_in_y.append(sequence_y[k])
            match_line1.append(' ')

    seq1 = end_gap_in_y + seq1[::-1]
    seq2 = end_gap_in_x + seq2[::-1]
    match_line = match_line1 + match_line[::-1]

    # Gaps in C or 3' terminal
    if optloc[1] == dimension_of_column - 1:
        for i in range(optloc[0] + 1, dimension_of_row):
            seq2.append('-')
            seq1.append(sequence_y[i - 2])
            match_line.append(' ')
    elif optloc[0] == dimension_of_row - 1:
        for i in range(optloc[1] + 1, dimension_of_column):
            seq1.append('-')
            seq2.append(sequence_x[i - 2])
            match_line.append(' ')
    alignment_matrix = [seq2, match_line, seq1]

    print("The best score is: %d" % best)
    return alignment_matrix, score_matrix

def print_matrix_to_file(matrix, fileName):
    """Write a matrix into file.

    matrix: a list of list in Python, storing an alignment or a score
    matrix.
    fileName: str, a file name (with a path) to store the matrix.
    It is not recommended to tinker with this function.
    """
    with open(fileName, 'w') as f:
        for row in matrix:
            print(''.join(map(str,row)), file=f)


def print_matrix_on_screen(matrix, width=5):
    """Print a matrix on the screen.

    matrix: a list of list in Python, storing an alignment or a score
    matrix.
    width: that of the space one cell occupies.
    This will facilitate your testing.
    """
    for row in matrix:
        print(''.join(['{0:>{w}}'.format(item, w=width) for item in row]))

def main():
    """Main function.

    Please change it accordingly to make the program work.
    """
    # get command line options
    args = get_args()

    # set substitution matrix:
    if args.exchange_matrix == "pam250":
        exchangeMatrix = pam250
    elif args.exchange_matrix == "blosum62":
        exchangeMatrix = blosum62
    else:
        exchangeMatrix = identity

    # read sequences from fasta file, and catch error reading file
    try:
        sequences = readSequences(open(args.fasta))
    except OSError as e:
        print("ERROR: cannot open or read fasta input file:", e.filename)

    for seq in sequences:
        print(seq)

    # call alignment routine(s):
    if args.align_global:
        alignment, score_matrix = do_global_alignment(
                sequences, exchangeMatrix, args.gap_penalty)
    elif args.align_local:
        alignment, score_matrix = do_local_alignment(
                sequences, exchangeMatrix, args.gap_penalty)
    elif args.align_semiglobal:
        alignment, score_matrix = do_semiglobal_alignment(
                sequences, exchangeMatrix, args.gap_penalty)
    else:
        sys.exit("BUG! this should not happen.")

    if args.alignment:
        print_matrix_to_file(alignment, args.alignment)
    if args.score_matrix:
        print_matrix_to_file(score_matrix, args.score_matrix)
    if args.print_on_screen:
        print_matrix_on_screen(alignment, width=5)
        print('\n')
        print_matrix_on_screen(score_matrix, width=5)


if __name__ == "__main__":
    main()


# last line
