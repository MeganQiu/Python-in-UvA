import random, subprocess
from os import listdir, system
from os.path import isfile, join

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def generate_tests():
    max_seq1_len = 20
    max_seq2_len = 20
    tests_num = 10
    for i in range(0, tests_num):
        allowed_chars = "ARNDCQEGHILKMFPSTWYV"
        seq1_len = random.randint(10, max_seq1_len)
        seq2_len = random.randint(10, max_seq2_len)

        seq1 = ''.join([allowed_chars[random.randint(0, len(allowed_chars) - 1)] for _ in range(0, seq1_len)])
        seq2 = ''.join([allowed_chars[random.randint(0, len(allowed_chars) - 1)] for _ in range(0, seq2_len)])

        with open('tests\test%d.in' % i, 'w') as file:
            file.write('>seq1\n%s\n>seq2\n%s\n' % (seq1, seq2))


def generate_correct_results(input_file):
    base_name = input_file.replace('.in', '')

    for exchange_matrix in ['pam250', 'blosum62', 'identity']:
        for alignment_method in ['global', 'semi_global', 'local']:
            output_align = 'results\%s_%s_%s.align.correct' % (base_name, exchange_matrix, alignment_method)
            output_score = 'results\%s_%s_%s.score.correct' % (base_name, exchange_matrix, alignment_method)

            subprocess.call('python36 align.py -f %s -e %s --%s -p 2 -o %s -m %s ' %
                   ('tests\%s' % input_file, exchange_matrix, alignment_method, output_align, output_score))
            
            print (base_name, alignment_method, exchange_matrix)


def check_results(input_file):
    base_name = input_file.replace('.in', '')

    for exchange_matrix in ['pam250', 'blosum62', 'identity']:
        for alignment_method in ['global', 'semi_global', 'local']:
            output_align = 'results\%s_%s_%s.align' % (base_name, exchange_matrix, alignment_method)
            output_score = 'results\%s_%s_%s.score' % (base_name, exchange_matrix, alignment_method)

            subprocess.call('python36 align.py -f %s -e %s --%s -p 2 -o %s -m %s ' %
                   ('tests\%s' % input_file, exchange_matrix, alignment_method, output_align, output_score))

            if not isfile(output_align):
                exit(bcolors.BOLD + bcolors.FAIL + 'Fail: alignment file not found.\n'
                                                   'Check if your script is writing to the file' + bcolors.ENDC +
                     'pointed to by the --output (-o) argument.')

            if not isfile(output_align):
                exit(bcolors.BOLD + bcolors.FAIL + 'Fail: scoring matrix file not found.\n' + bcolors.ENDC +
                     'Check if your script is writing to the file'
                     'pointed to by the --score_matrix (-m) argument.')

            runs = [['align', 'alignment', output_align], ['score', 'scoring', output_score]]

            for extension, run_type, filename in runs:
                with open(filename, 'r') as check_file:
                    with open(filename.replace('.%s' % extension, '.%s.correct' % extension), 'r') as correct_file:
                        # Ignore leading and trailing whitespace (.strip()), empty lines, and capitalisation
                        check_raw = check_file.readlines()
                        check_lines = ''.join(check_raw)
                        check_lines_clean = ''.join([line.strip(' \t').lower() for line in check_raw
                                                     if line.strip(' \t') != '\n'])
                        correct_raw = correct_file.readlines()
                        correct_lines = ''.join(correct_raw)
                        correct_lines_clean = ''.join([line.strip(' \t').lower() for line in correct_raw
                                                       if line.strip(' \t') != '\n'])

                        if check_lines_clean != correct_lines_clean:
                            print(bcolors.BOLD + bcolors.FAIL + 'Fail' + bcolors.ENDC + bcolors.FAIL +
                                  ' - [%s] %s does not match (%s, %s)' %
                                  (base_name, run_type, exchange_matrix, alignment_method) + bcolors.ENDC)
                            print('Correct:\n' + correct_lines)
                            print('Yours:\n' + check_lines)
                        else:
                            print(bcolors.BOLD + bcolors.OKGREEN + 'Pass' + bcolors.ENDC + ' - [%s] %s (%s, %s)' %
                                  (base_name, run_type, exchange_matrix, alignment_method) + bcolors.ENDC)


if __name__ == '__main__':
    for test_file in listdir('tests'):
        if isfile(join('tests', test_file)):
            # generate_correct_results(test_file)
            check_results(test_file)
