"""
convolution cyclopeptide sequencing.py file containing the implementation of the algorithms for the 'Part 1 - Programming' section

Part 1 of HW 5 for CISC 471, Computational Biology.

By: Hershil Devnani (20001045)

There are two ways to run the program, outlined below. Both methods run the unittests.
Sample Usage:
  $ python -m unittest unittests.py
  $ python -m main main.py

You can also modify the code in this file, main.py, in the main function to try different data sets and
parameters for each function by uncommenting and modifying the commented out function calls.
"""
import ast
import re
from collections import Counter
from copy import deepcopy

AMINO_ACID_MASS = {
    "": 0, "G": 57, "A": 71, "S": 87, "P": 97, "V": 99, "T": 101, "C": 103, "I": 113, "L": 113, "N": 114, "D": 115,
    "K": 128, "Q": 128, "E": 129, "M": 131, "H": 137, "F": 147, "R": 156, "Y": 163, "W": 186
}

MASS_TO_AMINO_ACID = {0: '', 57: 'G', 71: 'A', 87: 'S', 97: 'P', 99: 'V', 101: 'T', 103: 'C', 113: 'L', 114: 'N',
                      115: 'D', 128: 'Q', 129: 'E', 131: 'M', 137: 'H', 147: 'F', 156: 'R', 163: 'Y', 186: 'W'}


DEFAULT_EXPAND_SET = [57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]


# Function for convolution
def spectral_convolution(experimental_spectrum):
    convolution = dict()
    experimental_spectrum_temp = deepcopy(experimental_spectrum)
    experimental_spectrum_temp.sort(reverse=True)
    reversed_spectrum = deepcopy(experimental_spectrum)
    reversed_spectrum.sort(reverse=True)
    for amino_mass in reversed_spectrum:
        experimental_spectrum_temp.remove(amino_mass)
        for amino_sub in experimental_spectrum_temp:
            diff = amino_mass - amino_sub
            if diff != 0:
                if diff in convolution:
                    convolution[diff] += 1
                else:
                    convolution[diff] = 1
            # convolution.append(amino_mass - amino_sub)

    # sort the convolution
    sorted_vals = sorted(convolution.values(), reverse=True)
    sorted_mass = list()
    sorted_convolution = sorted(convolution)
    sorted_convolution = sorted(convolution, key=convolution.get, reverse=True)
    for mass in sorted_convolution:
        for nums in range(convolution[mass]):
            sorted_mass.append(mass)
    # for key in convolution.items():
    #     if c
    # for i in sorted_vals:
    #     sorted_mass.append()

    # return convolution
    return sorted_mass, convolution


def linear_score(peptide, spectrum):
    theoretical_spectrum = linear_spectrum(peptide)
    spectrum_copy = deepcopy(spectrum)
    score = 0
    for i in theoretical_spectrum:
        if i in spectrum_copy:
            spectrum_copy.remove(i)
            score += 1
    return score


def linear_spectrum(peptide):
    theoretical_spectrum = [""]
    theoretical_spectrum_masses = list()
    for i in range(len(peptide)):
        sub = peptide[i]
        theoretical_spectrum.append(sub)
        for k in range(i + 1, len(peptide)):
            sub += peptide[k]
            theoretical_spectrum.append(sub)
    theoretical_spectrum.remove(peptide)
    theoretical_spectrum.append(peptide)
    # print(theoretical_spectrum)

    for seq in theoretical_spectrum:
        mass = 0
        for amino_acid in seq:
            mass += AMINO_ACID_MASS[amino_acid]
        theoretical_spectrum_masses.append(mass)

    return sorted(theoretical_spectrum_masses)


# Function for scoring
def score(peptide, spectrum):
    theoretical_spectrum = cyclo_spectrum(peptide)
    spectrum_copy = deepcopy(spectrum)
    score = 0
    for i in theoretical_spectrum:
        if i in spectrum_copy:
            spectrum_copy.remove(i)
            score += 1
    return score


# Function for scoring
def score_mass_seq(peptide_mass_seq, spectrum):
    # generate peptide seq theoretical spectrum
    theo_spectrum = list()
    peptide_mass_seq_copy = deepcopy(peptide_mass_seq)
    peptide_mass_seq_copy = sorted(peptide_mass_seq_copy, reverse=True)
    # peptide_mass_seq_copy.append(0)
    for i in range(len(peptide_mass_seq_copy)):
        for k in range(i + 1, len(peptide_mass_seq_copy)):
        # for k in range(len(peptide_mass_seq_copy[i+1:])):
            diff = peptide_mass_seq_copy[i] - peptide_mass_seq_copy[k]
            theo_spectrum.append(diff)

    score = 0
    # peptide_mass_seq_copy = deepcopy(peptide_mass_seq)
    spectrum_copy = deepcopy(spectrum)
    for mass in theo_spectrum:
        if mass in spectrum_copy:
            # peptide_mass_seq_copy.remove(mass)
            spectrum_copy.remove(mass)
            score += 1
    return score


def linear_score_mass_seq(peptide_mass_seq, spectrum):
    if not peptide_mass_seq:
        return 0

    theoretical_spectrum = linear_spectrum_mass_seq(peptide_mass_seq)
    spectrum_copy = deepcopy(spectrum)
    score = 1
    for i in theoretical_spectrum:
        if i in spectrum_copy:
            spectrum_copy.remove(i)
            score += 1
    return score


def temp_score(pep, spec):
    if not pep:
        return 0
    theoretical_spectrum = linear_spectrum_mass_seq(pep)
    theoretical_spectrum = Counter(theoretical_spectrum)
    spec_temp = Counter(spec)
    return sum((theoretical_spectrum & spec_temp).values()) - sum((theoretical_spectrum - spec_temp).values())


def score_mass_seq_new(peptide_mass_seq, spectrum):
    theoretical_spectrum = cyclo_spectrum_mass_seq(peptide_mass_seq)
    spectrum_copy = deepcopy(spectrum)
    score = 0
    for i in theoretical_spectrum:
        if i in spectrum_copy:
            spectrum_copy.remove(i)
            score += 1
    return score


def cyclo_spectrum_mass_seq(peptide):
    theoretical_spectrum = [0]
    theoretical_spectrum_masses = list()
    for i in range(len(peptide)):
        sub = peptide[i % len(peptide)]
        theoretical_spectrum.append(sub)
        for k in range(i + 1, i + len(peptide) - 1):
            sub += peptide[k % len(peptide)]
            theoretical_spectrum.append(sub)
    theoretical_spectrum.append(sum(peptide))
    # print(theoretical_spectrum)

    # for seq in theoretical_spectrum:
    #     mass = 0
    #     for amino_acid in seq:
    #         mass += AMINO_ACID_MASS[amino_acid]
    #     theoretical_spectrum_masses.append(mass)

    return sorted(theoretical_spectrum)


def linear_spectrum_mass_seq(peptide):
    theoretical_spectrum = [0]
    theoretical_spectrum_masses = list()
    for i in range(len(peptide)):
        sub = peptide[i]
        theoretical_spectrum.append(sub)
        for k in range(i + 1, len(peptide)):
            sub += peptide[k]
            theoretical_spectrum.append(sub)
    # theoretical_spectrum.remove(peptide)
    # theoretical_spectrum.append(peptide)
    # theoretical_spectrum = sorted(theoretical_spectrum)
    # print(theoretical_spectrum)

    # for seq in theoretical_spectrum:
    #     mass = 0
    #     for amino_acid in seq:
    #         mass += AMINO_ACID_MASS[amino_acid]
    #     theoretical_spectrum_masses.append(mass)

    return sorted(theoretical_spectrum)


def cyclo_spectrum(peptide):
    theoretical_spectrum = [""]
    theoretical_spectrum_masses = list()
    for i in range(len(peptide)):
        sub = peptide[i % len(peptide)]
        theoretical_spectrum.append(sub)
        for k in range(i + 1, i + len(peptide) - 1):
            sub += peptide[k % len(peptide)]
            theoretical_spectrum.append(sub)
    theoretical_spectrum.append(peptide)
    # print(theoretical_spectrum)

    for seq in theoretical_spectrum:
        mass = 0
        for amino_acid in seq:
            mass += AMINO_ACID_MASS[amino_acid]
        theoretical_spectrum_masses.append(mass)

    return sorted(theoretical_spectrum_masses)


def trim(leaderboard, spectrum, n):
    # build dict of scores of leaderboard
    score_set = dict()
    new_leaderboard = list()
    for seq in leaderboard:
        amino_acid_str = mass_to_peptide(seq)
        str_score = score(amino_acid_str, spectrum)
        score_set[amino_acid_str] = str_score
    # return top scores
    sorted_score_set = sorted(score_set, key=score_set.get, reverse=True)
    counter = 0
    prev_score = 0
    for peptide in sorted_score_set:
        curr_score = score_set[peptide]
        # print(peptide)
        # print(score_set[peptide])
        if curr_score != prev_score:
            counter += 1
        prev_score = score_set[peptide]
        peptide_mass_seq = list()
        for i in peptide:
            peptide_mass_seq.append(AMINO_ACID_MASS[i])
        new_leaderboard.append(peptide_mass_seq)
        if counter == n:
            break

    # add and account for ties

        # new_leaderboard.append(peptide_to_mass_seq(i))
        # if len(new_leaderboard) == n:
        #     if
        #     break

    return new_leaderboard


def trim_new(leaderboard, spectrum, n):
    if len(leaderboard) < n:
        return leaderboard
    # build dict of scores of leaderboard
    score_set = dict()
    new_leaderboard = list()
    for seq in leaderboard:
        # amino_acid_str = mass_to_peptide(seq)
        # str_score = score(amino_acid_str, spectrum)
        # score = score_mass_seq(seq, spectrum)
        score = linear_score_mass_seq(seq, spectrum)
        # mass_seq_score = score_mass_seq(seq, spectrum)
        # score_set[amino_acid_str] = str_score
        score_set[format_peptide_mass_seq(seq)] = score

    score_tuple_set = [(k, v) for k, v in score_set.items()]
    score_tuple_set.sort(key=lambda x: x[1], reverse=True)
    temp = deepcopy(score_tuple_set)

    trimmed_set = list()
    curr = tuple()
    for i in score_tuple_set:
        curr = i
        trimmed_set.append(i[0])
        temp.remove(curr)
        if len(trimmed_set) == n:
            break

    for i in temp:
        score_i = i[1]
        if score_i == curr[1]:
            trimmed_set.append(i[0])

    for i in range(len(trimmed_set)):
        trimmed_set[i] = formatted_to_mass_seq(trimmed_set[i])

    return trimmed_set

    # # return top scores
    # sorted_score_set = sorted(score_set, key=score_set.get, reverse=True)
    # counter = 0
    # prev_score = 0
    # for peptide in sorted_score_set:
    #     curr_score = score_set[peptide]
    #     # print(peptide)
    #     # print(score_set[peptide])
    #     if curr_score != prev_score:
    #         counter += 1
    #     # counter += 1
    #     prev_score = score_set[peptide]
    #     peptide_mass_seq = list()
    #     peptide_mass_seq = formatted_to_mass_seq(peptide)
    #     # for i in peptide:
    #     #     peptide_mass_seq.append(AMINO_ACID_MASS[i])
    #     # new_leaderboard.append(ast.literal_eval(peptide))
    #     new_leaderboard.append(peptide_mass_seq)
    #     if len(new_leaderboard) >= n:
    #         break
    #     # if counter == n:
    #     #     break
    #
    # # add and account for ties
    #
    #     # new_leaderboard.append(peptide_to_mass_seq(i))
    #     # if len(new_leaderboard) == n:
    #     #     if
    #     #     break
    #
    # return new_leaderboard


def mass_to_peptide(mass_list):
    peptide = ""
    for i in mass_list:
        peptide += MASS_TO_AMINO_ACID[i]
    return peptide


def peptide_to_mass_seq(peptide):
    mass_list = []
    for i in peptide:
        mass_list.append(AMINO_ACID_MASS[i])
    return peptide


# Function for leaderboard
def leaderboard_cyclopeptide_sequencing(expands, spectrum, n):
    leaderboard = [[]]
    leader_peptide = [0]

    print(spectrum)

    while leaderboard:
    # for i in range(3):
        leaderboard = expand(leaderboard, expands)
        leaderboard_iter = deepcopy(leaderboard)
        # print("------------------------------------------------------------")
        # print([99, 71, 137, 57, 72, 57] in leaderboard)
        print(leaderboard)
        # print(f"expanded: {leaderboard}")
        print(f"len after expand: {len(leaderboard)}")
        # print(len(leaderboard[0]))
        for peptide_mass_seq in leaderboard_iter:
            # print(leaderboard)
            # print(peptide_mass_seq)
            # print(peptide_mass(peptide_mass_seq))
            # print(parent_mass(spectrum))

            # mass_to_peptide(peptide)
            if peptide_mass(peptide_mass_seq) == parent_mass(spectrum):
                # if linear_score(mass_to_peptide(peptide_mass_seq), spectrum) > linear_score(mass_to_peptide(leader_peptide), spectrum):
                #     leader_peptide = peptide_mass_seq
                # if score_mass_seq(peptide_mass_seq, spectrum) > score_mass_seq(leader_peptide, spectrum):
                #     leader_peptide = peptide_mass_seq
                print(peptide_mass_seq)
                print(linear_score_mass_seq(peptide_mass_seq, spectrum))
                # if linear_score_mass_seq(peptide_mass_seq, spectrum) > linear_score_mass_seq(leader_peptide, spectrum):
                #     leader_peptide = peptide_mass_seq
                if temp_score(peptide_mass_seq, spectrum) > temp_score(leader_peptide, spectrum):
                    leader_peptide = peptide_mass_seq
                # if score_mass_seq_new(peptide_mass_seq, spectrum) > score_mass_seq_new(leader_peptide, spectrum):
                #     leader_peptide = peptide_mass_seq
                    # print(f"NEXT LEADER: {peptide_mass_seq}")
            elif peptide_mass(peptide_mass_seq) > parent_mass(spectrum):
                leaderboard.remove(peptide_mass_seq)
                # print(leaderboard)
                # print("REMOVED")
                # leaderboard = []
                # break
        leaderboard = trim_new(leaderboard, spectrum, n)
        print(f"len after trim: {len(leaderboard)}")
        # print(f"TRIMMED: {leaderboard}")
        # print(leader_peptide)
        # break
    # print(leaderboard[:5])
    # return leader_peptide
    return leader_peptide


def expand_set(convolutions, m):
    # remove_dups = list(set(convolutions))
    remove_dups = list()
    for i in convolutions:
        if i not in remove_dups:
            remove_dups.append(i)

    subset = list()
    for i in remove_dups:
        if i > 57 or i < 200:
            subset.append(i)
        if len(subset) == m:
            break

    return subset


def expand_set_with_ties(conv_tuple_set, m):


    temp = deepcopy(conv_tuple_set)
    sublist = list()
    subset = set()
    pair = tuple()

    for i in conv_tuple_set:
        pair = i
        mass = pair[0]
        if mass >= 57 and mass <= 200:
            subset.add(mass)
            temp.remove(pair)
        if len(subset) == m:
            break

    for i in temp:
        mass = i[0]
        if mass >= 57 and mass <= 200:
            if i[1] == pair[1]:
                subset.add(mass)


    # remove_dups = list(set(convolutions))
    # counter = 0
    # subset = list()
    # prev_count = 0
    #
    # for i in convolutions:
    #     subset.append(i)
    # indx = convolutions.index(subset[-1])
    # if
        # curr_count = convolutions.count(i)
        # if curr_count != prev_count:
        # # if i not in subset:
        #     counter += 1
        # if i not in subset:
        #     subset.append(i)
        # prev_count = convolutions.count(i)
        # if counter == m:
        #     break

    # remove_dups = list()
    # for i in convolutions:
    #     if i not in remove_dups:
    #         remove_dups.append(i)
    #
    # subset = list()
    # for i in remove_dups:
    #     if i > 57 or i < 200:
    #         subset.append(i)
    #     if len(subset) == m:
    #         break

    return subset


def expand(leaderboard, expands):
    new_leaderboard = list()
    for i in range(len(leaderboard)):
        for option in expands:
            new_leaderboard.append(leaderboard[i] + [option])
    return new_leaderboard


# Function for getting amino acid mass
def amino_acid_mass(amino_acid):
    return AMINO_ACID_MASS[amino_acid]


def peptide_mass(peptide):
    mass = 0
    for amino_acid in peptide:
        mass += amino_acid
    return mass


# Function for parent mass
def parent_mass(experimental_spectrum):
    return max(experimental_spectrum)


def format_peptide_mass_seq(mass_seq):
    seq_str = ""
    for i in range(len(mass_seq) - 1):
        seq_str += str(mass_seq[i]) + "-"
    seq_str += str(mass_seq[-1])
    return seq_str


def formatted_to_mass_seq(formatted_str):
    mass_seq = formatted_str.split("-")
    mass_seq = [int(i) for i in mass_seq]
    return mass_seq


def parse_conv_data(filename):
    """Read in the text file and build the graph

    :param filename: text files will graph data
    :return: graph in the form of a dictionary
    """

    # grab text file
    # parse via delimiter
    # store as dictionary, key = node, values = list of next nodes
    graph_edges = list()

    with open(filename) as file:
        data_set = file.readlines()

    for line in data_set:
        # print(line.strip('\n').split(' -> '))
        edge = re.split(' ', line.strip('\n'))  # parse line input into edge + next nodes
        # edge = [int(i) for i in edge]  # convert data of each edge to integers
        # print(re.split(' ', line.strip('\n')))
        # graph_edges[edge[0]] = edge[1:]
        edge = [int(i) for i in edge]
        return edge
    return graph_edges


def parse_scoring_data(filename):
    """Read in the text file and build the graph

    :param filename: text files will graph data
    :return: graph in the form of a dictionary
    """

    with open(filename) as file:
        data_set = file.readlines()

    # parse k and t
    peptide = re.split(' ', data_set[0].strip('\n'))
    peptide = peptide[0]

    spectrum = re.split(' ', data_set[1].strip('\n'))
    spectrum = [int(i) for i in spectrum]

    # for i in range(len(data_set)):
    #     data_set[i] = data_set[i].strip('\n')

    return peptide, spectrum


def parse_data(filename):
    """Read in the text file and build the graph

    :param filename: text files will graph data
    :return: graph in the form of a dictionary
    """

    with open(filename) as file:
        data_set = file.readlines()

    # parse k and t
    m = re.split(' ', data_set[0].strip('\n'))
    n = re.split(' ', data_set[1].strip('\n'))

    m = int(m[0])
    n = int(n[0])

    spectrum = re.split(' ', data_set[2].strip('\n'))
    spectrum = [int(i) for i in spectrum]

    # for i in range(len(data_set)):
    #     data_set[i] = data_set[i].strip('\n')

    return m, n, spectrum


if __name__ == '__main__':
    x = [0, 137, 186, 323]
    # test_conv = spectral_convolution(x)
    # print(test_conv)

    # filename = "sample_convolution_data.txt"
    # filename = "rosalind_ba4h_8.txt"
    # x = parse_conv_data(filename)
    # # print(x)
    # test_conv, _ = spectral_convolution(x)
    # # print(test_conv)
    # # print()
    # for i in test_conv:
    #     print(i, end=" ")
    # print()

    # filename = "sample_convolution_solution.txt"
    # solution = parse_data(filename)
    # solution.sort()
    # test_conv.sort()

    # print(test_conv == solution)

    # filename = "sample_convolution_data.txt"
    # conv_data = parse_conv_data(filename)
    # convolutions = spectral_convolution(conv_data)
    # expands = expand_set(convolutions, 3)
    # print(expands)

    # filename = "sample_scoring_data_2.txt"
    # filename = "rosalind_ba4f_2.txt"
    # peptide, spectrum = parse_scoring_data(filename)
    # print(peptide)
    # print(spectrum)
    # theoretical_spectrum = cyclo_spectrum(peptide)
    # print(theoretical_spectrum)
    # for i in theoretical_spectrum:
    #     print(i, end=" ")
    # print()
    # peptide_score = score(peptide, spectrum)
    # print(peptide_score)
    #
    # linear_spec = linear_spectrum("NQEL")
    # print(linear_spec)
    #
    # x = formatted_to_mass_seq("113-147-71-129")
    # print(x)

    # filename = "sample_scoring_data.txt"
    # peptide, spectrum = parse_scoring_data(filename)
    #
    # theo_new = linear_spectrum_mass_seq([114, 128, 129, 113])
    # print(theo_new)
    # scr = linear_score_mass_seq([114, 128, 129, 113], spectrum)

    theo_cyclo = cyclo_spectrum_mass_seq([113, 129, 128, 114])
    print(theo_cyclo)

    filename = "../sample_data.txt"
    # filename = "rosalind_ba4i_2.txt"
    m, n, spectrum = parse_data(filename)

    scr = linear_score_mass_seq([71, 99, 129], spectrum)
    scr = linear_score_mass_seq([113, 194, 186], spectrum)
    print(scr)
    scr = linear_score_mass_seq([99, 71, 137], spectrum)
    print(scr)
    # scr = linear_score("NQEL", spectrum)
    # print(scr)


    filename = "../sample_data.txt"
    # filename = "rosalind_ba4i_18.txt"
    m, n, spectrum = parse_data(filename)
    print(spectrum)
    convolutions, conv_dict = spectral_convolution(spectrum)

    conv_tuple_set = [(k, v) for k, v in conv_dict.items()]
    conv_tuple_set.sort(key=lambda x: x[1], reverse=True)

    print(f"Convolved: {convolutions}")
    # expands = expand_set(convolutions, 5)
    edited_convolutions = list()
    for i in convolutions:
        if i not in edited_convolutions:
            edited_convolutions.append(i)
    # for i in convolutions:
    #     if i in MASS_TO_AMINO_ACID:
    #         edited_convolutions.append(i)
    print(edited_convolutions)
    # expands = expand_set(convolutions, m)
    expands = expand_set_with_ties(conv_tuple_set, m)
    print(f"Expands: {expands}")
    leader = leaderboard_cyclopeptide_sequencing(expands, spectrum, n)
    print(f"leader: {leader}")
    formatted = format_peptide_mass_seq(leader)
    print(formatted)

    # filename = "leaderboard_sample.txt"
    # m, n, spectrum = parse_data(filename)
    # leader = leaderboard_cyclopeptide_sequencing(DEFAULT_EXPAND_SET, spectrum, n)
    # print(f"leader: {leader}")
    # formatted = format_peptide_mass_seq(leader)
    # print(formatted)

    # filename = "rosalind_ba4g_3.txt"
    # m, n, spectrum = parse_data(filename)
    # leader = leaderboard_cyclopeptide_sequencing(DEFAULT_EXPAND_SET, spectrum, n)
    # print(f"leader: {leader}")
    # formatted = format_peptide_mass_seq(leader)
    # print(formatted)
