"""
convolution_cyclopeptide_sequencing.py file containing the implementation of the algorithms for the 'Part 1 - Programming' section

Part 1 of HW 5 for CISC 471, Computational Biology.

By: Hershil Devnani (20001045)

There are two ways to run the program, outlined below. Both methods run the unittests.
Sample Usage:
  $ python -m unittest unittests.py
  $ python -m main main.py

You can also modify the code in this file, main.py, in the main function to try different data sets and
parameters for each function by uncommenting and modifying the commented out function calls.
"""

from parsers import *
from helpers import *
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
    # print(theoretical_spectrum)
    # print(spectrum)
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
    # print(theoretical_spectrum)
    theoretical_spectrum = Counter(theoretical_spectrum)
    spec_temp = Counter(spec)
    return sum((theoretical_spectrum & spec_temp).values()) - sum((theoretical_spectrum - spec_temp).values())


def score_mass_seq_new(peptide_mass_seq, spectrum):
    theoretical_spectrum = cyclo_spectrum_mass_seq(peptide_mass_seq)
    # print(theoretical_spectrum)
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
    
    # trimmed_set = list()
    #
    # for i in temp:
    #     trimmed_set.append(i[0])
    #     if len(trimmed_set) == m:
    #         break
    #
    # last_score = trimmed_set[-1][1]
    # for i in temp[len(trimmed_set):]:
    #     if i[1] == last_score:
    #         trimmed_set.append(i[0])
    #
    # for i in range(len(trimmed_set)):
    #     trimmed_set[i] = formatted_to_mass_seq(trimmed_set[i])

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


def trim_new_new(leaderboard, spectrum, n):
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
        # score = temp_score(seq, spectrum)
        # score = cyclo_spectrum_mass_seq(seq)
        # mass_seq_score = score_mass_seq(seq, spectrum)
        # score_set[amino_acid_str] = str_score
        score_set[format_peptide_mass_seq(seq)] = score

    score_tuple_set = [(k, v) for k, v in score_set.items()]
    score_tuple_set.sort(key=lambda x: x[1], reverse=True)
    temp = deepcopy(score_tuple_set)

    trimmed_set = list()
    last_score = 0

    for i in temp:
        trimmed_set.append(i[0])
        last_score = i[1]
        if len(trimmed_set) == n:
            break

    # last_score = trimmed_set[-1]
    for i in temp[len(trimmed_set):]:
        if i[1] == last_score:
            trimmed_set.append(i[0])

    for i in range(len(trimmed_set)):
        trimmed_set[i] = formatted_to_mass_seq(trimmed_set[i])

    return trimmed_set


# Function for leaderboard
def leaderboard_cyclopeptide_sequencing(expands, spectrum, n):
    leaderboard = [[]]
    leader_peptide = [0]

    print(spectrum)

    # base case
    if n == 0 or spectrum == [0]:
        return [0]

    while leaderboard:
    # for i in range(3):
        leaderboard = expand(leaderboard, expands)
        leaderboard_iter = deepcopy(leaderboard)
        sortd = sorted(leaderboard)
        # print("------------------------------------------------------------")
        # print([99, 71, 137, 57, 72, 57] in leaderboard)
        # print(leaderboard)
        # print(f"expanded: {leaderboard}")
        print(f"len after expand: {len(leaderboard)}")
        highest = list()
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
                # print(peptide_mass_seq)
                # print(linear_score_mass_seq(peptide_mass_seq, spectrum))
                # if linear_score_mass_seq(peptide_mass_seq, spectrum) > linear_score_mass_seq(leader_peptide, spectrum):
                    # if len(peptide_mass_seq) <= 8:
                    #     highest.append(peptide_mass_seq)
                    #     print(f"MASS SEQ NEW: {peptide_mass_seq}")
                    # highest.append((peptide_mass_seq, linear_score_mass_seq(peptide_mass_seq, spectrum)))
                    # leader_peptide = peptide_mass_seq
                # if temp_score(peptide_mass_seq, spectrum) > temp_score(leader_peptide, spectrum):
                #     leader_peptide = peptide_mass_seq
                if score_mass_seq_new(peptide_mass_seq, spectrum) > score_mass_seq_new(leader_peptide, spectrum):
                    leader_peptide = peptide_mass_seq
                    # print(f"NEXT LEADER: {peptide_mass_seq}")
            elif peptide_mass(peptide_mass_seq) > parent_mass(spectrum):
                leaderboard.remove(peptide_mass_seq)
                # print(leaderboard)
                # print("REMOVED")
                # leaderboard = []
                # break
        leaderboard_old = trim_new(leaderboard, spectrum, n)
        leaderboard = trim_new_new(leaderboard, spectrum, n)
        sortd = sorted(leaderboard)
        print(leaderboard)
        print(f"len after trim: {len(leaderboard)}")
        print(f"len after trim old: {len(leaderboard_old)}")
        print()
        # print(f"TRIMMED: {leaderboard}")
        # print(leader_peptide)
        # break
    # print(leaderboard[:5])
    # return leader_peptide
    print(f"Highest: {highest}")
    return leader_peptide


# TODO: optimize scoring functions by allowing options for seq type and linear or cyclic


if __name__ == '__main__':

    # filename = "sample_data.txt"
    # # filename = "rosalind_ba4i_2.txt"
    # m, n, spectrum = parse_data(filename)
    #
    # scr = linear_score_mass_seq([71, 99, 129], spectrum)
    # scr = linear_score_mass_seq([113, 194, 186], spectrum)
    # print(scr)
    # scr = linear_score_mass_seq([99, 71, 137], spectrum)
    # print(scr)
    # # scr = linear_score("NQEL", spectrum)
    # # print(scr)


    filename = "../sample_data.txt"
    # filename = "rosalind_ba4i_test.txt"
    # filename = "rosalind_ba4i_test_4.txt"
    # # 99-71-137-57-72-57
    # # filename = "rosalind_ba4i_5.txt"
    # # filename = "rosalind_ba4i_5_dataset.txt"
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
    # leader = leaderboard_cyclopeptide_sequencing(expands, spectrum, n)
    leader = leaderboard_cyclopeptide_sequencing(expands, spectrum, n)
    print(f"leader: {leader}")
    formatted = format_peptide_mass_seq(leader)
    print(formatted)

    x = [[71, 99, 129, 57, 79, 58], [99, 71, 137, 57, 58, 71], [71, 99, 129, 57, 57, 80], [71, 99, 129, 57, 80, 57],
         [99, 71, 137, 57, 57, 72], [99, 71, 137, 57, 72, 57], [57, 129, 99, 71, 66, 71], [57, 137, 71, 99, 58, 71],
         [57, 129, 99, 71, 57, 80], [57, 129, 99, 71, 65, 72], [57, 137, 71, 99, 57, 72], [71, 99, 129, 57, 66, 71],
         [57, 129, 99, 71, 80, 57], [57, 137, 71, 99, 72, 57], [71, 99, 129, 57, 58, 79], [71, 99, 129, 57, 65, 72],
         [71, 99, 129, 57, 71, 66], [71, 99, 129, 57, 72, 65], [99, 71, 137, 57, 71, 58], [57, 129, 99, 71, 66, 57],
         [57, 137, 71, 99, 58, 57], [57, 129, 99, 71, 58, 79], [57, 129, 99, 71, 71, 66], [57, 129, 99, 71, 72, 65],
         [57, 129, 99, 71, 79, 58], [57, 137, 71, 99, 71, 58], [71, 99, 129, 57, 79, 57], [99, 71, 137, 57, 58, 57],
         [57, 129, 99, 71, 57, 66], [57, 129, 99, 71, 65, 58], [57, 137, 71, 99, 57, 58], [57, 129, 99, 71, 66, 58],
         [57, 129, 99, 71, 66, 65], [57, 129, 99, 71, 66, 66], [57, 137, 71, 99, 58, 58], [57, 137, 71, 99, 58, 65],
         [57, 137, 71, 99, 58, 66], [99, 71, 137, 57, 58, 58], [99, 71, 137, 57, 58, 65], [99, 71, 137, 57, 58, 66],
         [57, 129, 99, 71, 65, 57], [71, 99, 129, 57, 66, 57], [57, 129, 99, 71, 58, 65], [57, 129, 99, 71, 57, 57],
         [57, 129, 99, 71, 57, 58], [57, 129, 99, 71, 57, 65], [57, 129, 99, 71, 57, 71], [57, 129, 99, 71, 57, 72],
         [57, 129, 99, 71, 57, 79], [57, 129, 99, 71, 65, 65], [57, 129, 99, 71, 65, 66], [57, 129, 99, 71, 65, 71],
         [57, 137, 71, 99, 57, 57], [57, 137, 71, 99, 57, 65], [57, 137, 71, 99, 57, 66], [57, 137, 71, 99, 57, 71],
         [71, 99, 129, 57, 57, 57], [71, 99, 129, 57, 57, 58], [71, 99, 129, 57, 57, 65], [71, 99, 129, 57, 57, 66],
         [71, 99, 129, 57, 57, 71], [71, 99, 129, 57, 57, 72], [71, 99, 129, 57, 57, 79], [71, 99, 129, 57, 66, 58],
         [71, 99, 129, 57, 66, 65], [71, 99, 129, 57, 66, 66], [99, 71, 137, 57, 57, 57], [99, 71, 137, 57, 57, 58],
         [99, 71, 137, 57, 57, 65], [99, 71, 137, 57, 57, 66], [99, 71, 137, 57, 57, 71], [57, 129, 99, 71, 58, 57],
         [57, 129, 99, 71, 71, 57], [57, 129, 99, 71, 72, 57], [57, 129, 99, 71, 79, 57], [57, 137, 71, 99, 65, 57],
         [57, 137, 71, 99, 66, 57], [57, 137, 71, 99, 71, 57], [71, 99, 129, 57, 58, 57], [71, 99, 129, 57, 65, 57],
         [71, 99, 129, 57, 71, 57], [71, 99, 129, 57, 72, 57], [99, 71, 137, 57, 65, 57], [99, 71, 137, 57, 66, 57],
         [99, 71, 137, 57, 71, 57]]

    temp = list()
    for i in x:
        temp.append((i, linear_score_mass_seq(i, spectrum)))

    print(temp)

    # 71 - 99 - 129 - 57 - 137
    # 99-71-137-57-72-57

    a = [71, 99, 129, 57, 79, 58]
    # a = [1, 2, 3]
    b = [99, 71, 137, 57, 72, 57]

    a_theo_spec = linear_spectrum_mass_seq(a)
    b_theo_spec = linear_spectrum_mass_seq(b)
    print(f"Linear Speac: {a_theo_spec}")
    a_theo_spec = cyclo_spectrum_mass_seq(a)
    print(f"Cyclo Speac: {a_theo_spec}")

    print(b_theo_spec)
    print(spectrum)

    a_score = linear_score_mass_seq(a, spectrum)
    b_score = linear_score_mass_seq(b, spectrum)

    print(a_score)
    print(b_score)

    a_score = temp_score(a, spectrum)
    b_score = temp_score(b, spectrum)

    print(a_score)
    print(b_score)

    # compare temp score and linear score mass seq
    peptide_sample = [99, 71, 137, 57, 72, 57]
    spectrum_sample = [57, 57, 71, 99, 129, 137, 170, 186, 194, 208, 228, 265, 285, 299, 307, 323, 356, 364, 394, 422, 493]
    print(spectrum_sample)
    print(temp_score(peptide_sample, spectrum_sample))
    print(linear_score_mass_seq(peptide_sample, spectrum_sample))
    # print(score(peptide_to_mass_seq(peptide_sample), spectrum_sample))

    # sums = sum(theo_spec) + sum(spectrum_sample)
    # diffs = sum(theo_spec) - sum(spectrum_sample)
    # print(sums)
    # print(diffs)
    # for i in theo_spec:
    #     for k in spectrum_sample:
    #         sums.append(i + k)
    #
    # for i in theo_spec:
    #     for k in spectrum_sample:
    #         diffs.append(i - k)
    # print(sum(sums) - sum(diffs))

    # sample = [[128, 137, 156, 168, 113, 186, 114, 163, 128, 114, 113],
    #           [128, 137, 156, 97, 71, 113, 186, 114, 163, 128, 114, 113],
    #           [114, 163, 128, 114, 113, 128, 137, 156, 97, 71, 113, 73, 113],
    #           [113, 71, 97, 156, 137, 128, 113, 114, 128, 58, 105, 114, 128, 58],
    #           [113, 71, 97, 156, 137, 128, 113, 114, 128, 58, 105, 114, 57, 71, 58]]
    #
    # for i in sample:
    #     print(i)
    #     print(linear_score_mass_seq(i, spectrum))
# 71-99-129-57-79-58
# 71-99-129-57-79-58
# 57-131-71-186-57-186-87-101-129-57-80-71-156-87-131