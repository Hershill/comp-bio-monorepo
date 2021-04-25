"""
convolution_cyclopeptide_sequencing.py file containing the implementation of the algorithms for the 'Part 1 - Programming' section

Part 1 of HW 5 for CISC 471, Computational Biology.

By:
    - Hershil Devnani (20001045)
    - Rayan Shaikli (20059806)

There are two ways to run the program, outlined below. Both methods run the unittests.
Sample Usage:
  $ python -m unittest convolution_cyclopeptide_sequencing_tests.py
  $ python -m main main.py

You can also modify the code in this file, main.py, in the main function to try different data sets and
parameters for each function by uncommenting and modifying the commented out function calls.
"""

from parsers import *
from helpers import *
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
    """Perform convolution of the experimental spectrum

    :param experimental_spectrum: the given experimental spectrum
    :return: list of elements in the convolution of the experimental spectrum
    """

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

    # sort the convolution
    sorted_mass = list()
    sorted_convolution = sorted(convolution, key=convolution.get, reverse=True)
    for mass in sorted_convolution:
        for nums in range(convolution[mass]):
            sorted_mass.append(mass)

    return sorted_mass, convolution


def cyclo_spectrum(peptide):
    """Generate the cyclic mass spectrum for a peptide

    :param peptide: the peptide as a string
    :return: the cyclic spectrum of the peptide
    """

    theoretical_spectrum = [""]
    theoretical_spectrum_masses = list()

    for i in range(len(peptide)):
        sub = peptide[i % len(peptide)]
        theoretical_spectrum.append(sub)
        for k in range(i + 1, i + len(peptide) - 1):
            sub += peptide[k % len(peptide)]
            theoretical_spectrum.append(sub)

    theoretical_spectrum.append(peptide)

    for seq in theoretical_spectrum:
        mass = 0
        for amino_acid in seq:
            mass += AMINO_ACID_MASS[amino_acid]
        theoretical_spectrum_masses.append(mass)

    return sorted(theoretical_spectrum_masses)


def linear_spectrum(peptide):
    """Generate the linear mass spectrum for a peptide

    :param peptide: the peptide as a string
    :return: the linear spectrum of the peptide
    """

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

    for seq in theoretical_spectrum:
        mass = 0
        for amino_acid in seq:
            mass += AMINO_ACID_MASS[amino_acid]
        theoretical_spectrum_masses.append(mass)

    return sorted(theoretical_spectrum_masses)


def cyclo_spectrum_mass_seq(peptide):
    """Return the cyclic spectrum of a given peptide

    :param peptide: the peptide as a list of amino acid masses
    :return: theoretical cyclic spectrum of the peptide
    """

    theoretical_spectrum = [0]

    for i in range(len(peptide)):
        sub = peptide[i % len(peptide)]
        theoretical_spectrum.append(sub)
        for k in range(i + 1, i + len(peptide) - 1):
            sub += peptide[k % len(peptide)]
            theoretical_spectrum.append(sub)
    theoretical_spectrum.append(sum(peptide))

    return sorted(theoretical_spectrum)


def linear_spectrum_mass_seq(peptide):
    """Return the linear spectrum of a given peptide

    :param peptide: the peptide as a list of amino acid masses
    :return: theoretical linear spectrum of the peptide
    """

    theoretical_spectrum = [0]

    for i in range(len(peptide)):
        sub = peptide[i]
        theoretical_spectrum.append(sub)
        for k in range(i + 1, len(peptide)):
            sub += peptide[k]
            theoretical_spectrum.append(sub)

    return sorted(theoretical_spectrum)


def linear_score_mass_seq(peptide_mass_seq, spectrum):
    """Score a peptide given a list of amino acid masses

    :param peptide_mass_seq: the peptide as a list of it's amino acidmasses
    :param spectrum: the spectrum to score the peptide against
    :return: the score of the peptide
    """

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


def expand(leaderboard, expands):
    """Expands the peptide in the leaderboard

    :param leaderboard: the current leaderboard
    :param expands: set of possible amino acids the peptides on the leaderboard are allowed to be expanded by
    :return: the new leaderboard
    """

    new_leaderboard = list()

    for i in range(len(leaderboard)):
        for option in expands:
            new_leaderboard.append(leaderboard[i] + [option])

    return new_leaderboard


# Function for scoring
def score_peptide(peptide, spectrum, scoring_type):
    """Score the peptide given a spectrum

    :param peptide: the peptide, as a string of amino acids or a list of amino acid masses
    :param spectrum: the spectrum to score against
    :param scoring_type: the method for scoring, either "cyclic" or "linear" scoring method options
    :return: the score of the peptide
    """

    # determine input type and select functions accordingly
    if type(peptide) == list:
        function_mapping = {
            "linear": linear_spectrum_mass_seq,
            "cyclic": cyclo_spectrum_mass_seq
        }
    elif type(peptide) == str:
        function_mapping = {
            "linear": linear_spectrum,
            "cyclic": cyclo_spectrum
        }
    else:
        return "Invalid peptide type"

    theoretical_spectrum = function_mapping[scoring_type](peptide)
    spectrum_copy = deepcopy(spectrum)
    score = 0

    for i in theoretical_spectrum:
        if i in spectrum_copy:
            spectrum_copy.remove(i)
            score += 1

    return score


def trim(leaderboard, spectrum, n):
    """Trim the leaderboard down to n elements with ties

    :param leaderboard: the current leaderboard
    :param spectrum: the current experimental spectrum
    :param n: max number of elements on the leaderboard
    :return: trimmed leaderboard
    """

    if len(leaderboard) < n:
        return leaderboard
    # build dict of scores of leaderboard
    score_set = dict()

    for seq in leaderboard:

        score = linear_score_mass_seq(seq, spectrum)
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
def leaderboard_cyclopeptide_sequencing(expands, spectrum, n, scoring="linear"):
    """Perform the cyclo peptide sequencing and find a best scoring amino acid sequence

    :param expands: the set of amino acid masses used to expand the leaderboard by
    :param spectrum: the experimental spectrum of the peptide
    :param n: the maximum number of elements allowed on the leaderboard
    :param scoring: the scoring parameter, can be "cyclic" or "linear"; defaults to "linear"
    :return: the higest scoring amino acid mass sequence given the experimental spectrum
    """

    leaderboard = [[]]
    leader_peptide = [0]

    # base case
    if n == 0 or spectrum == [0]:
        return [0]

    while leaderboard:
        leaderboard = expand(leaderboard, expands)
        leaderboard_iter = deepcopy(leaderboard)
        for peptide_mass_seq in leaderboard_iter:
            if peptide_mass(peptide_mass_seq) == parent_mass(spectrum):
                if score_peptide(peptide_mass_seq, spectrum, scoring_type=scoring) > \
                        score_peptide(leader_peptide, spectrum, scoring_type=scoring):
                    leader_peptide = peptide_mass_seq
            elif peptide_mass(peptide_mass_seq) > parent_mass(spectrum):
                leaderboard.remove(peptide_mass_seq)
        leaderboard = trim(leaderboard, spectrum, n)
    return leader_peptide


if __name__ == '__main__':

    filename = "sample_data.txt"

    m, n, spectrum = parse_data(filename)
    convolutions, conv_dict = spectral_convolution(spectrum)
    convolutions_set = convert_convolution_set_to_tuple(conv_dict)
    print(f"Convolved: {convolutions}")
    expands = expand_set_with_ties(convolutions_set, m)
    print(f"Expands: {expands}")

    # options for scoring include "cyclic" or "linear:
    leader = leaderboard_cyclopeptide_sequencing(expands, spectrum, n, scoring="cyclic")
    formatted = format_peptide_mass_seq(leader)
    print(f"leader: {formatted}")