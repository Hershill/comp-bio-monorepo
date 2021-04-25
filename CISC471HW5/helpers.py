"""
helpers.py file containing the implementation of many helper functions used across all three algorithms.

Part 1 of HW 5 for CISC 471, Computational Biology.

By:
    - Hershil Devnani (20001045)
    - Rayan Shaikli (20059806)

This class has many functions that are used in the algorithm to return values that can be easily accessible via a simple
function call
"""

import convolution_cyclopeptide_sequencing
from copy import deepcopy


# Function for getting amino acid mass
def amino_acid_mass(amino_acid):
    """Get the mss of an amino acid

    :param amino_acid: amino acid as a string
    :return: mass of the amino acid
    """

    return convolution_cyclopeptide_sequencing.AMINO_ACID_MASS[amino_acid]


def peptide_mass(peptide):
    """Get the mass of a peptide

    :param peptide: peptide as a list of amino acid integer masses
    :return: mass of the peptide
    """

    return sum(peptide)


# Function for parent mass
def parent_mass(experimental_spectrum):
    """Get the largest measured mass in the experimental spectrum

    :param experimental_spectrum: the experimental spectrum
    :return: the parent mass of the experimental spectrum
    """

    return max(experimental_spectrum)


def format_peptide_mass_seq(mass_seq):
    """Format a peptide into amino acid masses separated by a dah

    :param mass_seq: list of the seqeunce of masses of amino acids in a list
    :return: required formatting as a string of amino acid mass sequence
    """
    seq_str = ""

    for i in range(len(mass_seq) - 1):
        seq_str += str(mass_seq[i]) + "-"
    seq_str += str(mass_seq[-1])

    return seq_str


def formatted_to_mass_seq(formatted_str):
    """Convert peptide string to list containing sequence of masses of amino acids in the peptide

    :param formatted_str: peptide as a string
    :return: list of integers representing masses of amino acids in peptide
    """

    mass_seq = formatted_str.split("-")
    mass_seq = [int(i) for i in mass_seq]

    return mass_seq


def expand_set(convolutions, m):
    """Set of amino acids used to expand the leaderboard

    :param convolutions: convolution of spectrum
    :param m: max number of convolution set
    :return: set of amino acid masses to use as expand set for the leaderboard
    """

    remove_dups = list()
    subset = list()

    for i in convolutions:
        if i not in remove_dups:
            remove_dups.append(i)

    for i in remove_dups:
        if i > 57 or i < 200:
            subset.append(i)
        if len(subset) == m:
            break

    return subset


def expand_set_with_ties(conv_tuple_set, m):
    """Set of amino acids used to expand the leaderboard

    :param conv_tuple_set: convolution of spectrum
    :param m: max number of convolution set
    :return: set of amino acid masses to use as expand set for the leaderboard
    """

    temp = deepcopy(conv_tuple_set)
    subset = set()
    pair = tuple()

    for i in conv_tuple_set:
        pair = i
        mass = pair[0]
        if 57 <= mass <= 200:
            subset.add(mass)
            temp.remove(pair)
        if len(subset) == m:
            break

    for i in temp:
        mass = i[0]
        if 57 <= mass <= 200:
            if i[1] == pair[1]:
                subset.add(mass)

    return subset


def mass_to_peptide(mass_list):
    """Convert a peptide of mass sequences to a string

    :param mass_list: the peptide as a list of amino acid masses
    :return: peptide as a string of amino acids
    """

    peptide = ""

    for i in mass_list:
        peptide += convolution_cyclopeptide_sequencing.MASS_TO_AMINO_ACID[i]

    return peptide


def peptide_to_mass_seq(peptide):
    """Convert a peptide to a list of it's amino acid masses

    :param peptide: peptide as a string
    :return: peptide as a list of it's amino acid masses
    """

    mass_list = []

    for i in peptide:
        mass_list.append(convolution_cyclopeptide_sequencing.AMINO_ACID_MASS[i])

    return mass_list


def convert_convolution_set_to_tuple(convolutions_dict):
    """Convert set of convolution masses stored as a dictionary containing the number of occurrences of the mass into a
    sorted list of tuples with each element of the list representing a set (mass, number of occurrences in convolution)

    :param convolutions_dict: dictionary containing the convolution with mass and it's number of occurrences
    :return:
    """

    conv_tuple_set = [(k, v) for k, v in convolutions_dict.items()]
    conv_tuple_set.sort(key=lambda x: x[1], reverse=True)

    return conv_tuple_set

