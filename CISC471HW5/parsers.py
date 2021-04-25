"""
helpers.py file containing the implementation of many helper functions used across all three algorithms.

Part 1 of HW 5 for CISC 471, Computational Biology.

By:
    - Hershil Devnani (20001045)
    - Rayan Shaikli (20059806)

This class has many functions that are used in the algorithm to parse input data
accessible via a simple function call
"""

import re


def parse_spectrum_data(filename):
    """Read in the text file and return the required data

        :param filename: file containing RNA sequence
        :return: RNA sequence as a string
        """
    with open(filename) as file:
        data_set = file.read().replace("\n", "")
    return data_set


def parse_conv_data(filename):
    """Read in the text file and return the required data

    :param filename: text files convolution spectrum data
    :return: graph in the form of a dictionary
    """

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
    return


def parse_scoring_data(filename):
    """Read in the text file and return the required data

    :param filename: text files with scoring data
    :return: graph in the form of a dictionary
    """

    with open(filename) as file:
        data_set = file.readlines()

    peptide = re.split(' ', data_set[0].strip('\n'))
    peptide = peptide[0]

    spectrum = re.split(' ', data_set[1].strip('\n'))
    spectrum = [int(i) for i in spectrum]

    return peptide, spectrum


def parse_leaderboard_seqeuncing_data(filename):
    """Read in the text file and return the required data

    :param filename: text files with seqeuncing data
    :return: graph in the form of a dictionary
    """

    with open(filename) as file:
        data_set = file.readlines()

    # parse n
    n = re.split(' ', data_set[0].strip('\n'))
    n = int(n[0])

    spectrum = re.split(' ', data_set[1].strip('\n'))
    spectrum = [int(i) for i in spectrum]

    return n, spectrum


def parse_data(filename):
    """Read in the text file and return the required data

    :param filename: text files will graph data
    :return: graph in the form of a dictionary
    """

    with open(filename) as file:
        data_set = file.readlines()

    # parse k and t
    m = re.split(' ', data_set[0].strip('\n'))
    n = re.split(' ', data_set[1].strip('\n'))

    try:
        m = int(m[0])
    except:
        m = 0

    try:
        n = int(n[0])
    except:
        n = 0

    spectrum = re.split(' ', data_set[2].strip('\n'))
    try:
        spectrum = [int(i) for i in spectrum]
    except:
        spectrum = [0]

    return m, n, spectrum

