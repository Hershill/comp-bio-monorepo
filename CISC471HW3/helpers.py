"""
helpers.py file containing the implementation of many helper functions used across all three algorithms.

Part 1 of HW 3 for CISC 471, Computational Biology.

By: Hershil Devnani (20001045)

This class has many functions that are shared across all three algorithms.

"""
import re


def parse_data(filename):
    """Read in the text file and build the graph

    :param filename: text files will graph data
    :return: graph in the form of a dictionary
    """

    with open(filename) as file:
        data_set = file.readlines()

    # parse k and t
    k_and_t = re.split(' ', data_set[0].strip('\n'))
    k = int(k_and_t[0])
    t = int(k_and_t[1])

    data_set.pop(0)

    for i in range(len(data_set)):
        data_set[i] = data_set[i].strip('\n')

    return k, t, data_set


def motif_matrix_score(motifs, k):
    """

    :param motifs:
    :param k:
    :return:
    """
    # get occurrences of each base
    # occurrences = {"A": 0, "C": 0, "G": 0, "T": 0}
    # A 2
    # C 1
    # G 0
    # T 2

    index = {"A": 0, "C": 1, "G": 2, "T": 3}

    occur = [[0] * k for i in range(4)]

    # gives count matrix
    for motif in motifs:
        for i in range(len(motif)):
            occur[index[motif[i]]][i] += 1

    score = 0

    # sum each column then subtract the largest number
    for k in range(k):
        col = list()
        for i in range(len(occur)):
            col.append(occur[i][k])
        score += sum(col) - max(col)

    return score


def profile_most_probable(dna, k, profile):
    """

    :param dna:
    :param k:
    :param profile:
    :return:
    """
    k_mer_prob = -1
    selected_kmer = ""
    k_mers = list()

    index = {"A": 0, "C": 1, "G": 2, "T": 3}

    for i in range(len(dna) - k):
        k_mer = dna[i:i + k]
        k_mers.append(k_mer)
    # add the last k-mer to the list
    k_mers.append(dna[-k:])

    # for each k-mer calculate it's probability based on profile
    for k_mer in k_mers:
        score = 1
        for i in range(len(k_mer)):
            score = score * profile[index[k_mer[i]]][i]
        if score > k_mer_prob:
            selected_kmer = k_mer
            k_mer_prob = score

    # return k-mer with greatest probability
    return selected_kmer


def get_profile(dna_set, k):
    """

    :param dna_set:
    :param k:
    :return:
    """
    # profile 4 x k list, i.e.
    # A 0.2 0.2 0.3 0.2 0.3
    # C 0.4 0.3 0.1 0.5 0.1
    # G 0.3 0.3 0.5 0.2 0.4
    # T 0.1 0.2 0.1 0.1 0.2

    # dict/enum for easy lookup for adding to profile
    index = {"A": 0, "C": 1, "G": 2, "T": 3}

    # init the matrix
    profile = [[0.0] * k for i in range(4)]

    for dna in dna_set:
        for i in range(len(dna)):
            profile[index[dna[i]]][i] += 1/len(dna_set)

    return profile


def get_profile_pseudocounts(dna_set, k):
    """

    :param dna_set:
    :param k:
    :return:
    """
    # profile 4 x k list, i.e.
    # A 0.2 0.2 0.3 0.2 0.3
    # C 0.4 0.3 0.1 0.5 0.1
    # G 0.3 0.3 0.5 0.2 0.4
    # T 0.1 0.2 0.1 0.1 0.2

    # dict/enum for easy lookup for adding to profile
    index = {"A": 0, "C": 1, "G": 2, "T": 3}

    # init the matrix
    profile = [[0.0] * k for i in range(4)]

    for dna in dna_set:
        for i in range(len(dna)):
            profile[index[dna[i]]][i] += 1

    for counts in profile:
        for j in range(len(counts)):
            counts[j] += 1
            counts[j] = counts[j] / (len(dna_set) * 2)

    return profile
