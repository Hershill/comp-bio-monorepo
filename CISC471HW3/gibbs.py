"""
gibbs.py file containing the implementation of the Gibbs Sampler algorithms for the 'Part 1 - Programming' section

Part 1 of HW 3 for CISC 471, Computational Biology.

By: Hershil Devnani (20001045)

There are two ways to run the program, outlined below. Both methods run the unittests.
Sample Usage:
  $ python -m main gibbs.py

You can also modify the code in this file, main.py, in the main function to try different data sets and
parameters for each function by uncommenting and modifying the commented out function calls.
"""

import random
from helpers import *


def gibbs_sampler(dna, k, t, N):
    """The function that implements the Gibbs Sampler algorithm

    :param dna: set of DNA strings
    :param k: length of k-mer
    :param t: number of DNA strings in dna
    :param N: number of times to randomly select and replace motifs i.e. sample motifs
    :return: set of motifs found by the algorithm
    """
    if k ==0 or t == 0 or N == 0:
        return []

    random_starts = list()
    random_k_mers = list()

    # randomly select k-mers Motifs = (Motif_1, …, Motif_t) in each string from Dna
    for i in range(t):
        random_starts.append(random.randint(0, len(dna[0]) - k))
        # select random start index for k_mers

    for i in range(len(dna)):
        random_k_mers.append(dna[i][random_starts[i]:random_starts[i] + k])

    best_motifs = random_k_mers.copy()

    for j in range(0, N):
        # print(f"j: {j}")
        random_index = random.randint(0, t - 1)
        # removed_k_mer = best_motifs[random_index]
        removed_motif = dna[random_index]  # randomly selected sequence to remove
        best_motifs.pop(random_index)  # select a random sequence and remove it from set of motifs
        profile = get_profile_pseudocounts(best_motifs, k)  # create a profile from the remaining motifs

        # calculate probability of all k-mers in deleted string based on profile
        # results in n-k+1 probabilities (len - k + 1)
        probabilities = probabilities_k_mers(removed_motif, profile, k)

        # based on probabilities, generated a weighted random motif
        random_motif = random.choices(list(probabilities.keys()), weights=tuple(probabilities.values()), k=1)

        best_motifs = best_motifs[:random_index] + random_motif + best_motifs[random_index:]

        if motif_matrix_score(random_k_mers, k) < motif_matrix_score(best_motifs, k):
            best_motifs = random_k_mers.copy()

    return best_motifs


def probabilities_k_mers(dna_str, profile, k):
    """Calculate the probabilities of k-mers in a given DNA string

    :param dna_str: DNA string
    :param profile: profile to score the motifs on
    :param k: length of k-mer
    :return: the hightest scoring motif in the DNA string
    """
    k_mers = list()

    probabilities = dict()

    index = {"A": 0, "C": 1, "G": 2, "T": 3}

    for i in range(len(dna_str) - k):
        k_mer = dna_str[i:i + k]
        k_mers.append(k_mer)
    # add the last k-mer to the list
    k_mers.append(dna_str[-k:])

    # for each k-mer calculate it's probability based on profile
    for k_mer in k_mers:
        score = 1
        for i in range(len(k_mer)):
            score = score * profile[index[k_mer[i]]][i]
        probabilities[k_mer] = score  # probabilities for each k-mer in deleted sequence

    return probabilities


def gibbs_parse_data(filename):
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
    N = int(k_and_t[2])

    data_set.pop(0)

    for i in range(len(data_set)):
        data_set[i] = data_set[i].strip('\n')

    return k, t, N, data_set


def main():
    filename = "gibbs_positive_sample.txt"

    k, t, N, dna = gibbs_parse_data(filename)
    sample_motifs = gibbs_sampler(dna, k, t, N)

    for i in range(200):
        best_motifs_new = gibbs_sampler(dna, k, t, N)
        if not best_motifs_new:
            break
        if motif_matrix_score(best_motifs_new, k) < motif_matrix_score(sample_motifs, k):
            sample_motifs = best_motifs_new

    for i in sample_motifs:
        print(i)


if __name__ == '__main__':
    main()
