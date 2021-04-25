"""
randomized.py file containing the implementation of the Gibbs Sampler algorithms for the 'Part 1 - Programming' section

Part 1 of HW 3 for CISC 471, Computational Biology.

By: Hershil Devnani (20001045)

There are two ways to run the program, outlined below. Both methods run the unittests.
Sample Usage:
  $ python -m main randomized.py

You can also modify the code in this file, main.py, in the main function to try different data sets and
parameters for each function by uncommenting and modifying the commented out function calls.
"""

import random
from helpers import *


def randomized_motif_search(dna, k, t):
    """The function that implements the Randomized Motif Search algorithm

    :param dna: set of DNA strings
    :param k: length of k-mer
    :param t: number of DNA strings in dna
    :return: set of motifs found by the algorithm
    """
    # randomly select k-mers Motifs = (Motif_1, …, Motif_t) in each string from dna

    if k == 0 or t == 0:
        return []

    random_starts = list()
    random_k_mers = list()

    for i in range(t):
        random_starts.append(random.randint(0, len(dna[0]) - k))
        # select random start index for k_mers

    for i in range(len(dna)):
        random_k_mers.append(dna[i][random_starts[i]:random_starts[i]+k])

    best_motifs = random_k_mers

    while True:
        profile = get_profile_pseudocounts(best_motifs, k)
        new_motifs = motifs(profile, dna, k)
        if motif_matrix_score(new_motifs, k) < motif_matrix_score(best_motifs, k):
            best_motifs = new_motifs
        else:
            return best_motifs


def motifs(profile, dna_set, k):
    """Returns a set of profile most probable motifs in a set of DNA

    :param profile: profile used for scoring the motifs
    :param dna_set: set of DNA strings
    :param k: length of k-mer
    :return: set of motifs based on profile most probable motif in each DNA string
    """
    new_motifs = []
    for dna in dna_set:
        motif = profile_most_probable(dna, k, profile)
        new_motifs.append(motif)
    return new_motifs


def main():
    filename = "randomized_positive_sample.txt"
    k, t, dna = parse_data(filename)
    best_motifs = randomized_motif_search(dna, k, t)
    # print(best_motifs)
    # for i in best_motifs:
    #     print(i)

    for i in range(1000):
        best_motifs_new = randomized_motif_search(dna, k, t)
        if not best_motifs_new:
            break
        if motif_matrix_score(best_motifs_new, k) < motif_matrix_score(best_motifs, k):
            best_motifs = best_motifs_new

    for i in best_motifs:
        print(i)


if __name__ == '__main__':
    main()
