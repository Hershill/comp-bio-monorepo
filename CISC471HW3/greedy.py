"""
greedy.py file containing the implementation of the Gibbs Sampler algorithms for the 'Part 1 - Programming' section

Part 1 of HW 3 for CISC 471, Computational Biology.

By: Hershil Devnani (20001045)

There are two ways to run the program, outlined below. Both methods run the unittests.
Sample Usage:
  $ python -m main greedy.py

You can also modify the code in this file, main.py, in the main function to try different data sets and
parameters for each function by uncommenting and modifying the commented out function calls.
"""

from helpers import *


def greedy_motif_search(dna, k, t):
    """The function that implements the Greedy Motif Search algorithm

    :param dna: set of DNA strings
    :param k: length of k-mer
    :param t: number of DNA strings in dna
    :return: set of motifs found by the algorithm
    """

    if k == 0 or t == 0:
        return []

    # motif matrix formed by first k-mers in each string from Dna
    best_motifs = list()
    for i in range(len(dna)):
        best_motifs.append(dna[i][0:k])

    k_mers = list()  # all k_mers in first string from dna
    for i in range(len(dna[0]) - k):
        k_mers.append(dna[0][i:k+i])
    k_mers.append(dna[0][-k:])

    curr_motifs = list()

    # for each k-mer Motif in the first string from Dna
    for k_mer in k_mers:
        curr_motifs.append(k_mer)
        curr_motif = [k_mer]
        # form Profile from motifs Motif_1, …, Motif_i - 1
        for i in range(1, t):
            profile = get_profile(curr_motif, k)
            most_prob_next = profile_most_probable(dna[i], k, profile)
            curr_motif.append(most_prob_next)
        selected_motifs = curr_motif
        # get score of current motif set
        if motif_matrix_score(selected_motifs, k) < motif_matrix_score(best_motifs, k):
            best_motifs = selected_motifs

    return best_motifs


def main():
    filename = "greedy_positive_sample.txt"
    k, t, dna = parse_data(filename)
    motifs = greedy_motif_search(dna, k, t)

    for i in motifs:
        print(i)


if __name__ == '__main__':
    main()
