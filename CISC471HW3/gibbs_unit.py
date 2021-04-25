"""
main.py file that runs the unittests when the file is called or run using the python CLI.

Part 1 of HW 3 for CISC 471, Computational Biology.

By: Hershil Devnani (20001045)

There are two ways to run the program, outlined below. Both methods run the unittests.
Sample Usage:
  $ python -m unittest gibbs_unit.py
  $ python -m main main.py

You can also modify the code in this file, main.py, in the main function to try different data sets and
parameters for each function by uncommenting and modifying the commented out function calls.
"""

import unittest
from gibbs import *
from helpers import *


class GibbsUnitTests(unittest.TestCase):
    """Testing class for the Gibbs Sampler Algorithm
    """

    def test_gibbs_positive(self):
        k, t, N, dna = gibbs_parse_data("gibbs_positive_sample.txt")
        sample_motifs = gibbs_sampler(dna, k, t, N)
        mismatch_percentage = 0

        for i in range(200):
            best_motifs_new = gibbs_sampler(dna, k, t, N)
            if motif_matrix_score(best_motifs_new, k) < motif_matrix_score(sample_motifs, k):
                sample_motifs = best_motifs_new

        solution = ['ACGTCCACCGGCGTC', 'AAGCGCACCGGGGTG', 'ACCCTTACCGGGGTG', 'AAGTTCCTCGGGGTG', 'AAGTTTTATGGGGTG',
                    'AAGTTTACCGGGTGC', 'AAGTTTCGAGGGGTG', 'CTGTTTACCGGGGTA', 'AAGTTGCTCGGGGTG', 'AAACATACCGGGGTG',
                    'AAGTTTAGGAGGGTG', 'AAGGAAACCGGGGTG', 'AAGTTTACACAGGTG', 'TAGTTTACCGGGGAT', 'CCTTTTACCGGGGTG',
                    'AAGTGAGCCGGGGTG', 'AAGTCGTCCGGGGTG', 'AAGTTTACCGGACAG', 'AAGTTTACCAATGTG', 'AAGTTTACCGTCATG']

        solution_score = 70

        for i in sample_motifs:
            if i not in solution:
                mismatch_percentage += 1

        mismatch_percentage = mismatch_percentage / len(sample_motifs)

        self.assertLessEqual(motif_matrix_score(sample_motifs, k), solution_score)
        self.assertLessEqual(mismatch_percentage, 0.1)  # percentage of mismatches lass than 10%

    def test_gibbs_negative(self):
        k, t, N, dna = gibbs_parse_data("gibbs_negative_sample.txt")
        sample_motifs = gibbs_sampler(dna, k, t, N)

        for i in range(200):
            best_motifs_new = gibbs_sampler(dna, k, t, N)
            if motif_matrix_score(best_motifs_new, k) < motif_matrix_score(sample_motifs, k):
                sample_motifs = best_motifs_new

        solution = []

        self.assertEqual(sample_motifs, solution)

    def test_gibbs_rosalind_sample_a(self):
        k, t, N, dna = gibbs_parse_data("rosalind_ba2g_gibbs_1.txt")
        sample_motifs = gibbs_sampler(dna, k, t, N)
        mismatch_percentage = 0

        for i in range(200):
            best_motifs_new = gibbs_sampler(dna, k, t, N)
            if motif_matrix_score(best_motifs_new, k) < motif_matrix_score(sample_motifs, k):
                sample_motifs = best_motifs_new

        solution = ['GGAACGACGTGATAC', 'CCCGATAAGTGAACA', 'CCCTCGCGATGAACA', 'CCCTCGCGATGAACA', 'TCCTCGAAGTGAAGT',
                    'CCCTCGACAGGAACA', 'CCCTGTCAGTGAACA', 'CGTCCGAAGTGAACA', 'TAATCGAAGTGAACA', 'CCTCGGAAGTGAACA',
                    'CCCTCGAAGCATACA', 'CCCTCGAAGTCGCCA', 'CCCTTTCAGTGAACA', 'CCCTCGAATCTAACA', 'CCCCATAAGTGAACA',
                    'TTCTCGAAGTGAACT', 'CCCTCCCGGTGAACA', 'CCCTCTCGGTGAACA', 'CCCTCGAAGTGAGTG', 'CCCTCGAAGTGCGAA']

        solution_score = 70

        for i in sample_motifs:
            if i not in solution:
                mismatch_percentage += 1

        mismatch_percentage = mismatch_percentage / len(sample_motifs)

        self.assertLessEqual(motif_matrix_score(sample_motifs, k), solution_score)
        self.assertLessEqual(mismatch_percentage, 0.1)  # percentage of mismatches lass than 10%

    def test_gibbs_rosalind_sample_b(self):
        k, t, N, dna = gibbs_parse_data("rosalind_ba2g_gibbs_2.txt")
        sample_motifs = gibbs_sampler(dna, k, t, N)
        mismatch_percentage = 0

        for i in range(200):
            best_motifs_new = gibbs_sampler(dna, k, t, N)
            if motif_matrix_score(best_motifs_new, k) < motif_matrix_score(sample_motifs, k):
                sample_motifs = best_motifs_new

        solution = ['CCCACCACCAGTTAA', 'CCCACCACCAGTTAA', 'CGCAAATGTTGTTAA', 'CGCACCATCTGTTAA', 'CGCACCGGTTCAGAA',
                    'CGCACTCCTTGTTAA', 'CGCTAAGGTTGTTAA', 'CGCACGCCTTGTTAA', 'GGCACCGGTTGTTTT', 'CGCACCGCCAGTTAA',
                    'TTGACCGGTTGTTAA', 'CGCACCGGTTGCGTA', 'CGCATGAGTTGTTAA', 'CGCACCGGTCCGTAA', 'CGTTGCGGTTGTTAA',
                    'ATCACCGGTTGTTAT', 'CGCACCGGCATTTAA', 'CGCACCGGTTGTATT', 'CATTCCGGTTGTTAA', 'CGCACCTCCTGTTAA']

        solution_score = 68

        for i in sample_motifs:
            if i not in solution:
                mismatch_percentage += 1

        mismatch_percentage = mismatch_percentage / len(sample_motifs)

        self.assertLessEqual(motif_matrix_score(sample_motifs, k), solution_score)  # score matches
        self.assertLessEqual(mismatch_percentage, 0.1)  # percentage of mismatches lass than 10%

    def test_gibbs_rosalind_sample_c(self):
        k, t, N, dna = gibbs_parse_data("rosalind_ba2g_gibbs_3.txt")
        sample_motifs = gibbs_sampler(dna, k, t, N)
        mismatch_percentage = 0

        for i in range(200):
            best_motifs_new = gibbs_sampler(dna, k, t, N)
            if motif_matrix_score(best_motifs_new, k) < motif_matrix_score(sample_motifs, k):
                sample_motifs = best_motifs_new

        solution = ['GTATTTAGCCGAATT', 'TTACTTGGGCGAATT', 'TTAGAAGGAATAATT', 'ATAGAAGGGCGAAAC', 'CGTGAAGGGCGAATT',
                    'TCTAAAGGGCGAATT', 'TTAGAAACTCGAATT', 'TTAGAAGGGCTCGTT', 'TTAGTTTGGCGAATT', 'TTAGCCTGGCGAATT',
                    'TTAGATCTGCGAATT', 'TTATCGGGGCGAATT', 'CGACGAGGTGGATTT', 'TTAGAAGCAAGAATT', 'TTGTTAGGGCGAATT',
                    'TTAGAAAAACGAATT', 'TTAGAAGGGATCATT', 'AGAGAAGGGCGAATG', 'TTAGAAGGGCGACGA', 'TTAGAAGGGCGCCAT']

        solution_score = 68

        for i in sample_motifs:
            if i not in solution:
                mismatch_percentage += 1

        mismatch_percentage = mismatch_percentage / len(sample_motifs)

        self.assertLessEqual(motif_matrix_score(sample_motifs, k), solution_score)
        self.assertLessEqual(mismatch_percentage, 0.1)  # percentage of mismatches lass than 10%


if __name__ == '__main__':
    unittest.main()
