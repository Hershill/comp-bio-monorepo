"""
main.py file that runs the unittests when the file is called or run using the python CLI.

Part 1 of HW 3 for CISC 471, Computational Biology.

By: Hershil Devnani (20001045)

There are two ways to run the program, outlined below. Both methods run the unittests.
Sample Usage:
  $ python -m unittest randomized_unit.py
  $ python -m main main.py

You can also modify the code in this file, main.py, in the main function to try different data sets and
parameters for each function by uncommenting and modifying the commented out function calls.
"""

import unittest
from randomized import *
from helpers import *


class RandomizedUnitTests(unittest.TestCase):
    """Testing class for the Randomized Motif Search Algorithm
    """

    def test_randomized_positive(self):
        k, t, dna = parse_data("randomized_positive_sample.txt")
        best_motifs = randomized_motif_search(dna, k, t)
        mismatch_percentage = 0

        for i in range(1000):
            best_motifs_new = randomized_motif_search(dna, k, t)
            if motif_matrix_score(best_motifs_new, k) < motif_matrix_score(best_motifs, k):
                best_motifs = best_motifs_new

        solution = ['CATGGGGAAAACTGA', 'CCTCTCGATCACCGA', 'CCTATAGATCACCGA', 'CCGATTGATCACCGA', 'CCTTGTGCAGACCGA',
                    'CCTTGCCTTCACCGA', 'CCTTGTTGCCACCGA', 'ACTTGTGATCACCTT', 'CCTTGTGATCAATTA', 'CCTTGTGATCTGTGA',
                    'CCTTGTGATCACTCC', 'AACTGTGATCACCGA', 'CCTTAGTATCACCGA', 'CCTTGTGAAATCCGA', 'CCTTGTCGCCACCGA',
                    'TGTTGTGATCACCGC', 'CACCGTGATCACCGA', 'CCTTGGTTTCACCGA', 'CCTTTGCATCACCGA', 'CCTTGTGATTTACGA']

        for i in best_motifs:
            if i not in solution:
                mismatch_percentage += 1

        mismatch_percentage = mismatch_percentage / len(best_motifs)

        self.assertLessEqual(mismatch_percentage, 0.1)  # percentage of mismatches lass than 10%

    def test_randomized_negative(self):
        k, t, dna = parse_data("randomized_negative_sample.txt")
        best_motifs = randomized_motif_search(dna, k, t)

        for i in range(1000):
            best_motifs_new = randomized_motif_search(dna, k, t)
            if motif_matrix_score(best_motifs_new, k) < motif_matrix_score(best_motifs, k):
                best_motifs = best_motifs_new

        solution = []

        self.assertEqual(best_motifs, solution)

    def test_randomized_rosalind_sample_a(self):
        k, t, dna = parse_data("rosalind_ba2f_randomized_1.txt")
        best_motifs = randomized_motif_search(dna, k, t)
        mismatch_percentage = 0

        for i in range(1000):
            best_motifs_new = randomized_motif_search(dna, k, t)
            if motif_matrix_score(best_motifs_new, k) < motif_matrix_score(best_motifs, k):
                best_motifs = best_motifs_new

        solution = ['CATATACCGAGGACG', 'CCTAAAGCGTTGACC', 'CCTCTCCTTTTGACC', 'CCTCTCGCCCAGACC', 'CCTCTCCTTTTGACC',
                    'ACTCTCGCGTTGAGG', 'CCGTACGCGTTGACC', 'CCTCTCGTCCTGACC', 'CCTCGATCGTTGACC', 'CCTCGATCGTTGACC',
                    'GATCTCGCGTTGACT', 'AGGCTCGCGTTGACC', 'CCTCTATAGTTGACC', 'CCTCTCGCGTTGGGA', 'CCTCTGCAGTTGACC',
                    'CCTGCGGCGTTGACC', 'CCTCTCGCGGCCACC', 'CACTTCGCGTTGACC', 'CCTCTCGCGTTAGAC', 'CCTCTCGCGTGCGCC']

        for i in best_motifs:
            if i not in solution:
                mismatch_percentage += 1

        mismatch_percentage = mismatch_percentage / len(best_motifs)

        self.assertLessEqual(mismatch_percentage, 0.1)  # percentage of mismatches lass than 10%

    def test_randomized_rosalind_sample_b(self):
        k, t, dna = parse_data("rosalind_ba2f_randomized_2.txt")
        best_motifs = randomized_motif_search(dna, k, t)
        mismatch_percentage = 0

        for i in range(1000):
            best_motifs_new = randomized_motif_search(dna, k, t)
            if motif_matrix_score(best_motifs_new, k) < motif_matrix_score(best_motifs, k):
                best_motifs = best_motifs_new

        solution = ['TGGACTGCTTTCAAG', 'TGGCTGACCCACAAG', 'ATGACAACCCACAAC', 'TGGAGTTCCCACAAG', 'TGGACAGAGCACAAG',
                    'TGGACAGTGCACAAG', 'TGGACAACCCTATAG', 'GGGACAACCCACAGC', 'TGGACGCTCCACAAG', 'TGGACAACCCATCTG',
                    'TACTCAACCCACAAG', 'TGGACAAAAAACAAG', 'TGGACAACCCACCTA', 'TGGACCCGCCACAAG', 'TGGACAACTACCAAG',
                    'TGTTAAACCCACAAG', 'TGGCGTACCCACAAG', 'AATACAACCCACAAG', 'TGGATGGCCCACAAG', 'TGGACAACCATAAAG']

        for i in best_motifs:
            if i not in solution:
                mismatch_percentage += 1

        mismatch_percentage = mismatch_percentage / len(best_motifs)

        self.assertLessEqual(mismatch_percentage, 0.1)  # percentage of mismatches lass than 10%

    def test_randomized_rosalind_sample_c(self):
        k, t, dna = parse_data("rosalind_ba2f_randomized_3.txt")
        best_motifs = randomized_motif_search(dna, k, t)
        mismatch_percentage = 0

        for i in range(1000):
            best_motifs_new = randomized_motif_search(dna, k, t)
            if motif_matrix_score(best_motifs_new, k) < motif_matrix_score(best_motifs, k):
                best_motifs = best_motifs_new

        solution = ['TCGGTTCGGCCACGA', 'TGTTATCTGCAATGC', 'TGTGTTAGGCAATGC', 'TGTGCGATGCAATGC', 'TGTGTATCTCAATGC',
                    'TGTGTACTCGGATGC', 'TGTGTACTGCACAAC', 'TTGCTACTGCAATGC', 'TGTGTACTGCCGGGC', 'CATGTACTGCAATGA',
                    'TGTGTACTGAGGTGC', 'TGTGTACTGCAACCT', 'TGTGTACGCAAATGC', 'TGTGTAAACCAATGC', 'CGTGTACTGCAATCG',
                    'TGTGTTAAGCAATGC', 'TGTCAGCTGCAATGC', 'ATGGTACTGCAATGC', 'TGCTAACTGCAATGC', 'TGTGAGATGCAATGC']

        for i in best_motifs:
            if i not in solution:
                mismatch_percentage += 1

        mismatch_percentage = mismatch_percentage / len(best_motifs)

        self.assertLessEqual(mismatch_percentage, 0.1)  # percentage of mismatches lass than 10%


if __name__ == '__main__':
    unittest.main()
