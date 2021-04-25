"""
main.py file that runs the unittests when the file is called or run using the python CLI.

Part 1 of HW 3 for CISC 471, Computational Biology.

By: Hershil Devnani (20001045)

There are two ways to run the program, outlined below. Both methods run the unittests.
Sample Usage:
  $ python -m unittest greedy_unit.py
  $ python -m main main.py

You can also modify the code in this file, main.py, in the main function to try different data sets and
parameters for each function by uncommenting and modifying the commented out function calls.
"""

import unittest
from greedy import *
from helpers import *


class GreedyUnitTests(unittest.TestCase):
    """Testing class for the Greedy Motif Search Algorithm
    """

    def test_greedy_positive(self):
        k, t, dna = parse_data("greedy_positive_sample.txt")
        solution = ['CAG', 'CAG', 'CAA', 'CAA', 'CAA']
        motifs = greedy_motif_search(dna, k, t)

        self.assertEqual(motifs, solution)

    def test_greedy_negative(self):
        k, t, dna = parse_data("greedy_negative_sample.txt")
        solution = []
        motifs = greedy_motif_search(dna, k, t)

        self.assertEqual(motifs, solution)

    def test_greedy_rosalind_sample_a(self):
        k, t, dna = parse_data("rosalind_ba2d_greedy_1.txt")
        solution = ['TGGTGTAAGCAC', 'GGGGATGTGGAG', 'ATCTATCGCGTC', 'TTCATAACCTCG', 'TGGTAAGTGCTG', 'GTCTAAGACGCG',
                    'GGCGTTCTGTAC', 'TGGGATGCACGG', 'TGGTATGTACTG', 'TGCTTTGACGAG', 'TGGTATGTACGG', 'TGGTATGGACTG',
                    'TGGAGTGAAGTG', 'TTGAGTGTCTGG', 'TGGTATGGACAG', 'ACCCTTACGAAG', 'TGGTAAGAACAG', 'TGGTATGTACTG',
                    'TGGGAAGTACAG', 'TGGGAAGCACTG', 'TGCAATGAATCG', 'TGGGATGTACAG', 'TGGCAAGAACGG', 'TGGTATGTACGG',
                    'TGGCAAGGACAG']
        motifs = greedy_motif_search(dna, k, t)

        self.assertEqual(motifs, solution)

    def test_greedy_rosalind_sample_b(self):
        k, t, dna = parse_data("rosalind_ba2d_greedy_2.txt")
        solution = ['CGGGGCCTTTGC', 'CGTCGCTGGGTG', 'GCCGGCGTACAG', 'CCGGCTTGAATA', 'GGCCGCCGGTAA', 'GTCTCGATAGAA',
                    'CTTGGCATAGGG', 'CGTCCCGTAGTC', 'CCTTCGCAACGG', 'CTTCGCATGGGG', 'CTTGGCAATGGC', 'CTTGGCAAGGGC',
                    'CCCCGGAGGCGG', 'CTTCGCAAAGGC', 'GTCTGTAGAGGG', 'GGTGGTCTAGGC', 'CTTTGCATAGGC', 'CTTTGCAGGGGG',
                    'CCTTGGCAATGG', 'CTTCGCATTGGC', 'CTTCGCAAAGGG', 'CTTTGCATAGGA', 'CTTGGCAAAGGA', 'CCGGCGGTACGC',
                    'CTTGGCAGGGGC']
        motifs = greedy_motif_search(dna, k, t)

        self.assertEqual(motifs, solution)

    def test_greedy_rosalind_sample_c(self):
        k, t, dna = parse_data("rosalind_ba2d_greedy_3.txt")
        solution = ['ATTGTGACTAAC', 'AGTGTTTTGACA', 'CACTTAGGGGAG', 'GGAGAGGAGACA', 'ACCACAAGCGCA', 'TACGTGCCACTC',
                    'AATATGCGGACA', 'AATGCGCGAGTA', 'AGTGTAGGTATC', 'AAAGTGCCCAAC', 'AAAGTGCCTAAA', 'GGCGTGCAAATA',
                    'GGTGTGTCTGCG', 'AATGTGCCTACC', 'AATGTGCGTCTA', 'AATGTGCCGATC', 'AACGTGCCTATC', 'GAAGTTCGGGAC',
                    'ATAGTGGCAGTA', 'AACGTGCCTACC', 'AAAGTGCCGATA', 'AACATGCATATA', 'AAAGTGCCAACG', 'AATGTGCCTAAC',
                    'AACGCAAGTATA']
        motifs = greedy_motif_search(dna, k, t)

        self.assertEqual(motifs, solution)


if __name__ == '__main__':
    unittest.main()
