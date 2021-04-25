"""
main.py file that runs the unittests when the file is called or run using the python CLI.

Part 1 of HW 3 for CISC 471, Computational Biology.

By: Hershil Devnani (20001045)

There are two ways to run the program, outlined below. Both methods run the unittests.
Sample Usage:
  $ python -m unittest herlpers_unit.py
  $ python -m main main.py

You can also modify the code in this file, main.py, in the main function to try different data sets and
parameters for each function by uncommenting and modifying the commented out function calls.
"""

import unittest
from helpers import *


class HelpersUnitTests(unittest.TestCase):
    """Testing class for the Helpers class
    """

    def test_profile_build(self):
        d = ["TAAC", "GTCT", "ACTA", "AGGT"]
        profile = get_profile(d, 4)

        correct_profile = [[0.5, 0.25, 0.25, 0.25],
                           [0.0, 0.25, 0.25, 0.25],
                           [0.25, 0.25, 0.25, 0.0],
                           [0.25, 0.25, 0.25, 0.5]]

        self.assertEqual(profile, correct_profile)

    def test_profile_pseudocounts_build(self):
        d = ["TAAC", "GTCT", "ACTA", "AGGT"]
        profile = get_profile_pseudocounts(d, 4)

        correct_profile = [[0.375, 0.25, 0.25, 0.25],
                           [0.125, 0.25, 0.25, 0.25],
                           [0.25, 0.25, 0.25, 0.125],
                           [0.25, 0.25, 0.25, 0.375]]

        self.assertEqual(profile, correct_profile)

    def test_motif_matrix_score(self):
        sample = ["TCGGGGGTTTTT", "CCGGTGACTTAC", "ACGGGGATTTTC", "TTGGGGACTTTT", "AAGGGGACTTCC", "TTGGGGACTTCC",
                  "TCGGGGATTCAT", "TCGGGGATTCCT", "TAGGGGAACTAC", "TCGGGTATAACC"]
        score = motif_matrix_score(sample, 12)

        self.assertEqual(score, 30)

    def test_profile_most_probable_a(self):
        d = "ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT"
        k = 5
        profile = [[0.2, 0.2, 0.3, 0.2, 0.3],
                   [0.4, 0.3, 0.1, 0.5, 0.1],
                   [0.3, 0.3, 0.5, 0.2, 0.4],
                   [0.1, 0.2, 0.1, 0.1, 0.2]]

        profile_most_test = profile_most_probable(d, k, profile)
        self.assertEqual(profile_most_test, "CCGAG")

    def test_profile_most_probable_b(self):
        d = "ACTATCTAAGGATAGTCTGACGACGGATTGATGGATCTATTCTCGGTCTCGGCTAGGCTGATAAGGATATAGCATAAGGTGACGCTTTAACTATCTGGCTGAC" \
            "TGACGCGGTTCAGTGACGTTAAGCCTTATGTCGACGTGGCGCAATAGGGAAGTTTCAGCGATGAATTTGCCCCCGACGCGGAGTTTTATAGACCAGT"
        k = 7

        profile = [[0.25, 0.107, 0.25, 0.214, 0.143, 0.25, 0.214],
                   [0.107, 0.429, 0.25, 0.143, 0.321, 0.25, 0.286],
                   [0.321, 0.107, 0.25, 0.357, 0.357, 0.179, 0.286],
                   [0.321, 0.357, 0.25, 0.286, 0.179, 0.321, 0.214]]

        profile_most_test = profile_most_probable(d, k, profile)
        self.assertEqual(profile_most_test, "TCGGCTA")


if __name__ == '__main__':
    unittest.main()
