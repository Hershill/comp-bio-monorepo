"""
main.py file that runs the unittests when the file is called or run using the python CLI.

Part 1 of HW 4 for CISC 471, Computational Biology.

By: Hershil Devnani (20001045)

There are two ways to run the program, outlined below. Both methods run the unittests.
Sample Usage:
  $ python -m unittest unittests.py
  $ python -m main main.py
"""

import unittest
from rna_translate import *


class TestProgrammingPartOne(unittest.TestCase):
    """Testing class for the required unittests
    """

    def test_rna_translation_positive_a(self):
        rna_str = parse_data("sample_data_a.txt")
        solution = parse_data("solution_a.txt")
        amino_seq = map_to_amino_acids(rna_str)

        # make sure solution matches computed result
        self.assertEqual(amino_seq, solution)

    def test_rna_translation_positive_b(self):
        rna_str = parse_data("sample_data_b.txt")
        solution = parse_data("solution_b.txt")
        amino_seq = map_to_amino_acids(rna_str)

        # make sure solution matches computed result
        self.assertEqual(amino_seq, solution)

    def test_rna_translation_positive_c(self):
        rna_str = parse_data("sample_data_c.txt")
        solution = parse_data("solution_c.txt")
        amino_seq = map_to_amino_acids(rna_str)

        # make sure solution matches computed result
        self.assertEqual(amino_seq, solution)

    def test_rna_translation_negative(self):
        rna_str = parse_data("negative_sample_data_a.txt")
        solution = None
        amino_seq = map_to_amino_acids(rna_str)

        # make sure solution matches computed result
        self.assertEqual(amino_seq, solution)


if __name__ == '__main__':
    unittest.main()
