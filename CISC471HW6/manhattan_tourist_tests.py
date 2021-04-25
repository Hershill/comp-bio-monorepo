"""
manhattan_tourist_tests.py file that runs the unittests when the file is called or run using the python CLI.

Part 1 of HW 6 for CISC 471, Computational Biology.

By: Hershil Devnani (20001045)

There are two ways to run the program, outlined below. Both methods run the unittests.

Sample Usage:
  $ python -m unittest manhattan_tourist_tests.py
  $ python -m main main.py
"""

import unittest
from manhattan_tourist import *


class TestProgrammingPartOne(unittest.TestCase):
    """Testing class for the required unittests and all Rosalind problems solved for the assignment
    """

    def test_manhattan_tourist_positive_a(self):
        m, n, down, right = parse_data("manhattan_tourist_sample_data_positive.txt")
        longest_path_len = manhattan_tourist(m, n, down, right)

        solution = 34

        self.assertEqual(solution, longest_path_len)

    def test_manhattan_tourist_positive_b(self):
        m, n, down, right = parse_data("manhattan_tourist_sample_data_positive_2.txt")
        longest_path_len = manhattan_tourist(m, n, down, right)

        solution = 84

        self.assertEqual(solution, longest_path_len)

    def test_manhattan_tourist_positive_c(self):
        m, n, down, right = parse_data("rosalind_ba5b_1.txt")
        longest_path_len = manhattan_tourist(m, n, down, right)

        solution = 96

        self.assertEqual(solution, longest_path_len)

    def test_manhattan_tourist_positive_d(self):
        m, n, down, right = parse_data("rosalind_ba5b_2.txt")
        longest_path_len = manhattan_tourist(m, n, down, right)

        solution = 54

        self.assertEqual(solution, longest_path_len)

    def test_manhattan_tourist_positive_e(self):
        m, n, down, right = parse_data("rosalind_ba5b_3.txt")
        longest_path_len = manhattan_tourist(m, n, down, right)

        solution = 102

        self.assertEqual(solution, longest_path_len)

    def test_manhattan_tourist_negative(self):
        m, n, down, right = parse_data("manhattan_tourist_sample_data_negative.txt")
        longest_path_len = manhattan_tourist(m, n, down, right)

        solution = 0

        self.assertEqual(solution, longest_path_len)