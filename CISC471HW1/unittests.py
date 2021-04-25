"""
main.py file that runs the unittests when the file is called or run using the python CLI.

Part 1 of HW 1 for CISC 471, Computational Biology.

By: Hershil Devnani (20001045)

There are two ways to run the program, outlined below. Both methods run the unittests.
Sample Usage:
  $ python -m unittest unittests.py
  $ python -m main main.py

You can also modify the code in this file, main.py, in the main function to try different data sets and
parameters for each function by uncommenting and modifying the commented out function calls.
"""

import unittest
from algorithms import *


class TestProgrammingPartOne(unittest.TestCase):
    """Testing class for the required unittests
    """

    # def test_frequent_k_mer(self):
    #     self.assertEqual(frequent_k_mer("ACAACTATGCATCACTATCGGGAACTATCCT", 5), ['ACTAT'])

    def test_k_mer_positive(self):
        print("Output of test_k_mer_positive: " + str(frequent_k_mer("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4)))  # print the result
        self.assertEqual(frequent_k_mer("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4), ['GCAT', 'CATG'])  # assert the result

    def test_k_mer_mismatches_positive(self):
        print("Output of test_k_mer_mismatches_positive: " + str(frequent_k_mer_mismatch("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4, 2)))  # print the result
        self.assertEqual(frequent_k_mer_mismatch("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4, 2), ['ACGT', 'TGAT'])  # assert the result

    def test_k_mer_negative(self):
        print("Output of test_k_mer_negative: " + str(frequent_k_mer("ACGT", 5)))  # print the result
        self.assertEqual(frequent_k_mer("ACGT", 5), [])  # assert the result

    def test_k_mer_mismatches_negative(self):
        # print the result
        print("Output of test_k_mer_mismatches_negative: " + str(frequent_k_mer_mismatch("ATGTT", 6, 2)))
        self.assertEqual(frequent_k_mer_mismatch("ATGTT", 6, 2), [])  # assert the result


if __name__ == '__main__':
    unittest.main()
