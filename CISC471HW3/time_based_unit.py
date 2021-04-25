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
import time
from greedy import *
from randomized import *
from gibbs import *


class TimeBasedConstantKUnitTests(unittest.TestCase):
    """Testing class for the Real Time Analysis of the Algorithms against different value of T - the number of DNA strings
    """

    def setUp(self):
        self.startTime = time.time()

    def tearDown(self):
        t = time.time() - self.startTime
        print(f'K: {self.k}, T: {self.t}, {self.id()}: Avg. Time: {t/3:.3f}')

    def test_greedy_constant_k_1(self):
        for i in range(3):
            self.k, self.t, self.dna = parse_data("time_test_all_constant_k_1.txt")
            greedy_motif_search(self.dna, self.k, self.t)

    def test_greedy_constant_k_2(self):
        for i in range(3):
            self.k, self.t, self.dna = parse_data("time_test_all_constant_k_2.txt")
            greedy_motif_search(self.dna, self.k, self.t)

    def test_greedy_constant_k_3(self):
        for i in range(3):
            self.k, self.t, self.dna = parse_data("time_test_all_constant_k_3.txt")
            greedy_motif_search(self.dna, self.k, self.t)

    def test_random_constant_k_1(self):
        for i in range(3):
            self.k, self.t, self.dna = parse_data("time_test_all_constant_k_1.txt")
            best_motifs = randomized_motif_search(self.dna, self.k, self.t)

            for i in range(1000):
                best_motifs_new = randomized_motif_search(self.dna, self.k, self.t)
                if motif_matrix_score(best_motifs_new, self.k) < motif_matrix_score(best_motifs, self.k):
                    best_motifs = best_motifs_new

    def test_random_constant_k_2(self):
        for i in range(3):
            self.k, self.t, self.dna = parse_data("time_test_all_constant_k_2.txt")
            best_motifs = randomized_motif_search(self.dna, self.k, self.t)

            for i in range(1000):
                best_motifs_new = randomized_motif_search(self.dna, self.k, self.t)
                if motif_matrix_score(best_motifs_new, self.k) < motif_matrix_score(best_motifs, self.k):
                    best_motifs = best_motifs_new

    def test_random_constant_k_3(self):
        for i in range(3):
            self.k, self.t, self.dna = parse_data("time_test_all_constant_k_3.txt")
            best_motifs = randomized_motif_search(self.dna, self.k, self.t)

            for i in range(1000):
                best_motifs_new = randomized_motif_search(self.dna, self.k, self.t)
                if motif_matrix_score(best_motifs_new, self.k) < motif_matrix_score(best_motifs, self.k):
                    best_motifs = best_motifs_new

    def test_gibbs_constant_k_1(self):
        for i in range(3):
            self.k, self.t, self.N, self.dna = gibbs_parse_data("time_test_all_constant_k_1.txt")
            sample_motifs = gibbs_sampler(self.dna, self.k, self.t, self.N)

            for i in range(200):
                best_motifs_new = gibbs_sampler(self.dna, self.k, self.t, self.N)
                if motif_matrix_score(best_motifs_new, self.k) < motif_matrix_score(sample_motifs, self.k):
                    sample_motifs = best_motifs_new

    def test_gibbs_constant_k_2(self):
        for i in range(3):
            self.k, self.t, self.N, self.dna = gibbs_parse_data("time_test_all_constant_k_2.txt")
            sample_motifs = gibbs_sampler(self.dna, self.k, self.t, self.N)

            for i in range(200):
                best_motifs_new = gibbs_sampler(self.dna, self.k, self.t, self.N)
                if motif_matrix_score(best_motifs_new, self.k) < motif_matrix_score(sample_motifs, self.k):
                    sample_motifs = best_motifs_new

    def test_gibbs_constant_k_3(self):
        for i in range(3):
            self.k, self.t, self.N, self.dna = gibbs_parse_data("time_test_all_constant_k_3.txt")
            sample_motifs = gibbs_sampler(self.dna, self.k, self.t, self.N)

            for i in range(200):
                best_motifs_new = gibbs_sampler(self.dna, self.k, self.t, self.N)
                if motif_matrix_score(best_motifs_new, self.k) < motif_matrix_score(sample_motifs, self.k):
                    sample_motifs = best_motifs_new


class TimeBasedConstantTUnitTests(unittest.TestCase):
    """Testing class for the Real Time Analysis of the Algorithms against different value of K - the number of k-mer
    """

    def setUp(self):
        self.startTime = time.time()

    def tearDown(self):
        t = time.time() - self.startTime
        print(f'K: {self.k}, T: {self.t}, {self.id()}: Avg. Time: {t/3:.3f}')

    def test_greedy_constant_t_1(self):
        for i in range(3):
            self.k, self.t, self.dna = parse_data("time_test_all_constant_t_1.txt")
            greedy_motif_search(self.dna, self.k, self.t)

    def test_greedy_constant_t_2(self):
        for i in range(3):
            self.k, self.t, self.dna = parse_data("time_test_all_constant_t_2.txt")
            greedy_motif_search(self.dna, self.k, self.t)

    def test_greedy_constant_t_3(self):
        for i in range(3):
            self.k, self.t, self.dna = parse_data("time_test_all_constant_t_3.txt")
            greedy_motif_search(self.dna, self.k, self.t)

    def test_random_constant_t_1(self):
        for i in range(3):
            self.k, self.t, self.dna = parse_data("time_test_all_constant_t_1.txt")
            best_motifs = randomized_motif_search(self.dna, self.k, self.t)

            for i in range(1000):
                best_motifs_new = randomized_motif_search(self.dna, self.k, self.t)
                if motif_matrix_score(best_motifs_new, self.k) < motif_matrix_score(best_motifs, self.k):
                    best_motifs = best_motifs_new

    def test_random_constant_t_2(self):
        for i in range(3):
            self.k, self.t, self.dna = parse_data("time_test_all_constant_t_2.txt")
            best_motifs = randomized_motif_search(self.dna, self.k, self.t)

            for i in range(1000):
                best_motifs_new = randomized_motif_search(self.dna, self.k, self.t)
                if motif_matrix_score(best_motifs_new, self.k) < motif_matrix_score(best_motifs, self.k):
                    best_motifs = best_motifs_new

    def test_random_constant_t_3(self):
        for i in range(3):
            self.k, self.t, self.dna = parse_data("time_test_all_constant_t_3.txt")
            best_motifs = randomized_motif_search(self.dna, self.k, self.t)

            for i in range(1000):
                best_motifs_new = randomized_motif_search(self.dna, self.k, self.t)
                if motif_matrix_score(best_motifs_new, self.k) < motif_matrix_score(best_motifs, self.k):
                    best_motifs = best_motifs_new

    def test_gibbs_constant_t_1(self):
        for i in range(3):
            self.k, self.t, self.N, self.dna = gibbs_parse_data("time_test_all_constant_t_1.txt")
            sample_motifs = gibbs_sampler(self.dna, self.k, self.t, self.N)
            # mismatch_percentage = 0

            for i in range(200):
                best_motifs_new = gibbs_sampler(self.dna, self.k, self.t, self.N)
                if motif_matrix_score(best_motifs_new, self.k) < motif_matrix_score(sample_motifs, self.k):
                    sample_motifs = best_motifs_new

    def test_gibbs_constant_t_2(self):
        for i in range(3):
            self.k, self.t, self.N, self.dna = gibbs_parse_data("time_test_all_constant_t_2.txt")
            sample_motifs = gibbs_sampler(self.dna, self.k, self.t, self.N)
            # mismatch_percentage = 0

            for i in range(200):
                best_motifs_new = gibbs_sampler(self.dna, self.k, self.t, self.N)
                if motif_matrix_score(best_motifs_new, self.k) < motif_matrix_score(sample_motifs, self.k):
                    sample_motifs = best_motifs_new

    def test_gibbs_constant_t_3(self):
        for i in range(3):
            self.k, self.t, self.N, self.dna = gibbs_parse_data("time_test_all_constant_t_3.txt")
            sample_motifs = gibbs_sampler(self.dna, self.k, self.t, self.N)
            # mismatch_percentage = 0

            for i in range(200):
                best_motifs_new = gibbs_sampler(self.dna, self.k, self.t, self.N)
                if motif_matrix_score(best_motifs_new, self.k) < motif_matrix_score(sample_motifs, self.k):
                    sample_motifs = best_motifs_new


if __name__ == '__main__':
    unittest.main()
