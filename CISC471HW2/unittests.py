"""
main.py file that runs the unittests when the file is called or run using the python CLI.

Part 1 of HW 2 for CISC 471, Computational Biology.

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
from contigs import *


class TestProgrammingPartOne(unittest.TestCase):
    """Testing class for the required unittests
    """

    def test_euclerian_graph_positive_a(self):
        graph = parse_data("unittest_euclerian_cycle_positive_a.txt")
        cycle, unvisited = create_cycle(start_node(graph), graph)
        eulerian_cycle = find_euler_cycle(cycle, unvisited)

        for node in eulerian_cycle:
            vals = list(graph.values())
            vals_flat_list = [item for sublist in vals for item in sublist]
            keys = list(graph.keys())
            graph_count = vals_flat_list.count(node) + keys.count(node)
            cycle_count = eulerian_cycle.count(node)
            # self.assertTrue(cycle_count == graph_count)
            if eulerian_cycle[0] == node:
                self.assertTrue(cycle_count == graph_count)
            else:
                self.assertTrue((cycle_count + 1) == graph_count)

    def test_euclerian_graph_positive_b(self):
        graph = parse_data("unittest_euclerian_cycle_positive_b.txt")
        cycle, unvisited = create_cycle(start_node(graph), graph)
        eulerian_cycle = find_euler_cycle(cycle, unvisited)

        for node in eulerian_cycle:
            vals = list(graph.values())
            vals_flat_list = [item for sublist in vals for item in sublist]
            keys = list(graph.keys())
            graph_count = vals_flat_list.count(node) + keys.count(node)
            cycle_count = eulerian_cycle.count(node)
            # self.assertTrue(cycle_count == graph_count)
            if eulerian_cycle[0] == node:
                self.assertTrue(cycle_count == graph_count)
            else:
                self.assertTrue((cycle_count + 1) == graph_count)

    # def test_euclerian_graph_negative(self):
    #     graph = parse_data("unittest_euclerian_cycle_negative.txt")
    #     cycle, unvisited = create_cycle(start_node(graph), graph)
    #     eulerian_cycle = find_euler_cycle(cycle, unvisited)
    #
    #     for node in eulerian_cycle:
    #         vals = list(graph.values())
    #         vals_flat_list = [item for sublist in vals for item in sublist]
    #         keys = list(graph.keys())
    #         graph_count = vals_flat_list.count(node) + keys.count(node)
    #         cycle_count = eulerian_cycle.count(node)
    #         # self.assertTrue(cycle_count == graph_count)
    #         if eulerian_cycle[0] == node:
    #             self.assertTrue(cycle_count == graph_count)
    #         else:
    #             self.assertTrue((cycle_count + 1) == graph_count)

    def test_generate_contigs_positive_a(self):
        filename = "unittest_generate_contigs_positive_a.txt"
        graph = build_graph(filename)
        contigs = get_contigs(graph)

        CORRECT_CONTIGS = ['ATG', 'ATG', 'TGT', 'TGGA', 'CAT', 'GAT', 'AGA']

        self.assertCountEqual(contigs, CORRECT_CONTIGS)

    def test_generate_contigs_positive_b(self):
        filename = "unittest_generate_contigs_positive_b.txt"
        graph = build_graph(filename)
        contigs = get_contigs(graph)

        CORRECT_CONTIGS = ['TAAT', 'ATG', 'ATG', 'ATG', 'TGCCAT', 'TGG', 'TGTT', 'GGG', 'GGAT']

        self.assertCountEqual(contigs, CORRECT_CONTIGS)

    def test_generate_contigs_negative(self):
        filename = "unittest_generate_contigs_negative.txt"
        graph = build_graph(filename)
        contigs = get_contigs(graph)

        CORRECT_CONTIGS = []

        self.assertCountEqual(contigs, CORRECT_CONTIGS)


if __name__ == '__main__':
    unittest.main()
