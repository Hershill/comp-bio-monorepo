"""
main.py file that runs the unittests when the file is called or run using the python CLI.

Part 1 of HW 1 for CISC 471, Computational Biology.

By: Hershil Devnani (20001045)

Sample Usage to run the UnitTests:
  $ python -m main main.py

You can also modify the code in this file, main.py, in the main function to try different data sets and
parameters for each function by uncommenting and modifying the commented out function calls.
"""

from greedy_unit import GreedyUnitTests
from randomized_unit import RandomizedUnitTests
from gibbs_unit import GibbsUnitTests
from helpers_unit import HelpersUnitTests

import unittest

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(GreedyUnitTests)
    unittest.TextTestRunner(verbosity=2).run(suite)

    suite = unittest.TestLoader().loadTestsFromTestCase(RandomizedUnitTests)
    unittest.TextTestRunner(verbosity=2).run(suite)

    suite = unittest.TestLoader().loadTestsFromTestCase(GibbsUnitTests)
    unittest.TextTestRunner(verbosity=2).run(suite)

    suite = unittest.TestLoader().loadTestsFromTestCase(HelpersUnitTests)
    unittest.TextTestRunner(verbosity=2).run(suite)
