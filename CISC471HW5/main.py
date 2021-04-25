"""
main.py file that runs the unittests when the file is called or run using the python CLI.

Part 1 of HW 5 for CISC 471, Computational Biology.

By:
    - Hershil Devnani (20001045)
    - Rayan Shaikli (20059806)

There are two ways to run the program, outlined below. Both methods run the unittests.
Sample Usage:
  $ python -m unittest unittests.py
  $ python -m main main.py
"""

from convolution_cyclopeptide_sequencing_tests import TestProgrammingPartOne
import unittest

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestProgrammingPartOne)
    unittest.TextTestRunner(verbosity=2).run(suite)
