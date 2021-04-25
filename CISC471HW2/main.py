#!/usr/bin/python
# -*- coding: utf-8 -*-

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

from unittests import TestProgrammingPartOne
import unittest

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestProgrammingPartOne)
    unittest.TextTestRunner(verbosity=2).run(suite)
