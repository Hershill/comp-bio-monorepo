# CISC 471 HW4

This project contains the programming implementation for homework 4. The code was written using Python 3.8+.

## File structure

The project contains the following files:

```
.
├── main.py                                       # Runs the unittests the algorithms implemented for section one
├── convolution_cyclopeptide_sequencing.py        # Contains the implementations of the required algorithms for section one
├── convolution_cyclopeptide_sequencing_tests.py  # Contains the unit tests for convolution_cyclopeptide_sequencing_tests.py
├── helpers.py                                    # Contains helper functions
├── parsers.py                                    # Contains parsers to parse file input data
├── xxx_test_roslaind_ba4x.txt                    # Contain the sample data setsfor the unittests
├── solution_x.py                                 # Contains the solutions to the sample data for the unittests
└── README.md                                     # This file, contians information about the program
```

## Running the programs

The following commands can be used to run the programming question - The Peptide Encoding Problem:

- `python main.py`
- `python -m unittest convolution_cyclopeptide_sequencing_tests.py`
- `python convolution_cyclopeptide_sequencing_tests.py`

## 1 - Peptide Encoding Problem + Unit Tests

The solutions for the programming question are contained in `convolution_cyclopeptide_sequencing.py` and can be run from `main.py` or 
`convolution_cyclopeptide_sequencing_tests.py` depending on whether you would like to run the individual functions or the unit tests. See below for more detail.

*Note: Due to large number of sub-problems and requiring testing, the unit tests take ~ 3 minutes to execute*

### Running unittests and individual functions

Unittests:
- Running `main.py` or `convolution_cyclopeptide_sequencing.py` using the commands above will run the unittests written for the
required functions.

Individual Functions:
- There are some commented out lines in the main method of `convolution_cyclopeptide_sequencing.py` and these can be uncommented or
modified to run the functions individually.
