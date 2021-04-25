# CISC 471 HW6

This project contains the programming implementation for homework 6. The code was written using Python 3.8+.

## File structure

The project contains the following files:

```
.
├── main.py                         # Runs the unittests the algorithms implemented for section one
├── manhattan_tourist.py            # Contains the implementations of the required algorithms for section one
├── manhattan_tourist_tests.py      # Contains the unit tests for manhattan_tourist.py
├── parsers.py                      # Contains functions to parse file input data
├── roslaind_ba5b_x.txt             # Contain sample data sets from Rosalind for the unittests
└── README.md                       # This file, contians information about the program
```

## Running the programs

The following commands can be used to run the programming question - The Manhattan Tourist Problem:

- `python -m main main.py`
- `python -m unittest manhattan_tourist_tests.py`
- `python manhattan_tourist_tests.py`

## 1 - The Manhattan Tourist Problem + Unit Tests

The solutions for the programming question are contained in `manhattan_tourist.py` and can be run from `main.py` or 
`manhattan_tourist_tests.py` depending on whether you would like to run the individual functions or the unit tests. See below for more detail.

### Running unittests and individual functions

Unittests:
- Running `main.py` or `manhattan_tourist_tests.py` using the commands above will run the unittests written for the
required functions.

Individual Functions:
- There are some commented out lines in the main method of `manhattan_tourist.py` and these can be uncommented or
modified to run the functions individually.
