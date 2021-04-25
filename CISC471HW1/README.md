# CISC 471 HW1

This project contains the required implemented algorithms for HomeWork 1. The code was written using Python 3.9.

## File structure

The project contains the following files:

```
.
├── main.py                        # Can run the unitests or individual functions for the algorithms implements for section one
├── algorithms.py                  # Contains the implementations of the required algorithms for section one
├── unittests.py                   # Contains the unit tests for algorithms.py
├── indices.py                     # Contains the implementation of the algorithm for section 2.1 - Indices
├── dishonest_profs_pseudocode.py  # Contains the pseudocode for section 2.3 - Dishonest Professors
└── README.md
```

## Running the programs

The following commands can be used to run the first part of the homework, parts 1.1 to 1.3:

- `python main.py`
- `python -m unittest unittests.py`
- `python unittests.py`

The following command run the implementation of an algorithm for part 2.1 of the assignment:

- `python indices.py`

## 1.1 to 1.3 - Frequent k-mers + Unit Tests

The solutions for section 1 are contained in `algorithms.py` and can be run from `main.py` or `unittests.py`
depending on whether you would like to run the individual functions or the unit tests. See below for more detail.

### Running unittests and individual functions

Unittests:
- Running `main.py` or `unittests.py` using the commands above will run the unittests written for the
required functions.

Individual Functions:
- There are some commented out lines in the main method of `main.py` and these can be uncommented or
modified to run the functions individually.

## 2.1 - Indices

The parameters passed into the `run` function in `indices.py` can be modified to reflect the
desired input for d and N. The output will display the all output units as well as the total
number of outputs for any given combination of d and N passed to the `run` function.

## 2.3 - Dishonest Professors

This file `dishonest_profs_pseudocode.py` contains the pseudocode to section 2.3. There is a more detailed explanation
provided in `hw1.pdf` as well.
