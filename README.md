# comp-bio-monorepo
Monorepo for all assignments from CISC 471 - Computational Biology

This project contains the programming implementation for all the homeworks in CISC 471 - Computational Biology.
The code was written using Python 3.8+.

## File structure

The project contains the following files:

```
.
├── manhattan_tourist_tests.py      # Directories that contain the implementation of
├── parsers.py                      # Contains functions to parse file input data
├── roslaind_ba5b_x.txt             # Contain sample data sets from Rosalind for the unittests
└── README.md                       # This file, contians information about the program
```

## Running the programs

The following commands can be used to run the programming homework from within it's directory:

- `python -m main main.py`
- `python {algorithm}.py`

## 1 - The Manhattan Tourist Problem + Unit Tests

The solutions for the programming homeworks are contained in the corresponding homeworks subdirectory and can be run from `main.py` or 
`{algorithm}.py` depending on whether you would like to run the individual functions or the unit tests. See below for more detail.

### Running unittests and individual functions

Unittests:
- Running `main.py` using the commands above will run the unittests written for the required functions.

Individual Functions:
- There are some commented out lines in the main method of `{algorithm}.py` and these can be uncommented or
modified to run the functions individually.
