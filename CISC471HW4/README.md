# CISC 471 HW4

This project contains the programming implementation for homework 4. The code was written using Python 3.8+.

## File structure

The project contains the following files:

```
.
├── main.py                # Runs the unittests the algorithms implemented for section one
├── rna_translate.py       # Contains the implementations of the required algorithms for section one
├── rna_translate_unit.py  # Contains the unit tests for rna_translate.py
├── sample_data_x.py       # Contains the sample data for the unittests
├── solution_x.py          # Contains the solutions to the sample data for the unittests
└── README.md              # This file, contians information about the program
```

## Running the programs

The following commands can be used to run the programming question - The Peptide Encoding Problem:

- `python main.py`
- `python -m unittest rna_translate_unit.py`
- `python rna_translate_unit.py`

## 1 - Peptide Encoding Problem + Unit Tests

The solutions for the programming question are contained in `rna_translate.py` and can be run from `main.py` or 
`rna_translate_unit.py` depending on whether you would like to run the individual functions or the unit tests. See below for more detail.

### Running unittests and individual functions

Unittests:
- Running `main.py` or `rna_translate_unit.py` using the commands above will run the unittests written for the
required functions.

Individual Functions:
- There are some commented out lines in the main method of `rna_translate.py` and these can be uncommented or
modified to run the functions individually.
