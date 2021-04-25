"""
persers.py file containing the implementation of parsing the data files for the algorithm

Part 1 of HW 6 for CISC 471, Computational Biology.

By: Hershil Devnani (20001045)

This class has a couple functions that are used in the algorithm to parse input data,
accessible via a simple function call
"""

import re


def parse_data(filename):
    """Read in the text file and return the required data

    :param filename: text files will graph data
    :return: m, n and weight matrices down and right
    """

    with open(filename) as file:
        data_set = file.readlines()

    # parse m and n
    m = re.split(' ', data_set[0].strip('\n'))[1]
    n = re.split(' ', data_set[0].strip('\n'))[0]

    m = int(m)
    n = int(n)

    data_set.pop(0)

    for i in range(len(data_set)):
        data_set[i] = data_set[i].strip('\n')

    split_index = data_set.index("-")

    down_elems = data_set[:split_index]
    right_elems = data_set[split_index + 1:]

    down, right = get_weight_matrices(down_elems, right_elems)

    return m, n, down, right


def get_weight_matrices(down_elems, right_elems):
    """Get the weight matrices for traversing down and right

    :param down_elems: list of strings for the weights of the nodes going down a column
    :param right_elems: list of strings for the weights of the nodes going across a row
    :return: nested list of integers for weight matrices down and right
    """
    down_list = list()
    right_list = list()

    for i in down_elems:
        down_list.append(list(map(int, i.split())))

    for i in right_elems:
        right_list.append(list(map(int, i.split())))

    return down_list, right_list
