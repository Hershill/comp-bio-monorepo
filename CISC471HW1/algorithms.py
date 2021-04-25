"""
algorithms.py file containing the implementation of the algorithms for the 'Part 1 - Programming' section

Part 1 of HW 1 for CISC 471, Computational Biology.

By: Hershil Devnani (20001045)

There are two ways to run the program, outlined below. Both methods run the unittests.
Sample Usage:
  $ python -m unittest unittests.py
  $ python -m main main.py

You can also modify the code in this file, main.py, in the main function to try different data sets and
parameters for each function by uncommenting and modifying the commented out function calls.
"""


# 1.1 Frequent k-mer
def frequent_k_mer(text, k):
    """ Finds the set of k_mers that occurred most frequently in text

    :param text: sequence to search through
    :param k: length of k_mer to look for
    :return: list of k_mers most frequently appearing in text
    """

    # init vars
    k_mer_len = len(text)
    k_mers = dict()  # dictionary with all slices of k_mers and the count of how often they show up
    i = 0  # iterator

    # find the most frequent k-mer in a string
    while i <= k_mer_len - k:
        next_k_mer = text[i:i + k]
        # check if word in dict already
        if next_k_mer in k_mers.keys():
            k_mers[next_k_mer] += 1
        else:
            # otherwise add k letter word to dict
            k_mers[next_k_mer] = 1
        # increment i for loop
        i += 1
    return max_dict_vals(k_mers)


# 1.2 Frequent k-mer with mismatches
def frequent_k_mer_mismatch(text, k, d):
    """ Finds the set of k_mers that occurred most frequently in text with at most d mismatches

    :param text: sequence to search through
    :param k: length of k_mer to look for
    :param d: max allowance for  mismatches
    :return: list of k_mers most frequently appearing in text with at most d mismatches
    """

    # init vars
    k_mer_len = len(text)
    k_mers = dict()  # dictionary with all slices of k_mers and the count of how often they show up
    i = 0  # iterator

    # find the most frequent k-mer in a string
    while i <= k_mer_len - k:
        num_matches = 0
        next_k_mer = text[i:i + k]
        # check hamming distance to all existing keys in dict
        for k_mer in k_mers.keys():
            hamming_dist = hamming_distance(k_mer, next_k_mer)
            # if hamming dist to next_k_mer <= 1 add one to value
            if hamming_dist <= d:
                num_matches += 1
                k_mers[k_mer] += 1

        # add new k_mer and matched values num to dict
        if next_k_mer not in k_mers.keys():
            k_mers[next_k_mer] = num_matches

        # increment i for loop
        i += 1
    return max_dict_vals(k_mers)


def hamming_distance(k_mer_a, k_mer_b):
    """ Compare k_mer_b with k_mer_a, both should be the same length of k_mer

    :param k_mer_a: comparator a
    :param k_mer_b: comparator b
    :return: return hamming distance between k_mer_a and k_mer_b
    """

    assert len(k_mer_a) == len(k_mer_b)  # make sure k_mer lens are equals
    hamming_dist = 0
    for base_indx in range(len(k_mer_a)):
        if k_mer_a[base_indx] != k_mer_b[base_indx]:
            hamming_dist += 1
    return hamming_dist


def max_dict_vals(k_mers_dict):
    """ Find keys of the max value found in dictionary

    :param k_mers_dict: dictionary to find keys with max value from
    :return: keys that had the max value in k_mers_dict
    """

    found_k_mers = list()
    k_mer_values = k_mers_dict.values()  # returns object of all values in the dict
    max_k_mer_value = max(k_mer_values, default=None)  # value of k mer found the most amount of times
    for k_mer, k_mer_len in k_mers_dict.items():
        if k_mer_len == max_k_mer_value:
            # append the keys of the k_mers found the most times to a list
            found_k_mers.append(k_mer)
    return found_k_mers
