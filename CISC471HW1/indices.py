"""
indices.py file that contains the implementation of the algorithm for '2.1 - Indices'

Part 2 of HW 1 for CISC 471, Computational Biology.

By: Hershil Devnani (20001045)

Sample Usage:
  $ python indices.py

Please see the accompanied file, hw1.pdf, for an explanation of the algorithm as well as answers to the remaining
question for question 2.1.
"""


# 2.1 Indices

# Algorithm that iterates over every index from (0, 0,..., 0) to (n1, n2,..., nd), where n1, n2,..., nd = Nj

def indices(a, d, n):
    """Iterates over indices from lists a to n

    :param a: initial list of indices
    :param d: length of list of indices
    :param n: value of list of n's, all values of N are equal
    :return: whether there are more output remaining in the sequence
    """
    while d > 0 and a[d - 1] == n[d - 1]:
        a[d - 1] = 0
        d = d - 1
    if d == 0:
        print(f"reset index {d - 1}")
        return False
    else:
        a[d - 1] = a[d - 1] + 1

        return True


def run(d, N):
    """Gather outputs and number of outputs based on given d and N values

    :param d: length of list of indices
    :param N: value of N
    :return: list of outputs
    """
    a_list = list()
    n_list = list()
    sequence = list()
    for i in range(d):
        a_list.append(0)
        n_list.append(N)
    sequence.append(a_list)
    while indices(a_list, d, n_list):
        sequence.append(a_list.copy())  # create a copy of a_list elements instead of the global object reference

    print(f"d: {d}, N: {N}, #Outputs: {len(sequence)}")
    print(f"Sequence {sequence}")


# 1,2,4,6,9,12,16,20,25,30,36
# 2,2,3,3,4,4,5,6,6,

# Sum of N = 0, Multiply of N = 0, d = 2, #output = 1 -- 2^0 = 1
# Sum of N = 1, Multiply of N = 0, d = 2, #output = 2 -- 2^1 = 2
# Sum of N = 3, Multiply of N = 2, d = 2, #output = 4 -- 2^3 = 4

# Sum of N = 10, Multiply of N = 25, d = 2, #output = 36 -- (11/2)^2
# Sum of N = 2, Multiply of N = 1, d = 2, #output = 2 -- 2^2 = 4


if __name__ == '__main__':
    # Runs indices for a given d and N, where n1,n2,...nd all equal N
    run(3, 1)
    # run(2, 2)
    # run(2, 3)
    # run(2, 4)
    # run(3, 1)
    # run(3, 2)
    # run(3, 3)
    # run(3, 4)

    # rule for number of sequences is (N+1)^d
