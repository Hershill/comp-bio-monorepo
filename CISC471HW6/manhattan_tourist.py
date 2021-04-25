"""
manhattan_tourist.py file containing the implementation of the Manhattan Tourist Problem algorithm for 'Part 1 - Programming'

Part 1 of HW 6 for CISC 471, Computational Biology.

By: Hershil Devnani (20001045)

There are two ways to run the program, outlined below. Both methods run the unittests.

Sample Usage:
  $ python -m unittest manhattan_tourist_tests.py
  $ python -m main main.py

You can also modify the code in this file, manhattan_tourist.py, in the main function to try different data sets and
parameters for each function by uncommenting and modifying the commented out function calls.
"""


from parsers import parse_data


def manhattan_tourist(m, n, down, right):
    """get the length of the longest path from (0, 0) to (n, m)

    :param m: number of nodes across of the graph
    :param n: number of nodes per column of the graph
    :param down: weights of the nodes as we move down a column
    :param right: weights of the nodes as we move across a row
    :return: the length of the longest path from source ot sink
    """
    s = list()

    for i in range(n + 1):
        s.append(list())
        for k in range(m + 1):
            s[i].append(0)

    for i in range(1, n + 1):
        s[i][0] = s[i-1][0] + down[i-1][0]

    for j in range(1, m + 1):
        s[0][j] = s[0][j-1] + right[0][j-1]

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s[i][j] = max(s[i-1][j] + down[i-1][j], s[i][j-1] + right[i][j-1])

    return s[n][m]


if __name__ == '__main__':
    # filename = "manhattan_tourist_sample_data_positive.txt"
    filename = "manhattan_tourist_sample_data_positive_2.txt"
    # filename = "rosalind_ba5b_2.txt"
    # filename = "manhattan_tourist_sample_data_negative.txt"
    m, n, down, right = parse_data(filename)

    len_longest_path = manhattan_tourist(m, n, down, right)
    print(f"Length of Longest path from node (0, 0) -> ({n}, {m}) is: ", end="")
    print(f"\033[92m\033[1m\033[4m{len_longest_path}\033[0m")
    print(f"\n(Duplicated output in case colored formatting of above output does not display porperly:")
    print(f"Length of Longest path from node (0, 0) -> ({n}, {m}) is: ", end="")
    print(f"{len_longest_path})")
