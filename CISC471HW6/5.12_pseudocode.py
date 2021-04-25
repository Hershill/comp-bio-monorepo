"""
Part 2 of HW 6 for CISC 471, Computational Biology.
By: Hershil Devnani (20001045)


Lesson 5.12: Exercise Break:
Design an algorithm for computing optimal local (rather than global)
alignment with affine gap penalties.


def longest_path_sequence(v, w, score, global, local, affine):
    # get the sequence of the longest path for v an w

    backtrack = nested list

    # init array for s, the longest path to any node s[i][j]
    for i in range(0, len(v)):
        s[i][0] = 0

    for i in range(0, len(w)):
        s[0][j] = 0

    for i in range(1, len(v)):
        for i in range(1, len(w)):
            match = 0
            if v[i-1] = w[j-1]:
                match = 1

            # s[i][j] = max(s[i-1][j], s[i][j-1], s[i-1][j-1] + match)

            # sigma = indel penalty
            # mew = mismatch penalty
            # epsilon = gap penalty
            sigma, mew, epsilon <- score  # penalty values derived from a given score matrix

            if global:
                if affine:
                    lower[i][j] = max(
                                    lower[i-1][j] - epsilon,
                                    middle[i-1][j] - sigma
                                    )

                    upper[i][j] = max(
                                    upper[i][j-1] - epsilon,
                                    middle[i][j-1] - sigma
                                    )

                    middle[i][j] = max(
                                    lower[i][j],
                                    middle[i-1][j-1] + score(v[i], w[i])
                                    upper[i][j]
                                    )
                else:
                    s[i][j] = max(
                                s[i-1][j] - sigma,
                                s[i][j-1] - sigma,
                                s[i-1][j-1] + match if v[i] == w[j],
                                s[i-1][j-1] - mew if v[i] != w[j])

            else if local:
                if affine:
                    # Pseudocode for global alignment with affine gap penalties

                    # lower represents insertions and contains edges that are transversed vertically
                    lower[i][j] = max(
                                    0,
                                    lower[i-1][j] - epsilon,
                                    middle[i-1][j] - sigma
                                    )

                    # upper represents deletions and contains edges that are transversed horizontally
                    upper[i][j] = max(
                                    0,
                                    upper[i][j-1] - epsilon,
                                    middle[i][j-1] - sigma
                                    )

                    # middle represents matches or mismatches and contains edges that are transversed diagonally
                    middle[i][j] = max(
                                    0,
                                    lower[i][j],
                                    middle[i-1][j-1] + score(v[i], w[i])
                                    upper[i][j]
                                    )

                    # if affine store tuple for which of the three matrices the sequence was chosen from
                    s[i][j] <- select sequence from maximum combination of going from middle to lower to upper nodes

                else:
                    # for local alignment
                    s[i][j] = max(
                                0,  # weight of newly added edge from (0, 0) to (i, j)
                                s[i-1][j] + score(v[i], -),  # score to move down
                                s[i][j-1] + score(-, w[j]),  # score to move right
                                s[i-1][j-1] + score(v[i], w[j])  # score to move diagonal
                                )

            # backtrack to get alignment sequence
            # if affine, using s[i][j] to work backwards between the 3 matrices upper, middle, lower
            if s[i][j] = s[i-1][j]:
                backtrack[i][j] = down
            else if s[i][j] = s[i][j-1]:
                backtrack[i][j] = right
            else if s[i][j] = s[i-1][j-1] + match:
                backtrack[i][j] = diagonal

    if local:
        score = max of s[i][j] (for 0 <= i <=n and 0 <= j <= n)
    else:
        score = s[n][m]

    return backtrack, score


def output_lcs(backtrack_sequence, v, i, j):
    if i == 0 or j -- 0:
        return ""

    if backtrack[i][j] == down:
        return output_lcs(backtrack_sequence, v, i - 1, j)

    else if backtrack[i][j] == right:
        return output_lcs(backtrack_sequence, v, i, j - 1)

    else:
        return output_lcs(backtrack_sequence, v, i - 1, j - 1) + v[i]


def global_alignment(v, w, score_matrix):
    # find the longest path in alignment graph
    backtrack_sequence, score = longest_path_sequence(v, w, score_matrix, true, false)
    output = output_lcs(backtrack_sequence, v, len(v) ,len(w)
    print(output)

    return output, score


# Pseudocode for local alignment
def local_alignment(v, w, score_matrix):
    # setting affine selection to true for local affine alignment
    backtrack_sequence, score = longest_path_sequence(v, w, score_matrix, false, true, true)
    output = output_lcs(backtrack_sequence, v, len(v) ,len(w)
    print(output)

    return output, score

def main():
    local_alignment(v, w, score_matrix)

main()

"""
