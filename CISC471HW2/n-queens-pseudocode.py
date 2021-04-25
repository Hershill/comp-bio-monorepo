
"""
n-queens-pseudocode.py file that contains the pseudocode of the algorithm for '2.3 - Dishonest Professors'

Part 2 of HW 2 for CISC 471, Computational Biology.

By: Hershil Devnani (20001045)

Please see the accompanied file, hw2.pdf, for an explanation of the pseudocode.

2. Write a recursive algorithm that either places the n Queen’s or determines that no such placement is possible.

def place_queens(n):
    if n = 2 or n = 3:
        return no solution

    places the first queen at a random position in column 1
    for queen in range(len(n) - 1):
        current_col <- select column n+1
        calculate the number of conflicts in each position in current_col
        place queen in column with least conflicts

    while

def place_queens(n, col, queens_placement):

    if n = 2 or n = 3:
        return false

    if col == n:
        # we have placed every queen peacefully with no conflict
        return true


    for i in range(n):
        # for loop to iterate over possible row positions in current column

        current_row <- i
        queens_placement[col] = current_row  # keep track of where the queen is in the column
        conflict = check_conflicts()

        # if the selected position has a conflict, keep looping through each row until a conflict free row is found
        # in the selected column
        if not conflict:
            # if this returns false, our next queen was not peacefully palced, so we run the for loop again and change the
            # position of the previously placed queen
            next_queen = queens(n, col+1, queens_placement)

            if next_queen:
                return true

    return false  # all possible placements return a conflict


def check_conflicts():
    return true if conflict
    return false if no conflict

def main():
    n = 4  # for example
    queens_placement = -1*n
    place_queens(4,0,queens_placement)

3. Modify the algorithm so that it counts all peaceful placements.


def place_queens(n, col, queens_placement):

    if n = 2 or n = 3:
        return false

    if col == n:
        # we have placed every queen peacefully with no conflict
        return true


    for i in range(n):
        # for loop to iterate over possible row positions in current column

        current_row <- i
        queens_placement[col] = current_row  # keep track of where the queen is in the column
        conflict = check_conflicts()

        # if the selected position has a conflict, keep looping through each row until a conflict free row is found
        # in the selected column
        if not conflict or queens_placement[col] = -1:  # also check if placement is set to -1, which means iterate and find another solution
            # if this returns false, our next queen was not peacefully palced, so we run the for loop again and change the
            # position of the previously placed queen
            next_queen = queens(n, col+1, queens_placement)

            if next_queen:
                # go back and adjust the placement of the rows by resetting the final position of the queen
                queens_placement[col] = -1
                return true

    return false  # all possible placements return a conflict


def check_conflicts():
    return true if conflict
    return false if no conflict

def main():
    n = 4  # for example
    queens_placement = -1*n
    place_queens(4,0,queens_placement)

"""