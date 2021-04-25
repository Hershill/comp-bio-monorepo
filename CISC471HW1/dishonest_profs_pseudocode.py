"""
dishonest_profs_pseudocode.py file that contains the pseudocode of the algorithm for '2.3 - Dishonest Professors'

Part 2 of HW 1 for CISC 471, Computational Biology.

By: Hershil Devnani (20001045)

Sample Usage:
  $ python indices.py

Please see the accompanied file, hw1.pdf, for an explanation of the pseudocode.
"""

# 2.3 Dishonest Professors pseudocode


"""
for all professors:
    selected_prof <- select a random professor

l <- number of maximum dishonest professors (in this case 49)

# initialize the counters
num_honest <- 0
num_dishonest <- 0
set_profs <- set of every professor

while l != 0:
    for x in (set_profs - selected_prof):
        response <- ask x if the selected_prof is honest, true if yes, otherwise no
        if response:
            num_honest <- +=1 if x replies yes
        else:
            num_dishonest <- += 1

        if dishonest > honest:
            we know that the at least num_dishonest or num_honest + 1 number are lying
            l <- l - dishonest
            set_profs <- set of all profs not including those professors asked in the previous round
        else if honest == l:
            selected_prof is an honest prof

ask selected_prof a question about whether each of the remaining 99 profs is honest or dishonest
based on selected_prof's answer, append to separate lists to separate honest and dishonest profs
"""
