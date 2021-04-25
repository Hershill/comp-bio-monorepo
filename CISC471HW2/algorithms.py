"""
algorithms.py file containing the implementation of the algorithms for the 'Part 1 - Programming' section

Part 1 of HW 2 for CISC 471, Computational Biology.

By: Hershil Devnani (20001045)

There are two ways to run the program, outlined below. Both methods run the unittests.
Sample Usage:
  $ python -m unittest unittests.py
  $ python -m main algorithms.py

You can also modify the code in this file, main.py, in the main function to try different data sets and
parameters for each function by uncommenting and modifying the commented out function calls.
"""

import re
import random
import copy


def create_cycle(start, graph):
    """ Create an initial cycle from a given graph

    :param start: node to start cycle at
    :param graph: graph of nodes
    :return: a cycle starting at start using data from graph
    """

    cycle = list()
    unexplored_edges = copy.deepcopy(graph)
    next_node = next_edges(start, unexplored_edges)
    cycle.append(start)

    current_node = start
    # while unexplored edges in graph
    while next_node is not None:  # at worst complexity, will run for the number of edges in the graph
        next_node = next_node[random.randrange(len(next_node))]  # in the event there is more than one next node
        cycle.append(next_node)

        # remove edge, next combo from unexplored edges
        unexplored_edges = update_unexplored_edges(current_node, next_node, unexplored_edges)
        current_node = next_node

        # keep building graph
        next_node = next_edges(next_node, unexplored_edges)

    return cycle, unexplored_edges


def find_euler_cycle(cycle, unexplored_edges):
    """Finds the Eulerian cycle for a given starting cycle. Loops and finds the Eulerian cycle given a cycle
    that is no Eulerian but still has unexplored edges in it's graph.

    :param cycle: a cycle with unexplored edges
    :param unexplored_edges: the set of unexplored edges in the graph that cycle was based off of
    :return: the Eulerian cycle
    """
    
    # get unvisited nodes in graph
    graph_unexp = unexplored_in_graph(cycle, unexplored_edges)  # nodes that have an edge we have not traversed

    if graph_unexp == []:
        print("Already Eularian")
        return cycle

    count = 0
    while graph_unexp:  # reconstruct using new start node, runs at worst case the number of edges in the graph
        current_node = random.choice(graph_unexp)
        # continue building the graph
        new_cycle = cycle[cycle.index(current_node):-1] + cycle[:cycle.index(current_node)]

        # get new cycle, merge and check for more graph unexp
        new_cyc, unexplored_edges = create_cycle(current_node, unexplored_edges)
        cycle = new_cycle + new_cyc
        graph_unexp = unexplored_in_graph(cycle, unexplored_edges)
        count += 1

    return cycle


def unexplored_in_graph(cycle, unexplored):
    """Keeps track of the nodes in cycle that have edges we have not explored yet. Set of potential start node
    for the next iteration of finding the Eulerian cycle.

    :param cycle: a given cycle
    :param unexplored: all unexplored edges in the graph that the cycle was based off of
    :return: Set of potential start node for the next iteration of finding the Eulerian cycle.
    """

    unexp_graph = list()
    for node in cycle:
        if node in unexplored:
            unexp_graph.append(node)
    return unexp_graph


def start_node(graph):
    """ Select a random start node from a given graph.

    :param graph: a given graph
    :return: a random node in graph
    """

    start = random.choice(list(graph.values()))[0]
    return start


def next_edges(edge, graph):
    """ Get all the next edges that a node points tp

    :param edge: current node we are at
    :param graph: a given graph
    :return: the set of next nodes or neighbouring nodes
    """

    if edge in graph:
        return graph[edge]
    return None


def update_unexplored_edges(edge, next_vertex, unexplored_edges):
    """Update the set of unexplored edges in a graph while building a cycle

    :param edge: current node
    :param next_vertex: the next node
    :param unexplored_edges: dictionary of unexplored edges
    :return: update dictionary of unexplored edges with pair of (edge, next_vertex) removed
    """

    # grab the key values
    update = unexplored_edges[edge]

    for val in update:
        if val == next_vertex:
            update.remove(next_vertex)

    if not update:  # if update is empty after removing next node
        del unexplored_edges[edge]
    else:
        unexplored_edges[edge] = update

    return unexplored_edges


def parse_data(filename):
    """Read in the text file and build the graph

    :param filename: text files will graph data
    :return: graph in the form of a dictionary
    """

    # grab text file
    # parse via delimiter
    # store as dictionary, key = node, values = list of next nodes
    graph_edges = dict()

    with open(filename) as file:
        data_set = file.readlines()

    for line in data_set:
        # print(line.strip('\n').split(' -> '))
        edge = re.split(' -> |,', line.strip('\n'))  # parse line input into edge + next nodes
        edge = [int(i) for i in edge]  # convert data of each edge to integers
        # print(re.split(' -> |,', line.strip('\n')))
        graph_edges[edge[0]] = edge[1:]
    return graph_edges


def display_cycle(cycle):
    """Formats the output of a cycle in the required format

    :param cycle: a cycle
    :return: the cycle formatted as a string according to the required format
    """

    output_string = ''
    for i in range(len(cycle) - 1):
        output_string += str(cycle[i]) + "->"
    output_string += str(cycle[-1])
    return output_string


def main():
    filename = "sample_data.txt"
    # filename = "sample_data_2.txt"
    graph = parse_data(filename)
    # select random start node

    # get the first cycle
    cycle, unvisited = create_cycle(start_node(graph), graph)
    new_cycle = find_euler_cycle(cycle, unvisited)
    print(display_cycle(new_cycle))

    # [6, 8, 7, 9, 6, 5, 4, 2, 1, 0, 3, 2, 6]


if __name__ == '__main__':
    main()
