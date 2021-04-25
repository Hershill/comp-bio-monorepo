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

import re
import random
import copy


def create_cycle(start, graph):
    cycle = list()
    unexplored_edges = copy.deepcopy(graph)
    next_node = next_edges(start, unexplored_edges)
    cycle.append(start)
    # unexplored_edges = update_unexplored_edges(start, next_node, unexplored_edges)

    current_node = start
    # while unexplored edges in graph
    while next_node is not None:
        # if type(next_node) is list:
        next_node = next_node[random.randrange(len(next_node))]  # in the event there is more than one next node
        cycle.append(next_node)

        # remove edge, next combo from unexplored edges
        unexplored_edges = update_unexplored_edges(current_node, next_node, unexplored_edges)

        current_node = next_node

        # keep building graph
        next_node = next_edges(next_node, unexplored_edges)

    return cycle, unexplored_edges


def new_cylce(start, completed, unexplored):
    # add last and first node back to unexplored
    # given [6,8,7,9,6] continue the graph
    # add 9,6 back to unexplored
    # start becomes 6
    # append lists
    # return cycle, unexplored, repeat
    return


def merge_cycles(cycleA, cycleB):
    return


def find_euler_cycle(cycle, unexplored_edges):
    # get unvisited nodes in graph
    graph_unexp = unexplored_in_graph(cycle, unexplored_edges)  # nodes that have an edge we have not traversed

    if graph_unexp == []:
        print("Already eularian")
        return cycle

    # print(f"graph unexp: {graph_unexp}")
    # current_node = random.choice(graph_unexp)
    # cycle = list()
    # first = cycle.pop(0)
    # if first not in unexplored_edges:
    #     unexplored_edges[first] = current_node
    # else:
    #     unexplored_edges[first] = unexplored_edges[first].append(current_node)

    # next_node = unexplored_edges[current_node]
    # print(f"new start node: {current_node}")
    # print(f"next of start node: {next_node}")
    count = 0
    while graph_unexp:
        print(f"graph unexp: {graph_unexp}")
        current_node = random.choice(graph_unexp)
        print(f"current node: {current_node}")
        # continue building the graph
        # add back unexplored edges
        # list of unexplored edges to add back
        # add_back = cycle[:cycle.index(current_node)]
        # add_back = cycle[0]
        # print(f"add back: {add_back}")
        new_cycle = cycle[cycle.index(current_node):-1] + cycle[:cycle.index(current_node)]

        print(f"new cycle: {new_cycle}")

        # for node in add_back:
        # if add_back not in unexplored_edges:
        #     unexplored_edges[add_back] = [current_node]
        # else:
        #     unexplored_edges[add_back] = unexplored_edges[add_back].append(current_node)

        # get new cycle, merge and check for more graph unexp
        new_cyc, unexplored_edges = create_cycle(current_node, unexplored_edges)
        print(f"append cycle: {new_cyc}")
        cycle = new_cycle + new_cyc
        print(f"next cycle: {cycle}")
        graph_unexp = unexplored_in_graph(cycle, unexplored_edges)
        print(f"graph unexp end of loop: {graph_unexp}")
        count += 1
        # if count == 5:
        #     break
        # break

    # while graph_unexp:
    #     # select random start node from nodes with unexplored in current graph
    #     # current_node = random.choice(graph_unexp)
    #     # start new cycle traversal from selected node
    #     # print(f"new start node: {current_node}")
    #     # create_cycle(current_node, )
    #     # break
    #     # graph_unexp = unexplored_in_graph(cycle, unexplored_edges)
    #     new_cyc, unexp = create_cycle(current_node, unexplored_edges)
    #     new_cyc = cycle + new_cyc
    #     graph_unexp = unexplored_in_graph(new_cyc, unexp)  # nodes that have an edge we have not traversed
    #     print(f"next graph unexp: {graph_unexp}")
    #     break

    return cycle


def unexplored_in_graph(cycle, unexplored):
    unexp_graph = list()
    for node in cycle:
        if node in unexplored:
            unexp_graph.append(node)
    return unexp_graph


def start_node(graph):
    start = random.choice(list(graph.values()))[0]
    return start


def next_edges(edge, graph):
    if edge in graph:
        return graph[edge]
    return None


def eulerian_cycle(graph_edges):
    cycle = list()
    unexplored_edges = copy.deepcopy(graph_edges)

    # unexplored edges in graph
    # take current cycles, for node in current graph, if node and next in unused edges select new unused edge

    # select random key
    # select random number in range of len of keys
    current_edge = random.choice(list(graph_edges.values()))[0]
    # random_edge_index = random.randrange(len(graph_edges[random_start[0]]))
    # next_edge = graph_edges[random_start[0]][random_edge_index]

    cycle.append(current_edge)
    # cycle.append(next_edge)  # break into function to get set of next nodes

    next_vertex = next_edge(current_edge, unexplored_edges)

    while next_vertex is not None:
        cycle.append(next_vertex)
        unexplored_edges = update_unexplored_edges(current_edge, next_vertex, unexplored_edges)
        current_edge = next_vertex
        peek = next_edge(current_edge, unexplored_edges)  # used to look ahead
        if peek is None and unexplored_edges != {}:
            # next vertex is none, i.e. has already been explored
            # and we have remaining unexplored edges

            for i in cycle:
                next_i = next_edge(i, unexplored_edges)
                if next_i is not None:
                    current_edge = next_i

            cycle = list()
            # current_edge = random.choice(list(unexplored_edges.values()))[0]
            cycle.append(current_edge)
            next_vertex = next_edge(current_edge, unexplored_edges)
            unexplored_edges = copy.deepcopy(graph_edges)
        elif peek is None and unexplored_edges == {}:
            # next vertex is none, i.e. has already been explored
            # and we have explored all edges
            return cycle
        else:
            next_vertex = peek

    # while
    # cycle = [random.choice(list(unexplored_edges.values()))[0]]
    # next_vertex = cycle[0]

    # next_vertex = next_edge(current_edge, unexplored_edges)

    # break

    # if not unexplored_edges:  # if unexplored_edges is empty
    #     return cycle
    # else:

    # while not unexplored_edges:  # while unexplored edges is not empty

    # get next edges of next edge and repeat until graph is formed and next is start edge
    # remove from unexplored_edges dict
    # first cycle found

    # loop, select edge in used edges, reset unused edges
    # continue to find cycles while we have unused edges

    return cycle


def next_edge(edge, unexplored_edges):
    if edge in unexplored_edges:
        # random_edge_index = random.randrange(len(unexplored_edges[edge]))
        random_edge_index = random.randrange(len(unexplored_edges[edge]))
        next_vertex = unexplored_edges[edge][random_edge_index]
        return next_vertex
    return None


def update_unexplored_edges(edge, next_vertex, unexplored_edges):
    # grab the key values
    update = unexplored_edges[edge]

    for val in update:
        if val == next_vertex:
            # remove = next_vertex
            # update = update.remove(next_vertex)
            update.remove(next_vertex)

    if not update:  # if update is empty after removing next node
        del unexplored_edges[edge]
    else:
        unexplored_edges[edge] = update

    return unexplored_edges


def parse_data(filename):
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
    output_string = ''
    for i in range(len(cycle) - 1):
        output_string += str(cycle[i]) + "->"
    output_string += str(cycle[-1])
    return output_string


def main():
    # filename = "sample_data.txt"
    filename = "sample_data_2.txt"
    # filename = "rosalind_ba3f.txt"
    # filename = "rosalind_ba3f_2.txt"
    graph = parse_data(filename)
    # print(graph)
    # cycle = eulerian_cycle(graph)
    # print(cycle)
    # print(display_cycle(cycle))

    # select random start node

    # get the first cycle
    cycle, unvisited = create_cycle(start_node(graph), graph)
    print(cycle)
    print(unvisited)
    new_cycle = find_euler_cycle(cycle, unvisited)
    print(f"Eucledian graph: {new_cycle}")
    # print(create_cycle(start_node(graph), graph))

    print(display_cycle(new_cycle))

    # [6, 8, 7, 9, 6, 5, 4, 2, 1, 0, 3, 2, 6]


if __name__ == '__main__':
    main()
