"""
algorithms.py file containing the implementation of the algorithms for the 'Part 1 - Programming' section

Part 1 of HW 2 for CISC 471, Computational Biology.

By: Hershil Devnani (20001045)

There are two ways to run the program, outlined below. Both methods run the unittests.
Sample Usage:
  $ python -m unittest unittests.py
  $ python -m main contigs.py

You can also modify the code in this file, main.py, in the main function to try different data sets and
parameters for each function by uncommenting and modifying the commented out function calls.
"""


def build_graph(filename):
    """Build graph from text file data

    :param filename: file to read from
    :return: built graph as a dictionary
    """
    graph = dict()

    with open(filename) as file:
        data_set = file.readlines()
        for line in data_set:
            line = line.strip('\n')
            # graph[line] = [prefix(line), suffix(line)]
            if prefix(line) not in graph:
                graph[prefix(line)] = list()
            graph[prefix(line)].append(suffix(line))
            # graph[prefix(line)] = suffix(line)

    return graph


def glue_nodes(graph):
    """Glue nodes together to find contigs

    :param graph: graph data
    :return: graph with overlapping sequences glued together
    """
    # glued_graph = copy.deepcopy(graph)
    glued_graph = dict()
    keys_del = list()

    # convert all keys to lists
    for edge in graph:
        graph[edge] = [graph[edge]]

    for pfx in graph:
        sfx_list = graph[pfx]
        for sfx in sfx_list:
            if sfx in graph:
                # if the suffix of the current node matches the prefix of another node
                glued_graph[pfx] = graph[pfx] + graph[sfx]
                keys_del.append(sfx)

    return glued_graph


def get_contigs(graph):
    """Find the contigs in the graph

    :param graph: graph data
    :return: all found contigs
    """
    # return list of maximal non-branching paths
    contigs = []
    del_list = []

    for node in graph:
        if inputs(node, graph) != 1 or outputs(node, graph) != 1:
            if outputs(node, graph) > 0:
                for outs in graph[node]:
                    nbp = []
                    nbp.append(node)
                    nbp.append(outs)
                    non_branching_path = node + outs[-1]
                    del_list.append(node)
                    w = outs
                    while inputs(w, graph) == 1 and outputs(w, graph) ==1:
                        for val in graph[w]:
                            nbp.append(val)
                            non_branching_path = non_branching_path + val[-1]
                            del_list.append(w)
                            w = val
                    contigs.append(non_branching_path)

    for dels in del_list:
        if dels in graph:
            del graph[dels]

    return contigs


def isolated_cycle(graph):
    """Find an isolated cycle, where all nodes are 1-in-1-out in the graph, which can be identified as contigs

    :param graph: graph data
    :return: isolated cycles as contigs
    """
    all_iso = list()
    for node in graph:
        iso_list = list()
        next_node = node
        if len(next_node) < 1:
            next_node = next_node[0]
            while outputs(next_node, graph) == 1:
                iso_list += graph[node]
                next_node = graph[node]
            iso_list += next_node
            all_iso.append(iso_list)

    return all_iso


def inputs(node, graph):
    """Counts the number of incoming edges of a node

    :param node: selected node
    :param graph: graph data
    :return: number  of incoming edges of node
    """
    count = 0
    for i in graph:
        if node in graph[i]:
            count += 1
    return count


def outputs(node, graph):
    """Counts the number of outgoing edges in a graph

    :param node: selected node
    :param graph: graph data
    :return: number of outgoing edges of node
    """
    if node in graph:
        return len(graph[node])
    return 0


def prefix(kmer):
    """Get the prefix of a k-mer

    :param kmer: selected k-mer
    :return: prefix of k-mer
    """
    return kmer[:len(kmer)-1]


def suffix(kmer):
    """Get the suffix of a k-mer

    :param kmer: selected k-mer
    :return: suffix of k-mer
    """
    return kmer[1:]


def main():
    filename = "unittest_generate_contigs_positive_b.txt"
    graph = build_graph(filename)
    contigs = get_contigs(graph)
    print(contigs)


if __name__ == '__main__':
    main()
