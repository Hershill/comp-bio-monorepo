"""
Part 2 of HW 6 for CISC 471, Computational Biology.
By: Hershil Devnani (20001045)


Lesson 5.13: Exercise Break:
Design a space-efficient algorithm for local sequence alignment.


# get local alignment
# isolate local alignment
# contain matrix of local alignment and apply space alignment


def alignmentPath(source, sink):
    middle_node = findMiddleNode
    alignmentPath(source, middle_node)
    alignmentPath(middle_node, sink)


def findMiddleNode(v, w, top, bottom, left, right):
    return highest scoring node in middle col of graph, index of middle node


def linearSpaceAlignment(v, w, top, bottom, left, right)
    if left == right
        # alignment of empty strings
        print path formed by bottom − top vertical edges
    if top == bottom
        print path formed by right − left horizontal edges

    middle_node, middle_index= findMiddleNode(v, w, top, bottom, left, right)

    if max score is 0 along any edge from top left node to middle node:
        ignore these nodes and remove them from teh graph

    linearSpaceAlignment(v, w, top, midNode, left, middle_index)
    print middle_node

    if middle_node = "↓" or middle_node ="↘"
        midNode = midNode + 1
    linearSpaceAlignment(v, w, midNode, bottom, middle_index, right)

"""