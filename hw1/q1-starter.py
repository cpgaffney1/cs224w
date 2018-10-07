################################################################################
# CS 224W (Fall 2018) - HW1
# Starter code for Problem 1
# Author: praty@stanford.edu
# Last Updated: Sep 27, 2018
################################################################################

import snap
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

# Setup
erdosRenyi = None
smallWorld = None
collabNet = None


# Problem 1.1
def genErdosRenyi(N=5242, E=14484):
    """
    :param - N: number of nodes
    :param - E: number of edges

    return type: snap.PUNGraph
    return: Erdos-Renyi graph with N nodes and E edges
    """
    ############################################################################
    # TODO: Your code here!
    pairs = []
    Graph = snap.TUNGraph.New()
    for i in range(N):
        Graph.AddNode(i)
        for j in range(i + 1, N):
            pairs.append((i, j))
    idx = np.array(np.random.choice(len(pairs), E, replace=False))
    pairs = np.array(pairs)
    chosen_pairs = pairs[idx]
    
    for pair in chosen_pairs:
        Graph.AddEdge(pair[0], pair[1])
        assert(Graph.IsEdge(pair[1], pair[0]))

    assert(Graph.GetEdges() == E)
    assert(Graph.GetNodes() == N)
    ############################################################################
    return Graph


def genCircle(N=5242):
    """
    :param - N: number of nodes

    return type: snap.PUNGraph
    return: Circle graph with N nodes and N edges. Imagine the nodes form a
        circle and each node is connected to its two direct neighbors.
    """
    ############################################################################
    # TODO: Your code here!
    Graph = snap.TUNGraph.New()
    Graph.AddNode(0)
    for i in range(1, N):
        Graph.AddNode(i)
        Graph.AddEdge(i-1, i)
    Graph.AddEdge(N-1, 0)
    assert(Graph.IsEdge(0, 1) and Graph.IsEdge(1, 0))
    for i in range(N):
        assert(Graph.GetNI(i).GetDeg() == 2)
    ############################################################################
    return Graph


def connectNbrOfNbr(Graph, N=5242):
    """
    :param - Graph: snap.PUNGraph object representing a circle graph on N nodes
    :param - N: number of nodes

    return type: snap.PUNGraph
    return: Graph object with additional N edges added by connecting each node
        to the neighbors of its neighbors
    """
    ############################################################################
    # TODO: Your code here!
    pairs = []
    for i in range(N):
        Graph.AddEdge(i, (i + 2) % N)
    for i in range(N):
        assert(Graph.GetNI(i).GetDeg() == 4)
    ############################################################################
    return Graph


def connectRandomNodes(Graph, M=4000):
    """
    :param - Graph: snap.PUNGraph object representing an undirected graph
    :param - M: number of edges to be added

    return type: snap.PUNGraph
    return: Graph object with additional M edges added by connecting M randomly
        selected pairs of nodes not already connected.
    """
    ############################################################################
    # TODO: Your code here!
    pairs = []
    for i in range(Graph.GetNodes()):
        for j in range(i+1, Graph.GetNodes()):
            if not Graph.IsEdge(i, j):
                pairs.append((i, j))
    idx = np.array(np.random.choice(len(pairs), M, replace=False))
    pairs = np.array(pairs)
    chosen_pairs = pairs[idx]
    for pair in chosen_pairs:
        Graph.AddEdge(pair[0], pair[1])
    ############################################################################
    return Graph


def genSmallWorld(N=5242, E=14484):
    """
    :param - N: number of nodes
    :param - E: number of edges

    return type: snap.PUNGraph
    return: Small-World graph with N nodes and E edges
    """
    Graph = genCircle(N)
    Graph = connectNbrOfNbr(Graph, N)
    Graph = connectRandomNodes(Graph, 4000)
    return Graph


def loadCollabNet(path):
    """
    :param - path: path to edge list file

    return type: snap.PUNGraph
    return: Graph loaded from edge list at `path and self edges removed

    Do not forget to remove the self edges!
    """
    ############################################################################
    # TODO: Your code here!
    Graph = snap.LoadEdgeList(snap.PUNGraph, path, 0, 1)
    ############################################################################
    return Graph


def getDataPointsToPlot(Graph):
    """
    :param - Graph: snap.PUNGraph object representing an undirected graph

    return values:
    X: list of degrees
    Y: list of frequencies: Y[i] = fraction of nodes with degree X[i]
    """
    ############################################################################
    # TODO: Your code here!
    X, Y = [], []
    degree_map = {}
    for node in Graph.Nodes():
        if node.GetDeg() not in degree_map.keys():
            degree_map[node.GetDeg()] = 0
        degree_map[node.GetDeg()] += 1
    for deg in degree_map.keys():
        X.append(deg)
        Y.append(float(degree_map[deg]) / Graph.GetNodes())
    print len(X)
    ############################################################################
    return X, Y


def Q1_1():
    """
    Code for HW1 Q1.1
    """
    global erdosRenyi, smallWorld, collabNet
    erdosRenyi = genErdosRenyi(5242, 14484)
    smallWorld = genSmallWorld(5242, 14484)
    collabNet = loadCollabNet("ca-GrQc.txt")
    assert(smallWorld is not None)
    assert(erdosRenyi is not None)
    x_erdosRenyi, y_erdosRenyi = getDataPointsToPlot(erdosRenyi)
    plt.loglog(x_erdosRenyi, y_erdosRenyi, color = 'y', label = 'Erdos Renyi Network')

    x_smallWorld, y_smallWorld = getDataPointsToPlot(smallWorld)
    plt.loglog(x_smallWorld, y_smallWorld, linestyle = 'dashed', color = 'r', label = 'Small World Network')

    x_collabNet, y_collabNet = getDataPointsToPlot(collabNet)
    plt.loglog(x_collabNet, y_collabNet, linestyle = 'dotted', color = 'b', label = 'Collaboration Network')

    plt.xlabel('Node Degree (log)')
    plt.ylabel('Proportion of Nodes with a Given Degree (log)')
    plt.title('Degree Distribution of Erdos Renyi, Small World, and Collaboration Networks')
    plt.legend()
    plt.savefig('q1_plots.png')


# Execute code for Q1.1
Q1_1()


# Problem 1.2 - Clustering Coefficient

def calcClusteringCoefficientSingleNode(Node, Graph):
    """
    :param - Node: node from snap.PUNGraph object. Graph.Nodes() will give an
                   iterable of nodes in a graph
    :param - Graph: snap.PUNGraph object representing an undirected graph

    return type: float
    returns: local clustering coeffient of Node
    """
    ############################################################################
    # TODO: Your code here!
    C = 0.0
    if Node.GetDeg() < 2:
        return C
    ei = 0.0
    for i in range(Node.GetDeg()):
        for j in range(i + 1, Node.GetDeg()):
            if Graph.IsEdge(Node.GetNbrNId(i), Node.GetNbrNId(j)):
                ei += 1
    
    assert(ei <= (Node.GetDeg() * (Node.GetDeg() - 1)) / 2)
    C = 2.0 * ei / (Node.GetDeg() * (Node.GetDeg() - 1))
    ############################################################################
    return C

def calcClusteringCoefficient(Graph):
    """
    :param - Graph: snap.PUNGraph object representing an undirected graph

    return type: float
    returns: clustering coeffient of Graph
    """
    ############################################################################
    # TODO: Your code here! If you filled out calcClusteringCoefficientSingleNode,
    #       you'll probably want to call it in a loop here
    C = 0.0
    for node in Graph.Nodes():
        Ci = calcClusteringCoefficientSingleNode(node, Graph)
        C += Ci
    C /= float(Graph.GetNodes())
    ############################################################################
    return C

def Q1_2():
    """
    Code for Q1.2
    """
    C_erdosRenyi = calcClusteringCoefficient(erdosRenyi)
    C_smallWorld = calcClusteringCoefficient(smallWorld)
    C_collabNet = calcClusteringCoefficient(collabNet)

    print('Clustering Coefficient for Erdos Renyi Network: %f' % C_erdosRenyi)
    print('Clustering Coefficient for Small World Network: %f' % C_smallWorld)
    print('Clustering Coefficient for Collaboration Network: %f' % C_collabNet)


# Execute code for Q1.2
Q1_2()
