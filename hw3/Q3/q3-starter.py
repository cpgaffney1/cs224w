###############################################################################
# CS 224W (Fall 2018) - HW2
# Starter code for Problem 3
###############################################################################

import snap
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np

# Setup
num_voters = 10000
decision_period = 10

def get_neighbors(graph, nodeId):
    n = set([])
    for i in range(graph.GetNI(nodeId).GetDeg()):
        nid = graph.GetNI(nodeId).GetNbrNId(i)
        n.add(nid)
    return n


def read_graphs(path1, path2):
    """
    :param - path1: path to edge list file for graph 1
    :param - path2: path to edge list file for graph 2

    return type: snap.PUNGraph, snap.PUNGraph
    return: Graph 1, Graph 2
    """
    ###########################################################################
    # TODO: Your code here!
    ###########################################################################
    Graph1 = snap.LoadEdgeList(snap.PUNGraph, path1, 0, 1)
    Graph2 = snap.LoadEdgeList(snap.PUNGraph, path2, 0, 1)
    return Graph1, Graph2


def initial_voting_state(Graph):
    """
    Function to initialize the voting preferences.

    :param - Graph: snap.PUNGraph object representing an undirected graph

    return type: Python dictionary
    return: Dictionary mapping node IDs to initial voter preference
            ('A', 'B', or 'U')

    Note: 'U' denotes undecided voting preference.

    Example: Some random key-value pairs of the dict are
             {0 : 'A', 24 : 'B', 118 : 'U'}.
    """
    voter_prefs = {}
    ###########################################################################
    # TODO: Your code here!
    ###########################################################################
    for node in Graph.Nodes():
        support = ''
        if node.GetId() % 10 <= 3:
            support = 'A'
        elif node.GetId() % 10 >= 4 and node.GetId() % 10 <= 7:
            support = 'B'
        else:
            support = 'U'
        voter_prefs[node.GetId()] = support
    
    assert(len(voter_prefs) == num_voters)

    return voter_prefs


def iterate_voting(Graph, init_conf):
    """
    Function to perform the 10-day decision process.

    :param - Graph: snap.PUNGraph object representing an undirected graph
    :param - init_conf: Dictionary object containing the initial voting
                        preferences (before any iteration of the decision
                        process)

    return type: Python dictionary
    return: Dictionary containing the voting preferences (mapping node IDs to
            'A','B' or 'U') after the decision process.

    Hint: Use global variables num_voters and decision_period to iterate.
    """
    curr_conf = init_conf.copy()
    curr_alternating_vote = 'A'
    ###########################################################################
    # TODO: Your code here!
    
    for i in range(decision_period):
        for vid in range(num_voters):
            if curr_conf[vid] == 'U':
                neighbors = get_neighbors(Graph, vid)
                acount = get_vote_count(neighbors, curr_conf, 'A')
                bcount = get_vote_count(neighbors, curr_conf, 'B')
                if acount > bcount:
                    curr_conf[vid] = 'A'
                elif acount < bcount:
                    curr_conf[vid] = 'B'
                else:
                    curr_conf[vid] = curr_alternating_vote
                    curr_alternating_vote = 'A' if curr_alternating_vote == 'B' else 'B'
                

    ###########################################################################
    return curr_conf


def get_vote_count(nodeset, prefs, letter):
        count = 0
        for nid in nodeset:
            if prefs[nid] == letter:
                count += 1
        return count



def sim_election(Graph):
    """
    Function to simulate the election process, takes the Graph as input and
    gives the final voting preferences (dictionary) as output.
    """
    init_conf = initial_voting_state(Graph)
    conf = iterate_voting(Graph, init_conf)
    return conf


def winner(conf):
    """
    Function to get the winner of election process.
    :param - conf: Dictionary object mapping node ids to the voting preferences

    return type: char, int
    return: Return candidate ('A','B') followed by the number of votes by which
            the candidate wins.
            If there is a tie, return 'U', 0
    """
    ###########################################################################
    # TODO: Your code here!
    acount = get_vote_count(conf.keys(), conf, 'A')
    bcount = get_vote_count(conf.keys(), conf, 'B')
    if acount == bcount:
        winner = 'U'
    elif acount > bcount:
        winner = 'A'
    elif acount < bcount:
        winner = 'B'
    else:
        assert(False)
    margin = acount - bcount
    ###########################################################################
    return winner, margin

def Q1():
    print ("\nQ1:")
    Gs = read_graphs('graph1.txt', 'graph2.txt')    # List of graphs

    # Simulate election process for both graphs to get final voting preference
    final_confs = [sim_election(G) for G in Gs]

    # Get the winner of the election, and the difference in votes for both
    # graphs
    res = [winner(conf) for conf in final_confs]

    for i in xrange(2):
        print "In graph %d, candidate %s wins by %d votes" % (
                i+1, res[i][0], res[i][1]
        )


def Q2sim(Graph, k):
    """
    Function to simulate the effect of advertising.
    :param - Graph: snap.PUNGraph object representing an undirected graph
             k: amount to be spent on advertising

    return type: int
    return: The number of votes by which A wins (or loses), i.e. (number of
            votes of A - number of votes of B)

    Hint: Feel free to use initial_voting_state and iterate_voting functions.
    """
    ###########################################################################
    # TODO: Your code here!
    initial_prefs = initial_voting_state(Graph)
    for i in range(3000, 3000 + k / 100 - 1 + 1):
        initial_prefs[i] = 'A'
    
    prefs = iterate_voting(Graph, initial_prefs)
    winner_letter, margin = winner(prefs)

#    print(winner_letter, margin)

    return margin
    ###########################################################################


def find_min_k(diffs):
    """
    Function to return the minimum amount needed for A to win
    :param - diff: list of (k, diff), where diff is the value by which A wins
                   (or loses) i.e. (A-B), for that k.

    return type: int
    return: The minimum amount needed for A to win
    """
    ###########################################################################
    # TODO: Your code here!
    k_prev = float('-inf')
    for k, diff in diffs:
        assert(k > k_prev)
        k_prev = k
    for k, diff in diffs:
        if diff > 0:
            return k


    ###########################################################################


def makePlot(res, title, name):
    """
    Function to plot the amount spent and the number of votes the candidate
    wins by
    :param - res: The list of 2 sublists for 2 graphs. Each sublist is a list
                  of (k, diff) pair, where k is the amount spent, and diff is
                  the difference in votes (A-B).
             title: The title of the plot
    """
    Ks = [[k for k, diff in sub] for sub in res]
    res = [[diff for k, diff in sub] for sub in res]
    ###########################################################################
    # TODO: Your code here!
    ###########################################################################
    print Ks
    print res
    plt.plot(Ks[0], [0.0] * len(Ks[0]), ':', color='black')
    plt.plot(Ks[0], res[0], color='green', label='Graph 1')
    plt.plot(Ks[1], res[1], color='red', label='Graph 2')
    plt.xlabel('Amount spent ($)')
    plt.ylabel('#votes for A - #votes for B')
    plt.title(title)
    plt.legend()
    plt.savefig('{}.png'.format(name))
    plt.clf()


def Q2():
    print ("\nQ2:")
    # List of graphs
    Gs = read_graphs('graph1.txt', 'graph2.txt')

    # List of amount of $ spent
    Ks = [x * 1000 for x in range(1, 10)]

    # List of (List of diff in votes (A-B)) for both graphs
    res = [[(k, Q2sim(G, k)) for k in Ks] for G in Gs]

    # List of minimum amount needed for both graphs
    min_k = [find_min_k(diff) for diff in res]

    formatString = "On graph {}, the minimum amount you can spend to win is {}"
    for i in xrange(2):
        print formatString.format(i + 1, min_k[i])

    makePlot(res, 'TV Advertising', 'tv_ad_q2')


def get_degs(graph):
    degs = []
    ids = []
    for node in graph.Nodes():
        ids.append(node.GetId())
        degs.append(node.GetDeg())
    ids = np.array(ids)
    degs = np.array(degs)
    indices = np.flip(np.argsort(ids))
    ids = ids[indices]
    degs = degs[indices]
    for i in range(1, len(ids)):
        assert(ids[i] <= ids[i-1])
  
    return ids, degs

def Q3sim(Graph, k):
    """
    Function to simulate the effect of a dining event.
    :param - Graph: snap.PUNGraph object representing an undirected graph
             k: amount to be spent on the dining event

    return type: int
    return: The number of votes by which A wins (or loses), i.e. (number of
            votes of A - number of votes of B)

    Hint: Feel free to use initial_voting_state and iterate_voting functions.
    """
    ###########################################################################
    # TODO: Your code here!
    ids, degs = get_degs(Graph)
    ids = np.array(ids)
    degs = np.array(degs)
    indices = np.argsort(degs, kind='stable')
    degs = degs[indices]
    ids = ids[indices]
    for i in range(1, len(degs)):
        assert(degs[i] >= degs[i-1])
        if degs[i] == degs[i-1]:
            assert(ids[i] < ids[i-1])
    
    initial_prefs = initial_voting_state(Graph)

    n_high_rollers = k / 1000
    high_rollers = ids[-n_high_rollers:]
    for vid in high_rollers:
        initial_prefs[vid] = 'A'
    
    prefs = iterate_voting(Graph, initial_prefs)
    winner_letter, margin = winner(prefs)

    return margin
    ###########################################################################


def Q3():
    print ("\nQ3:")
    # List of graphs
    Gs = read_graphs('graph1.txt', 'graph2.txt')

    # List of amount of $ spent
    Ks = [x * 1000 for x in range(1, 10)]

    # List of (List of diff in votes (A-B)) for both graphs
    res = [[(k, Q3sim(G, k)) for k in Ks] for G in Gs]

    # List of minimum amount needed for both graphs
    min_k = [find_min_k(diff) for diff in res]

    formatString = "On graph {}, the minimum amount you can spend to win is {}"
    for i in xrange(2):
        print formatString.format(i + 1, min_k[i])

    makePlot(res, 'Wining and Dining', 'wine_dine_q3')


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

    ############################################################################
    return X, Y



   
def Q4():
    """
    Function to plot the distributions of two given graphs on a log-log scale.
    """
    print ("\nQ4:")
    ###########################################################################
    # TODO: Your code here!
    
    Gs = read_graphs('graph1.txt', 'graph2.txt')
    
    x1, y1 = getDataPointsToPlot(Gs[0]) 
    plt.loglog(x1, y1, color = 'b', label = 'Graph 1')
    
    x2, y2 = getDataPointsToPlot(Gs[1]) 
    plt.loglog(x2, y2, color = 'r', label = 'Graph 2')

    plt.xlabel('Node Degree (log)')
    plt.ylabel('Proportion of Nodes with a Given Degree (log)')
    plt.title('Degree Distribution for Graphs 1 and 2')
    plt.legend()
    plt.savefig('deg-dist.png')


 
    ###########################################################################


def main():
    Q1()
    Q2()
    Q3()
    Q4()


if __name__ == "__main__":
    main()
