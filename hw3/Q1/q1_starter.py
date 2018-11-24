import snap
import numpy as np
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
import random

def get_neighbors(graph, nodeId):
    n = set([])
    for i in range(graph.GetNI(nodeId).GetDeg()):
        nid = graph.GetNI(nodeId).GetNbrNId(i)
        n.add(nid)
    return n


def load_graph(name):
    '''
    Helper function to load undirected graphs.
    Check that the respective .txt files are in the same folder as this script;
    if not, change the paths below as required.
    '''
    if name == "normal":
        G = snap.LoadEdgeList(snap.PUNGraph, "polblogs.txt", 0, 1)
    elif name == 'rewired':
        G = snap.LoadEdgeList(snap.PUNGraph, "polblogs-rewired.txt", 0, 1)
    elif name == 'sample':
        G = snap.LoadEdgeList(snap.PUNGraph, "q1-sample.txt", 0, 1)
    else:
        raise ValueError("Invalid graph: please use 'normal', 'rewired' or 'sample'.")
    return G

def get_adjacency_matrix(Graph):
    '''
    This function might be useful for you to build the adjacency matrix of a
    given graph and return it as a numpy array
    '''
    ##########################################################################
    A = np.zeros((Graph.GetNodes(), Graph.GetNodes()))
    nid_map = {}
    i = 0
    for node in Graph.Nodes():
        nid_map[node.GetId()] = i
        i += 1
    for edge in Graph.Edges():
        i = nid_map[edge.GetSrcNId()]
        j = nid_map[edge.GetDstNId()]
        A[i][j] = 1
        A[j][i] = 1
    
    return A, nid_map
    ##########################################################################

def get_sparse_degree_matrix(Graph, nid_map):
    '''
    This function might be useful for you to build the degree matrix of a
    given graph and return it as a numpy array
    '''
    ##########################################################################
    #TODO: Your code here
    D = np.zeros((Graph.GetNodes(), Graph.GetNodes()))
    for node in nid_map:
        i = nid_map[node]
        D[i][i] = Graph.GetNI(node).GetDeg()
    return D
    ##########################################################################

def normalized_cut_minimization(Graph, other_A=None):
    '''
    Implement the normalized cut minimizaton algorithm we derived in the last
    homework here
    '''
    A, nid_map = get_adjacency_matrix(Graph)
    if other_A is not None:
        assert(np.array_equal(A, other_A))
    D = get_sparse_degree_matrix(Graph, nid_map)
    ##########################################################################
    #TODO: Your code here
    d_inv_root = np.linalg.inv(np.sqrt(D))
    _, v = np.linalg.eigh(np.dot(np.dot(d_inv_root, D - A), d_inv_root))
    x = np.dot(d_inv_root, v[:, 1])
    for i in range(len(x)):
        if x[i] <= 0:
            x[i] = -1
        else:
            x[i] = 1
    return get_assignment_sets(Graph, x, nid_map)
    ##########################################################################

def invert_map(mapping):
    return {mapping[key]: key for key in mapping.keys()}

def get_assignment_sets(Graph, x, nid_map):
    c1 = set([])
    c2 = set([])
    for node in Graph.Nodes():
        if x[nid_map[node.GetId()]] == 1:
            c1.add(node.GetId())
        else:
            c2.add(node.GetId())
    return c1, c2

def modularity(Graph, c1, c2):
    '''
    This function might be useful to compute the modularity of a given cut
    defined by two sets c1 and c2. We would normally require sets c1 and c2
    to be disjoint and to include all nodes in Graph
    '''
    ##########################################################################
    #TODO: Your code here
    assert(len(c1 & c2) == 0)
    assert(len(c1 | c2) == Graph.GetNodes())
    print len(c1)
    print len(c2)
    modularity = 0.
    for inode in Graph.Nodes():
        for jnode in Graph.Nodes():
            i = inode.GetId()
            j = jnode.GetId()
            if (i in c1 and j in c1) or (i in c2 and j in c2):
                aij = 1 if j in get_neighbors(Graph, i) else 0
                modularity += aij - (inode.GetDeg() * jnode.GetDeg()) / (2. * Graph.GetEdges())

    modularity = modularity / (2. * Graph.GetEdges())        
    print 'Modularity = {}'.format(modularity)
    return modularity
    ##########################################################################

def q1_1():
    '''
    Main q1_1 workflow. All of the work can be done in gen_config_model_rewire
    but you may modify this function as needed.
    '''
    ##########################################################################
    #TODO: Your code here
    graph = load_graph('sample')
    c1, c2 = normalized_cut_minimization(graph)
    print 'Sample'
    modularity(graph, c1, c2)
    graph = load_graph('normal')
    c1, c2 = normalized_cut_minimization(graph)
    print 'Polblogs'
    modularity(graph, c1, c2)
    graph = load_graph('rewired')
    c1, c2 = normalized_cut_minimization(graph)
    print 'Polblogs rewired'
    modularity(graph, c1, c2)
    ##########################################################################

def SSBM(n,pi,a,b, name):
    '''

    '''
    ##########################################################################
    #TODO: Your code here

    def get_group_index(nindex, cum_npi):
        index = 0
        for i in range(len(cum_npi)-1):
            if nindex < cum_npi[i+1]:
                break
            index += 1
        return index

    g = snap.TUNGraph.New()
    for i in range(n):
        g.AddNode(i)
    pi = np.array(pi)
    assert(np.sum(pi) == 1)
    npi = np.zeros(len(pi))
    for i in range(n):
        assigment = np.random.choice(list(range(len(npi))), size=1, p=pi)
        npi[assigment] += 1
    assert(sum(npi) == n)
    cum_npi = np.zeros(len(npi))
    for i in range(1, len(npi)):
        cum_npi[i] = cum_npi[i-1] + npi[i-1]
    A = np.zeros((n, n))
    for r in range(n):
        rgroup = get_group_index(r, cum_npi)
        for c in range(n):
            if r == c:
                continue
            cgroup = get_group_index(c, cum_npi)
            if rgroup == cgroup:
                rand = random.random()
                A[r][c] = 1 if rand < a else 0
                A[c][r] = 1 if rand < a else 0
            else:
                rand = random.random()
                A[r][c] = 1 if rand < b else 0
                A[c][r] = 1 if rand < b else 0

    # get clusters
    clusters = [set([])]
    c = 0
    for nodei in range(n):
        if c == len(cum_npi) - 1 or nodei < cum_npi[c + 1]:
            pass
        else:
            clusters.append(set([]))
            c += 1
        clusters[-1].add(nodei)
        
    if name != '':
        plt.imsave('ssbm{}.png'.format(name), A, cmap='binary')
    return A, clusters
    
    ##########################################################################

def q1_2():
    '''
    You can probably just implement everything required for question 1.2 here,
    but feel free to create additional helper functions, should you need them
    '''
    ##########################################################################
    #TODO: Your code here
    n = 10
    pi = (0.5,0.5)
    tin = 0.9
    tout = 0.3
    A, clusters = SSBM(n, pi, tin, tout, '0')
    print A
    print clusters
    n = 100
    pi = (0.1,0.2,0.3,0.4)
    tin = 0.6
    tout = 0.15
    SSBM(n, pi, tin, tout, '1')
    n = 1000
    pi = (0.25, 0.25,0.25,0.25)
    tin = 0.2
    tout = 0.6
    SSBM(n, pi, tin, tout, '2')
    n = 1000
    pi = (0.25,0.35,0.1,0.3)
    tin = 0.7
    tout = 0.5
    SSBM(n, pi, tin, tout, '3')

    ##########################################################################

def get_accuracy(c1, c2, c1_hat, c2_hat):
    '''
    Compute the accuracy of an assignment here!
    '''
    ##########################################################################
    #TODO: Your code here


    def compute_acc(a, b, ahat, bhat):
        return 2. * ((len(a & ahat) + len(b & bhat)) / (float(len(a) + len(b))) - 0.5)

    return max(compute_acc(c1, c2, c1_hat, c2_hat), compute_acc(c1, c2, c2_hat, c1_hat))
    ##########################################################################

def generate_graph_from_matrix(A):
    g = snap.TUNGraph.New()
    for i in range(len(A)):
        g.AddNode(i)
    for i in range(len(A)):
        for j in range(len(A)):
            if A[i][j]:
                g.AddEdge(i, j)
    return g

def generate_acc_plots(vary_theta_in):
    n = 1000
    pi = (0.5,0.5)
    acc = []
    for theta in np.linspace(0.0, 1.0, 50):
        if vary_theta_in:
            tout = 0.3
            A, clusters = SSBM(n, pi, theta, tout, '')
        else:
            tin = 0.7
            A, clusters = SSBM(n, pi, tin, theta, '')
        assert(len(clusters) == 2)
        c1, c2 = clusters
        G = generate_graph_from_matrix(A)
        c1_hat, c2_hat = normalized_cut_minimization(G, other_A=A)
        pred_acc = get_accuracy(c1, c2, c1_hat, c2_hat)
        print(pred_acc)
        acc.append(pred_acc)
    
    x = np.linspace(0.0, 1.0, 50)
    plt.plot(x, acc)
    plt.ylabel('Accuracy')
    xlabel = 'Theta_in' if vary_theta_in else 'Theta_out'
    plt.xlabel(xlabel)
    a = pi[0] * 1000.
    b = pi[1] * 1000.
    exact_recovery = (a + b) / 2 - np.sqrt(a * b)
    weak_recovery = ((a - b) ** 2) / (2 * a * b)
    plt.axvline(x=exact_recovery, color='r', linestyle='--')
    plt.axvline(x=weak_recovery, color='g', linestyle='--')
    plt.savefig('acc_plot_{}.png'.format(xlabel))
    plt.clf()

def q1_3():
    '''
    You can probably just implement everything required for question 1.3 here,
    but feel free to create additional helper functions, should you need them
    '''
    ##########################################################################
    #TODO: Your code here

    n = 1000
    pi = (0.5,0.5)
    tin = 0.55
    tout = 0.5
    A, clusters = SSBM(n, pi, tin, tout, '_p3_graph')
    assert(len(clusters) == 2)
    c1, c2 = clusters
    G = generate_graph_from_matrix(A)
    c1_hat, c2_hat = normalized_cut_minimization(G, other_A=A)
    print c2_hat
    print 'First graph accuracy = {}'.format(get_accuracy(c1, c2, c1_hat, c2_hat))
    
    generate_acc_plots(True)
    generate_acc_plots(False)
    
    ##########################################################################

if __name__ == "__main__":
    # Questions
    q1_2()
    q1_3()
    q1_1()
    print "Done with Question 1!\n"
