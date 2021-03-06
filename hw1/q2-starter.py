from sets import Set
import snap
import matplotlib
matplotlib.use('agg')
import numpy as np
import matplotlib.pyplot as plt

def load_graph(name):
    '''
    Helper function to load graphs.
    Wse "epinions" for Epinions graph and "email" for Email graph.
    Check that the respective .txt files are in the same folder as this script;
    if not, change the paths below as required.
    '''
    if name == "epinions":
        G = snap.LoadEdgeList(snap.PNGraph, "soc-Epinions1.txt", 0, 1)
    elif name == 'email':
        G = snap.LoadEdgeList(snap.PNGraph, "email-EuAll.txt", 0, 1)   
    else: 
        raise ValueError("Invalid graph: please use 'email' or 'epinions'.")
    return G

def q2_1():
    '''
    You will have to run the inward and outward BFS trees for the 
    respective nodes and reason about whether they are in SCC, IN or OUT.
    You may find the SNAP function GetBfsTree() to be useful here.
    '''
    def in_out_bfs_for_node(graph, node):
        print 'Starting at node: ' , str(node)
        bfs_out = snap.GetBfsTree(graph, node, True, False)
        bfs_in = snap.GetBfsTree(graph, node, False, True)
        print 'Out size = ' , str(bfs_out.GetEdges()) , ', In size = ' , str(bfs_in.GetEdges())


    def in_out_bfs(graph, name):
        print 'In/Out BFS tree size for graph: ' , name
        if name == 'email':
            in_out_bfs_for_node(graph, 2018)
        else:
            in_out_bfs_for_node(graph, 224)
        
        
    ##########################################################################
    #TODO: Run outward and inward BFS trees from node 2018, compare sizes 
    #and comment on where node 2018 lies.
    G = load_graph("email")
    #Your code here:
    in_out_bfs(G, 'email')
    
    
    ##########################################################################
    
    ##########################################################################
    #TODO: Run outward and inward BFS trees from node 224, compare sizes 
    #and comment on where node 224 lies.
    G = load_graph("epinions")
    #Your code here:
    in_out_bfs(G, 'epinions')
    
    
    
    
    
    ##########################################################################

    print '2.1: Done!\n'


def q2_2():
    '''
    For each graph, get 100 random nodes and find the number of nodes in their
    inward and outward BFS trees starting from each node. Plot the cumulative
    number of nodes reached in the BFS runs, similar to the graph shown in 
    Broder et al. (see Figure in handout). You will need to have 4 figures,
    one each for the inward and outward BFS for each of email and epinions.
    
    Note: You may find the SNAP function GetRndNId() useful to get random
    node IDs (for initializing BFS).
    '''
    ##########################################################################
    #TODO: See above.
    #Your code here:
    
    def reachable_analysis(graph, name):
        reached_out_cum = []
        reached_in_cum = []
        for i in range(100):
            node = graph.GetRndNId()
            bfs_out = snap.GetBfsTree(graph, node, True, False)
            bfs_in = snap.GetBfsTree(graph, node, False, True)
            reached_out_cum.append(bfs_out.GetNodes())
            reached_in_cum.append(bfs_in.GetNodes())

        reached_out_cum = sorted(reached_out_cum)
        reached_in_cum = sorted(reached_in_cum)

        y = reached_out_cum
        if True:#name == 'epinions':
            plt.plot(y, label='Reachability using outlinks for %s' % name)
            plt.yscale('log')
        #else:
        #    plt.plot(y, label='Reachability using outlinks for %s' % name)
        plt.savefig('outlinks_%s.png' % name)
        plt.close()

        y = reached_in_cum
        if True:#name == 'email':
            plt.plot(y, label='Reachability using inlinks for %s' % name)
            plt.yscale('log')
        #else:
        #    plt.plot(y, label='Reachability using inlinks for %s' % name)
        plt.savefig('inlinks_%s.png' % name)
        plt.close()

    reachable_analysis(load_graph('email'), 'email')
    reachable_analysis(load_graph('epinions'), 'epinions')
    ##########################################################################
    print '2.2: Done!\n'

def q2_3():
    '''
    For each graph, determine the size of the following regions:
        DISCONNECTED
        IN
        OUT
        SCC
        TENDRILS + TUBES
        
    You can use SNAP functions GetMxWcc() and GetMxScc() to get the sizes of 
    the largest WCC and SCC on each graph. 
    '''
    ##########################################################################
    #TODO: See above.
    #Your code here:
    def per_graph(graph, name):
        mxWcc = snap.GetMxWcc(graph)
        mxScc = snap.GetMxScc(graph)
        print ''
        print 'Size analysis on {}'.format(name)
        print 'Disconnected size = {}'.format(graph.GetNodes() - mxWcc.GetNodes())
        print 'SCC size = {}'.format(mxScc.GetNodes())
        
        trials = 200
        avg_reached_out = 0
        avg_reached_in = 0
        for _ in range(trials):
            nodeId = mxScc.GetRndNId()
            avg_reached_out += snap.GetBfsTree(graph, nodeId, True, False).GetNodes()
            avg_reached_in += snap.GetBfsTree(graph, nodeId, False, True).GetNodes()

        scc_out = float(avg_reached_out) / trials
        scc_in = float(avg_reached_in) / trials

        out_sz = scc_out - mxScc.GetNodes()
        in_sz = scc_in - mxScc.GetNodes()
        print 'OUT size = {}'.format(out_sz)
        print 'IN size = {}'.format(in_sz)
        print 'Tendrils/Tubes size = {}'.format(mxWcc.GetNodes() - mxScc.GetNodes() - out_sz - in_sz)
    

    per_graph(load_graph('email'), 'email')
    per_graph(load_graph('epinions'), 'epinions')
    
    
    
    
    
    
    
    
    
    
    
    
    
    ##########################################################################
    print '2.3: Done!\n' 

def q2_4():
    '''
    For each graph, calculate the probability that a path exists between
    two nodes chosen uniformly from the overall graph.
    You can do this by choosing a large number of pairs of random nodes
    and calculating the fraction of these pairs which are connected.
    The following SNAP functions may be of help: GetRndNId(), GetShortPath()
    '''
    ##########################################################################
    #TODO: See above.
    #Your code here
    
    def calc_path_prob(graph):
        trials = 1000
        reachable_count = 0
        for _ in range(trials):
            node1 = graph.GetRndNId()
            node2 = graph.GetRndNId()
            shortestPath = snap.GetShortPath(graph, node1, node2, True)
            if shortestPath > 0:
                reachable_count += 1
        return float(reachable_count) / float(trials)
    
    print 'email reachable prob = %f' % calc_path_prob(load_graph('email'))
    print 'epinions reachable prob = %f' % calc_path_prob(load_graph('epinions'))
    ##########################################################################
    print '2.4: Done!\n'
    
if __name__ == "__main__":
    q2_1()
    q2_2()
    q2_3()
    q2_4()
    print "Done with Question 2!\n"
