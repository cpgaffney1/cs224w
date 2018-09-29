import matplotlib
matplotlib.use("agg")
import snap
import matplotlib.pyplot as plt
import numpy as np

def stackoverflow():
    g = snap.LoadEdgeList(snap.PNGraph, "stackoverflow-Java.txt", 0, 1)
    components = snap.TCnComV()
    snap.GetWccs(g, components)
    print "Num connected comp = ", components.Len()
    mxwcc = snap.GetMxWcc(g)
    print "Num edges in largest = ", mxwcc.GetEdges()
    print "Num nodes in largest = ", mxwcc.GetNodes()
    rank = snap.TIntFltH()
    snap.GetPageRank(g, rank)
    rank.SortByDat(False)
    count = 0
    for node in rank:
        if count >= 3:
            break
        count += 1
        print "largest page rank score nodes = ", node, " (score = ", rank[node]

    hubs = snap.TIntFltH()
    auths = snap.TIntFltH()
    snap.GetHits(g, hubs, auths)
    
    hubs.SortByDat(False)
    count = 0
    for node in hubs:
        if count >= 3:
            break
        count += 1
        print "largest hub score nodes = ", node, " (score = ", hubs[node]

    auths.SortByDat(False)
    count = 0
    for node in auths:
        if count >= 3:
            break
        count += 1
        print "largest auth score nodes = ", node, " (score = ", auths[node]




def make_plots():
    gwiki = snap.LoadEdgeList(snap.PNGraph, "wiki-Vote.txt", 0, 1)
    out_deg_counts = {}
    for node in gwiki.Nodes():
        if node.GetOutDeg() not in out_deg_counts.keys():
            out_deg_counts[node.GetOutDeg()] = 0
        out_deg_counts[node.GetOutDeg()] += 1
    x = []
    y = []
    for key in sorted(out_deg_counts.keys()):
        if key == 0 or out_deg_counts[key] == 0:
            continue
        x.append(key)
        y.append(out_deg_counts[key])

    x = np.log10(x)
    y = np.log10(y)
    plt.plot(x, y, linestyle="", marker="o")
    plt.savefig('out-degree-hist.png')
    
    coef = np.polyfit(x, y, 1)
    print coef
    print ""



def wiki_analysis(gwiki=None):
    if gwiki is None:
        gwiki = snap.LoadEdgeList(snap.PNGraph, "wiki-Vote.txt", 0, 1)
    print("n nodes = ", gwiki.GetNodes()) 
    n_self_edges = 0
    n_zero_out_deg = 0
    n_zero_in_deg = 0
    n_large_out_deg = 0
    n_small_in_deg = 0
    for node in gwiki.Nodes():
        if gwiki.IsEdge(node.GetId(), node.GetId()):
            n_self_edges += 1
        if node.GetOutDeg() == 0:
            n_zero_out_deg += 1
        if node.GetInDeg() == 0:
            n_zero_in_deg += 1
        if node.GetOutDeg() > 10:
            n_large_out_deg += 1
        if node.GetInDeg() < 10:
            n_small_in_deg += 1
    print("n self edges = ", n_self_edges)
    print("n 0 out = ", n_zero_out_deg)
    print("n 0 in = ", n_zero_in_deg)
    print("n large out = ", n_large_out_deg)
    print("n small in = ", n_small_in_deg)
    n_directed = 0
    n_undirected = 0
    n_reciprocated = 0
    for edge in gwiki.Edges():
        if edge.GetSrcNId() != edge.GetDstNId():
            n_directed += 1
            if gwiki.IsEdge(edge.GetDstNId(), edge.GetSrcNId()):
                n_reciprocated += 0.5
                n_undirected += 0.5
            else: 
                n_undirected += 1
    print("n directed = ", n_directed)
    print("n undirected = ", n_undirected)
    print("n recip = ", n_reciprocated)
    print ""

stackoverflow()

make_plots()

# TODO fix undirected count
gtest = snap.TNGraph.New()
gtest.AddNode(1)
gtest.AddNode(2)
gtest.AddNode(3)
gtest.AddEdge(1,2)
gtest.AddEdge(2,1)
gtest.AddEdge(1,3)
gtest.AddEdge(1,1)
wiki_analysis(gwiki=gtest)

wiki_analysis()
