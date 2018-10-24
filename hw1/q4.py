import snap
from sets import Set
import os
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

def build_hdn(graph):
    if os.path.exists('hdn.graph'):
        fin = snap.TFIn('hdn.graph')
        dgraph = snap.TUNGraph.Load(fin)
    else:
        dgraph = snap.TUNGraph.New()
        ggraph = snap.TUNGraph.New()
        for node in graph.Nodes():
            if is_gene(node):
                ggraph.AddNode(node.GetId())
            else:
               dgraph.AddNode(node.GetId())
        for i in dgraph.Nodes():
            print 'node {}'.format(i.GetId())
            j = dgraph.GetNI(i.GetId())
            j.Next()
            while j < dgraph.EndNI():
                if i.GetId() == j.GetId():
                    continue
                if dgraph.IsEdge(i.GetId(), j.GetId()):
                    continue
                i_neighbors = get_neighbors(graph, i.GetId())
                j_neighbors = get_neighbors(graph, j.GetId())
                if len(i_neighbors.intersection(j_neighbors)) > 0:
                    dgraph.AddEdge(i.GetId(), j.GetId())
                j.Next()

        fout = snap.TFOut('hdn.graph')
        dgraph.Save(fout)
        fout.Flush()

    print ''
    print 'HDN graph'
    print 'nodes = {}'.format(dgraph.GetNodes())
    print 'edges = {}'.format(dgraph.GetEdges())
    print 'density = {}'.format(2.0 * dgraph.GetEdges() / (dgraph.GetNodes() * (dgraph.GetNodes() - 1.0)))
    print 'clust coef = {}'.format(snap.GetClustCf(dgraph, 1000))
    return dgraph

def load_graph():
    '''
    Helper function to load graphs.
    Wse "epinions" for Epinions graph and "email" for Email graph.
    Check that the respective .txt files are in the same folder as this script;
    if not, change the paths below as required.
    '''    
    G = snap.LoadEdgeList(snap.PUNGraph, 'all_gene_disease_associations.csv', 0, 2, ',')
    assert not G.Empty()
    assert G.IsNode(200001418)
    assert G.IsNode(10)
    return G

def is_gene(node):
    return node.GetId() < 200E6

def make_plots(graph):
    gene_deg_counts = {}
    disease_deg_counts = {}
    for node in graph.Nodes():
        if is_gene(node):
            if node.GetDeg() not in gene_deg_counts.keys():
                gene_deg_counts[node.GetDeg()] = 0
            gene_deg_counts[node.GetDeg()] += 1
        else:
            if node.GetDeg() not in disease_deg_counts.keys():
                disease_deg_counts[node.GetDeg()] = 0
            disease_deg_counts[node.GetDeg()] += 1
        
    deg_counts_plot(disease_deg_counts, 'blue')
    deg_counts_plot(gene_deg_counts, 'green')
    plt.savefig('deg_counts_plot.png')


def deg_counts_plot(deg_counts, color):
    x = []
    y = []
    for key in sorted(deg_counts.keys()):
        if key == 0 or deg_counts[key] == 0:
            continue
        x.append(key)
        y.append(deg_counts[key])

    x = np.log10(x)
    y = np.log10(y)
    plt.plot(x, y, color=color, linestyle="", marker="o")

def get_neighbors(graph, nodeId, othergraph=None):
    n = set([])
    for i in range(graph.GetNI(nodeId).GetDeg()):
        nid = graph.GetNI(nodeId).GetNbrNId(i)
        if othergraph is not None:
            if not othergraph.IsNode(nid):
                continue
        n.add(nid)
    return n
   
def max_clique(graph):
    mxcl = 0
    for node in graph.Nodes():
        if is_gene(node):
            sz = node.GetDeg()
            if sz > mxcl:
                mxcl = sz
    print ''
    print 'HDN max clique size = {}'.format(mxcl)


def disease_sim(graph):
    crohnid = 200000000 + 10346
    leukid = 200000000 + 23418

    def top_five_sim(nid, metric):
        if nid == crohnid:
            name = 'crohn'
        else:
            name = 'leuk'
        print ''
        print 'Scores for {} using {}'.format(name, metric)
        top_five = []
        top_five_scores = []
        for node in graph.Nodes():
            if is_gene(node):
                continue
            if nid == crohnid and node.GetId() == crohnid:
                continue
            if nid == leukid and node.GetId() == leukid:
                continue
            nid_nb = get_neighbors(graph, nid)
            nb = get_neighbors(graph, node.GetId())
            if metric == 'CN':
                score = len(nid_nb.intersection(nb))
            else:
                score = len(nid_nb.intersection(nb)) / (1.0 * len(nid_nb.union(nb)))
            if len(top_five) < 5:
                top_five.append(node.GetId())
                top_five_scores.append(score)
            else:
                min_ind = np.argmin(top_five_scores)
                if score > top_five_scores[min_ind]:
                    top_five[min_ind] = node.GetId()
                    top_five_scores[min_ind] = score
        print top_five
        print top_five_scores

    assert(graph.IsNode(crohnid))
    assert(graph.IsNode(leukid))
    top_five_sim(crohnid, 'CN')
    top_five_sim(crohnid, 'JA')
    top_five_sim(leukid, 'CN')
    top_five_sim(leukid, 'JA')

def contraction(graph, hdn):
    
    supernodes = set([])

    def contract_clique(node):
        clique = get_neighbors(graph, node.GetId(), othergraph=hdn)
        supernodeId = clique.pop()
        assert(hdn.IsNode(supernodeId))
        supernodes.add(supernodeId)
        for nodeId in clique:
            if nodeId in supernodes:
                continue
            for nbrId in get_neighbors(hdn, nodeId):
                if not hdn.IsNode(nbrId):
                    continue
                assert(hdn.IsNode(supernodeId))
                assert(hdn.IsNode(nbrId))
                hdn.AddEdge(supernodeId, nbrId)
            assert(not supernodeId == nodeId)
            assert(nodeId not in supernodes)
            hdn.DelNode(nodeId)
            #graph.DelNode(nodeId)

    node = graph.BegNI()
    while node < graph.EndNI():
        if is_gene(node) and node.GetDeg() > 250:
            contract_clique(node)
        node.Next()
        
    contracted = hdn
    clust_cf = snap.GetClustCf(contracted, int(1000))
    print ''
    print 'Contracted graph stats'
    print 'clust cf = {}'.format(clust_cf)
    print 'density = {}'.format(2.0 * contracted.GetEdges() / (contracted.GetNodes() * (contracted.GetNodes() - 1.0)))
    print 'nodes = {}'.format(contracted.GetNodes())
    print 'edges = {}'.format(contracted.GetEdges())

def main():
    graph = load_graph()
    print_graph_stats(graph)
    make_plots(graph)    
    hdn = build_hdn(graph)
    contraction(graph, hdn)
    max_clique(graph)
    disease_sim(graph)
   
def print_graph_stats(graph):
    print 'nodes = {}'.format(graph.GetNodes())
    
    disease_set = set([])
    gene_set = set([])
    for node in graph.Nodes():
        if not is_gene(node):
            disease_set.add(node.GetId())
        else:
            gene_set.add(node.GetId())
    print 'num diseases = {}'.format(len(disease_set))
    print 'num genes = {}'.format(len(gene_set))
    print 'edges = {}'.format(graph.GetEdges())

main()
