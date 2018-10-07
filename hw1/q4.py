import snap
from sets import Set
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
            print 'node {}'.format(i)
            for j in dgraph.Nodes():
                if i.GetId() == j.GetId():
                    continue
                if graph.IsEdge(i.GetId(), j.GetId()):
                    continue
                for y in ggraph.Nodes():
                    if graph.IsEdge(i.GetId(), y.GetId()) and graph.IsEdge(y.GetId(), j.GetId()):
                        dgraph.AddEdge(i.GetId(), j.GetId())
        fout = snap.TFOut('hdn.graph')
        dgraph.Save(fout)
        fout.Flush()

    print ''
    print 'HDN graph'
    print 'nodes = {}'.format(dgraph.GetNodes())
    print 'edges = {}'.format(dgraph.GetEdges())
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
    return node.GetId() > 200E6

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


def main():
    graph = load_graph()
    print_graph_stats(graph)
    make_plots(graph)    
    build_hdn(graph)

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
