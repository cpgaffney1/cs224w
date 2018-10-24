import snap
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

def get_neighbors(graph, nodeId):
    n = set([])
    for i in range(graph.GetNI(nodeId).GetDeg()):
        nid = graph.GetNI(nodeId).GetNbrNId(i)
        n.add(nid)
    return n
 

def load_graph(q):
    '''
    Helper function to load graphs.
    Wse "epinions" for Epinions graph and "email" for Email graph.
    Check that the respective .txt files are in the same folder as this script;
    if not, change the paths below as required.
    '''    
    if q == 1:
        G = snap.TUNGraph.Load(snap.TFIn('Q1/hw2-q1.graph'))
    elif q == 3:
        G = snap.TNGraph.Load(snap.TFIn('Q3/USpowergrid_n4941.txt'))
    return G

def sim(x, y):
    if np.linalg.norm(x) == 0 or np.linalg.norm(y) == 0:
        return 0
    return np.dot(x, y) / (np.linalg.norm(x) * np.linalg.norm(y))

def compute_topk_feature_sim(features, node=9, k=5):
    all_nodes = []
    all_sims = []
    most_similar_nodes = []
    similarity = []
    for key in features.keys():
        if key == node:
            continue
        score = sim(features[node], features[key])
        all_nodes.append(key)
        all_sims.append(score)
        if len(similarity) < k:
            similarity.append(score)
            most_similar_nodes.append(key)
        else:
            if score > min(similarity):
                ind = np.argmin(similarity)
                similarity[ind] = score
                most_similar_nodes[ind] = key
    similarity = np.array(similarity)
    most_similar_nodes = np.array(most_similar_nodes)
    sorted_ind = np.flip(np.argsort(similarity))
    similarity = similarity[sorted_ind]
    most_similar_nodes = most_similar_nodes[sorted_ind]
    print 'Finding most similar nodes to {}'.format(node) 
    print most_similar_nodes
    print similarity
    return all_nodes, all_sims


def collect_features(graph):
    print '\nCollecting features'
    features = {}
    for node in graph.Nodes():
        deg = node.GetDeg()
        in_edges = set([])
        out_edges = set([])
        neighbors = get_neighbors(graph, node.GetId())
        for i in neighbors:
            for j in get_neighbors(graph, i):
                if j in neighbors and ((i, j) not in in_edges and (j, i) not in in_edges):
                    in_edges.add((i, j))
                if j not in neighbors and ((i, j) not in out_edges and (j, i) not in out_edges):
                    out_edges.add((i, j))
        features[node.GetId()] = np.array([deg, len(in_edges), len(out_edges)])
    print 'Features of 9'
    print features[9]
    
    compute_topk_feature_sim(features)

    return features

def plot_hist(x, bins):
    plt.hist(x, bins=bins)
    plt.savefig('q1_hist.png')

def collect_recursive_features(graph, features, k=2):
    for i in range(k):
        n_starting_features = len(features[graph.GetRndNId()])
        for node in graph.Nodes():
            s = np.zeros(n_starting_features)
            nbrs = get_neighbors(graph, node.GetId())
            for n in nbrs:
                s += features[n][:n_starting_features]
            new_features = np.zeros(3 * len(s))
            new_features[:len(s)] = features[node.GetId()]
            if len(nbrs) > 0:
                new_features[len(s):2*len(s)] = s
                new_features[2*len(s):] = s / len(nbrs)
            features[node.GetId()] = new_features
    
    all_nodes, all_sims = compute_topk_feature_sim(features)
    buckets = np.linspace(np.min(all_sims), np.max(all_sims), num=20)
    print buckets
    bucket_ind = np.digitize(all_sims, buckets)
    plot_hist(all_sims, buckets)
    

    ind = np.argwhere(bucket_ind == 1)[0][0]
    print ''
    print ''
    print all_nodes[ind]
    print all_sims[ind]
    print features[all_nodes[ind]]
    print len(get_neighbors(graph, all_nodes[ind]))

    ind = np.argwhere(bucket_ind == 19)[1][0]
    print ''
    print ''
    print all_nodes[ind]
    print all_sims[ind]
    print features[all_nodes[ind]]
    print len(get_neighbors(graph, all_nodes[ind])) 
    vis_subgraph(graph, all_nodes[ind])

    ind = np.argwhere(bucket_ind == 14)[0][0]
    print ''
    print ''
    print all_nodes[ind]
    print all_sims[ind]
    print features[all_nodes[ind]]
    print len(get_neighbors(graph, all_nodes[ind]))
    vis_subgraph(graph, all_nodes[ind])



def vis_subgraph(graph, node):
    nodeSet = set([])
    nodeSet.add(node)
    for n in get_neighbors(graph, node):
        nodeSet.add(n)
        for n2 in get_neighbors(graph, n):
            nodeSet.add(n2)
            
    nodeV = snap.TIntV(len(nodeSet))
    i = 0
    for n in nodeSet:
        nodeV[i] = n
        i += 1
    sg = snap.ConvertSubGraph(snap.PUNGraph, graph, nodeV, False)
    snap.DrawGViz(sg, snap.gvlNeato, 'subgraph{}.png'.format(node), "Subgraph around node {}".format(node), True) 

def main():
    #q1
    graph = load_graph(1)
    features = collect_features(graph)
    collect_recursive_features(graph, features)

main()
