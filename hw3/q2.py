import networkx as nx
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np

h = nx.read_edgelist(path='graph/karate.edgelist', delimiter=' ')

nx.draw(h, with_labels=True)
plt.savefig('karate.png')

vecs = {}
with open('emb/karate.emb') as f:
    lines = f.readlines()
    num_nodes, dim = lines[0].split(' ')
    num_nodes = int(num_nodes)
    dim = int(dim)
    for line in lines[1:]:
        sp = line.split(' ')
        node = int(sp[0])
        vec = [float(feature) for feature in sp[1:]]
        assert(len(vec) == dim)
        vecs[node] = np.array(vec)

ref_node = 34
ref = vecs[ref_node]
closest_nodes = []
closest_dists = []
for node in vecs:
    if node == ref_node:
        continue
    dist = np.linalg.norm(ref - vecs[node])
    if len(closest_nodes) < 5:
        closest_nodes.append(node)
        closest_dists.append(dist)
    else:
        i = np.argmax(closest_dists)
        if dist < closest_dists[i]:
            closest_nodes[i] = node
            closest_dists[i] = dist

print(closest_nodes)

