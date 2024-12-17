import numpy as np
import networkx as nx
G = nx.read_edgelist("network.dat")
S = [G.subgraph(c).copy() for c in nx.connected_components(G)]
P = [c.number_of_nodes() for c in S]
print(np.max(P))
#for c in S:
#    print(c.size())

