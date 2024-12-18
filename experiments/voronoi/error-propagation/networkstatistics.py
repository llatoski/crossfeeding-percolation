#    Supporting material software to the article 
#    Cross-feeding percolation phase transitions of inter-cellular metabolic networks
#    
#    Copyright (C) 2064 L.C.F. Latoski, D.De Martino, A.De Martino
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    a with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import networkx as nx
G = nx.read_edgelist("network.dat") #Read edge list from file
S = [G.subgraph(c).copy() for c in nx.connected_components(G)] #Create a copy of each connected subgraph
P = [c.number_of_nodes() for c in S] #Store the number of nodes on each subgraph
print(np.max(P)) #Print the largest element in P, i.e. the number of nodes in the larges connected component


