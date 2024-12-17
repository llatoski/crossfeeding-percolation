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
import math
from scipy.spatial import Voronoi,voronoi_plot_2d

def dist(p1,p2):
  distance = math.sqrt( pow(p1[0]-p2[0],2) + pow(p1[1]-p2[1],2) )
  return distance
A = np.loadtxt("pos_Acells.dat",usecols=(0,1))
A0 = [i[0] for i in A]
A1 = [i[1] for i in A]
Fa = np.loadtxt("pos_Acells.dat",usecols=(2))
E = np.loadtxt("pos_Ecells.dat",usecols=(0,1))
Fe = np.loadtxt("pos_Ecells.dat",usecols=(2))
vor = Voronoi(A,incremental=True)

#Adding adjacent matrices
def copy_system_2d():
  Aux = [ [A0[i]+500,A1[i]] for i in range(len(A)) ]
  vor.add_points(Aux,restart=True)
  Aux = [ [A0[i],A1[i]-500] for i in range(len(A)) ]
  vor.add_points(Aux,restart=True)
  Aux = [ [A0[i],A1[i]+500] for i in range(len(A)) ]
  vor.add_points(Aux,restart=True)
  Aux = [ [A0[i]-500,A1[i]] for i in range(len(A)) ]
  vor.add_points(Aux,restart=True)
  Aux = [ [A0[i]-500,A1[i]+500] for i in range(len(A)) ]
  vor.add_points(Aux,restart=True)
  Aux = [ [A0[i]-500,A1[i]-500] for i in range(len(A)) ]
  vor.add_points(Aux,restart=True)
  Aux = [ [A0[i]+500,A1[i]+500] for i in range(len(A)) ]
  vor.add_points(Aux,restart=True)
  Aux = [ [A0[i]+500,A1[i]-500] for i in range(len(A)) ]
  vor.add_points(Aux,restart=True)
  return

def periodic_boundaries():
  for j in range(8):
    for i in range(len(A)):
      k=i+len(A)*(j+1)
      regions_copy[k] = i
  return

copy_system_2d()
regions_copy = np.arange(len(vor.points))
periodic_boundaries()

La=len(A)
Le=len(E)
box=np.zeros(shape=(Le))
G=np.zeros(shape=(La,len(vor.point_region)))
size = np.zeros(shape=(La,len(vor.point_region)))
norm = np.zeros(shape=La)
pe = np.zeros(shape=La)
neigh = np.zeros(shape=La)
dt = 1

def print_emiters():
  for n in range(Le):
    closer = -1
    distance = math.inf
    for m in range(La):
      if(dist(E[n],A[m])<distance):
        distance=dist(E[n],A[m])
        closer = m
    print(n,closer,( 1 - math.exp(-Fe[n]*dt) ) )
  return

def prob_matrix(): 
  for n in range(La):
    for m in range(len(vor.points)):
      pos = set(vor.regions[vor.point_region[n]]).intersection(vor.regions[vor.point_region[m]])
      G[n][m] = len(pos)
      if(vor.point_region[n]==vor.point_region[m]):
        G[n][m]=0 
      else:
        if(len(pos)>1): 
          size[n][m] = dist(vor.vertices[list(pos)[0]],vor.vertices[list(pos)[1]])
          norm[n]+=size[n][m]
          neigh[n]+=1
  return

def print_absorbers():  
  for n in range(La):
    pe[n] = norm[n]/(norm[n] + 2*math.pi*10)
    print(n,int(neigh[n]),pe[n],Fa[n])
  return

def neighborhood_open():
  for n in range(La):
    # a=0
    for m in range(len(vor.points)):
      if(G[n][m]>1): 
        if vor.point_region[m] not in vor.point_region[La:]: 
          print(n,m,size[n][m]/norm[n])
        else:
          print(n,-1,size[n][m]/norm[n])
        # a+=size[n][m]/norm[n]
    # print(a)
  return
      


def neighborhood_periodic():
  for n in range(La):
    for m in range(len(vor.points)):
      if(G[n][m]>1): 
          print(n,regions_copy[m],size[n][m]/norm[n])
    
  return

print_emiters()
prob_matrix()
print_absorbers()
neighborhood_open()
# neighborhood_periodic()