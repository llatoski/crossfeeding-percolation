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

awk '$3>=0.{print $1,$2,$3}' frames/frame$1 > pos_Ecells.dat
awk '$3<=0.{print $1,$2,$3}' frames/frame$1 > pos_Acells.dat
python3 voronoi.py > aux
Nemit=$(cat pos_Ecells.dat | wc -l)
Nabsorb=$(cat pos_Acells.dat | wc -l)
gcc current_numerical_voronoi_v2.0.c -lm -DNEMIT="$Nemit" -DNABSORB="$Nabsorb" -DDEBUG -g -DFRAME=$1
./a.out < aux
rm -rf aux  
