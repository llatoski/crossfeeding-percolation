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

awk '$3<=0.{print $1,$2,$3}' frames/frame$1 > pos_Acells$1.dat
awk '$3>0.{print $1,$2,$3}' frames/frame$1 > pos_Ecells$1.dat
Nemit=$(cat pos_Ecells$1.dat | wc -l)
Nabsorb=$(cat pos_Acells$1.dat | wc -l)
python3 voronoi.py $1 > aux$1
gcc current_numerical_voronoi_v2.0.c -lm -DNEMIT="$Nemit" -DNABSORB="$Nabsorb" -DFRAME="$1" -o "$1".out
./"$1".out < aux$1 
rm -rf pos_Ecells$1.dat
rm -rf pos_Acells$1.dat
rm aux$1
rm "$1".out