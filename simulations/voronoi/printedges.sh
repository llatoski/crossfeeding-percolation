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

cat networks/frame"$1"_Matrix.dat > aux"$1"
Nemit=$(cat aux"$1" | wc -l)
Nabsorb=$(awk '{print NF}' aux"$1" | sort -nu | tail -n 1)
cat frames/frame"$1" >> aux"$1"
gcc networkstatistics.c -lm -g -DEDGES -DNEMIT="$Nemit" -DNABSORB="$Nabsorb" -DTHRESHOLD="$2" -o "$1".out
./"$1".out < aux"$1"
rm aux"$1"
rm "$1".out