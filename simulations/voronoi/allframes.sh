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

THRESHOLD=0.025
for j in $(seq $1 21 3100);
do
    bash framing.sh $j
    bash single_frame.sh $j
    bash printedges.sh $j $THRESHOLD > edges$j.dsf
    count=$(cat edges$j.dsf | wc -l)
    if [ "$count" == 0 ]; 
    then
        printf '%s 0\n' $j >> "$2"aux
        # python3 average.py $j >> $2
    else
        python3 networkstatistics.py $j >> "$2"aux 
    fi
    rm edges$j.dsf
done