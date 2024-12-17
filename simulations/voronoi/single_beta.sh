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

STARTO=-4
ENDO=6
DELTA=$(bc <<< $STARTO-$ENDO)
betao=$(python3 -c "print($STARTO -"$1"*"$DELTA"/20.0)")
input=$(printf 'betao%s.dat' $betao)
# betao=1.0
# betag=-0.01
# input=$(printf 'betao%s_betag%s.dat' $betao $betag)
cat $input > cells.dat
rm frames/*
rm networks/*
# output=$(printf 'measures/clusters_betao%s_betag%s.dsf' $betao $betag)
output=$(printf 'measures/%s.dsf' $1)
rm "$output"aux
for i in $(seq 1 1 21);
do
    nohup bash allframes.sh $i $output &
done
