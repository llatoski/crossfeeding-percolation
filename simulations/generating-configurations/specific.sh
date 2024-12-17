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

betao="$1"
betag="$2"
output=$(printf 'betao%s_betag%s.dat' $betao $betag)
rm $output
g++ teste.cpp -D BETA=$betao -D BETA2=$betag -D TIMEWAIT=0 -D SAMPLETIME=100 -o exec"$betao"_"$betag".out
./exec"$betao"_"$betag".out >> $output
rm exec"$betao"_"$betag".out