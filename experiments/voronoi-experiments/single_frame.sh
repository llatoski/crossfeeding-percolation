awk '$3>=0.{print $1,$2,$3}' frames/frame$1 > pos_Ecells.dat
awk '$3<=0.{print $1,$2,$3}' frames/frame$1 > pos_Acells.dat
python3 voronoi.py > aux
Nemit=$(cat pos_Ecells.dat | wc -l)
Nabsorb=$(cat pos_Acells.dat | wc -l)
gcc current_numerical_voronoi_v2.0.c -lm -DNEMIT="$Nemit" -DNABSORB="$Nabsorb" -DDEBUG -g -DFRAME=$1
./a.out < aux
rm -rf aux  
