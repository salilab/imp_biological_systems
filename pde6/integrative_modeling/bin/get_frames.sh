#!/bin/bash

mkdir -p best_pdb
rm -f best_pdb/*.pdb

bin/process_output.py -f output/models_stats.out --soft -s Nframe -s Total_Score -s xlms_Score -s em_all_0 -s Domain_Linking_0 -s Excluded_Volume_0 -s Excluded_Volume_Symm_0 -s Template_Restraint_Inter -s Template_Restraint_Intra | sort -r -n -k 3 | tail -101 > best_pdb/scores.dat
frames=`cat best_pdb/scores.dat | grep -v '#' | awk '{print $2}'`

for i in $frames ;
do
echo $i
awk 'BEGIN{nf=0};($1=="MODEL"){nf=nf+1};($1!="MODEL"){if(nf=='"$i"'){print $0}}' output/trajectory.0.pdb > best_pdb/$i.A.pdb
awk 'BEGIN{nf=0};($1=="MODEL"){nf=nf+1};($1!="MODEL"){if(nf=='"$i"'){print $0}}' output/trajectory.0.symmetry.0.pdb | sed 's/ A / B /g' > best_pdb/$i.B.pdb
cat best_pdb/$i.A.pdb best_pdb/$i.B.pdb | grep -v "ENDMDL" > best_pdb/$i.pdb
done

mkdir -p clustering
rm -f best_pdb/*.?.pdb
rm -f clustering/clus*.pdb
rm -f clustering/clus*.score.dat
rm -f clustering/cluster.dat
rm -f clustering/coor.xyz
rm -f clustering/scores.dat

onepdb=`ls best_pdb/*.pdb | awk '(NR==1){print }'`

nconf=`ls best_pdb/*.pdb | wc | awk '{print $1}'`
echo $nconf $onepdb
ncoord=`cat $onepdb | grep ' CA ' | wc | awk '{print $1}'`

echo $onepdb $ncoord

((k=0)) ; for i in best_pdb/*.pdb ; do ((k++)) ; echo "XYZ" $k $i ; grep ' CA ' $i | awk '{print $7,$8, $9}' ; done > clustering/coor.xyz
for i in best_pdb/*.pdb ; do  echo $i ; done > clustering/xyz.pdblist

echo 'clustering/coor.xyz' > clustering/clus.in
echo 'clustering/xyz.tmp' >> clustering/clus.in
echo $nconf >> clustering/clus.in
echo $ncoord >> clustering/clus.in
echo 0 >> clustering/clus.in
echo 5 >> clustering/clus.in
echo 15.0 >> clustering/clus.in
echo clustering/cluster.dat >> clustering/clus.in
echo clustering/ccenter.dat >> clustering/clus.in

bin/cluster.x < clustering/clus.in > clustering/clus.out

rm -f clustering/xyz.tmp
rm -f clustering/clus.*.score.dat

((k=0))
for cmember in `cat clustering/cluster.dat | grep -v '#' | awk '{print $2}'`
do
((k++))

pdbname=`awk '(NR=='"$k"'){print }' clustering/xyz.pdblist`
framen=`echo $pdbname | tr -s '.' ' ' | tr -s '/' ' ' | awk '{print $2}'`
echo $framen
awk '($2=="'$framen'"){print}' best_pdb/scores.dat >> clustering/clus.$cmember.score.dat

echo $cmember $pdbname $k

echo "MODEL        1" >> clustering/clus.$cmember.pdb
cat $pdbname >> clustering/clus.$cmember.pdb
echo "ENDMDL" >> clustering/clus.$cmember.pdb


done
