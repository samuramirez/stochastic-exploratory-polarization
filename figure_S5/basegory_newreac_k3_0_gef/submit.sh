
#!/bin/bash

gcc main.c functions.c parameters.c -lm -o run

interval=3600
samp=60
sigma=0.02
Dm=0.0045
for Dc in 10.0
do
for N in 80
do
for nsim in {1..50} 
do
for L in 8.0
do
for tot42  in 5000
do
for totgef  in 15 25 50 100 200 300 500
do
for kmi0  in 0.1 #k1a  
do
kmi1=10   #k1b
for kmi2  in 0.032  #k2a
do
for kmi3  in 0.63  #k2b
do
for kmi4  in 0 #k3 
do
for kmi5  in 2  #k4a 
do
for kmi6  in 10  #k4b
do
for kmi7  in 4  #k5a
do
for kmi8  in 6.5  #k5b
do
for kmi9  in 0.2 #k7
do
for kmi10  in  0.5 #k10
do

sbatch -t 58:00:00 -n 1 -o gef_"$totgef".out --wrap="./run $nsim $N $L $interval $samp $sigma $Dm $Dc $tot42 $totgef $kmi0 $kmi1 $kmi2 $kmi3 $kmi4 $kmi5 $kmi6 $kmi7 $kmi8 $kmi9 $kmi10 gef $totgef"

#echo $nsim $N $L $interval $samp $sigma $Dm $Dc $tot42 $totgef $kmi0 $kmi1 $kmi2 $kmi3 $kmi4 $kmi5 $kmi6 $kmi7 $kmi8 $kmi9 $kmi10 k8  $kmi10

done
done
done	
done	
done	
done	
done	
done	
done	
done	
done	
done	
done	
done	
done	
done
