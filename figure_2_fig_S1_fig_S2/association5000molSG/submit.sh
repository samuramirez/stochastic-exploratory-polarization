
#!/bin/bash


gcc main.c functions.c parameters.c -lm -o run

#for N in 10 
for N in  20 40 60 
do

for ktype in  4 
do

for nA0 in   5000
do

for kaDeff  in  50 
do

for kd  in 10 1 
do

#./run $N $ktype $nA0 $kaDeff $kd
sbatch -t 12:00:00 -n 1  --wrap="./run $N $ktype $nA0 $kaDeff $kd"

done
done	
done	
done	
done	
