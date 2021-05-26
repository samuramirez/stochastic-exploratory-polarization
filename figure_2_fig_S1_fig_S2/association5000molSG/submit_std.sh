
#!/bin/bash


gcc main_std.c functions.c parameters.c -lm -o run_std

#for N in 10 
for N in  20 40 60 
do

for ktype in  0 1 4   
do

for nA0 in   5000
do

for kaDeff  in  50 
do

for kd  in 0 
do

#./run $N $ktype $nA0 $kaDeff $kd
sbatch -t 48:00:00 -n 1  --wrap="./run_std $N $ktype $nA0 $kaDeff $kd"

done
done	
done	
done	
done	
