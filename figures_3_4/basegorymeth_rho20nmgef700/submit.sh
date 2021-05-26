
#!/bin/bash


gcc main.c functions.c parameters.c -lm -o run

#arguments: N,  nsim, ktype


for N in 20 40 80 160 
do

for nsim in {1..10}
do

for ktype in 0 1 4   
do
sbatch -t 100:00:00 -n 1 -o N"$N"ktype"$ktype".out --wrap="./run $N $nsim $ktype"

done
done	
done	
