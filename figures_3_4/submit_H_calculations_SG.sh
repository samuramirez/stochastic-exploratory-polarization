
#!/bin/bash

for gef in 100 200 300 400 500 600 700
#for gef in  50
do
for Nmeth in 40 80 120 160 
#for Nmeth in  80 
do
#for ktype in 1 
#do


sbatch -t 24:00:00 -N 1 -n 10 -o GEF"$gef"Nmeth"$Nmeth"ktype$ktype.out --wrap="matlab -nodisplay -nosplash -r figureHsims_grid\($gef,$Nmeth\)"

done
done
#done
