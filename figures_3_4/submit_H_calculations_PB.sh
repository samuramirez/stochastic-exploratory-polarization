
#!/bin/bash

for gef in  50 100
#for gef in  500
do
#for Nmeth in 40 80 120 


sbatch -t 24:00:00 -N 1 -n 10 -o GEF"$gef".out --wrap="matlab -nodisplay -nosplash -r figureHsims_PBsmoldyn\($gef\)"

done
