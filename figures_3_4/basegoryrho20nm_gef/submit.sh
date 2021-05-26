
#!/bin/bash



for nsim in {1..30}
do
#for gef in 700 600 500 400 300 200 
for gef in 50 100 
do

sbatch -t 100:00:00 -n 1 -o simulation"$nsim".out --wrap="~/cmake/bin/smoldyn 2dpolarity.cfg  --define output_mol_count=output_count_gef${gef}n${nsim}.txt --define output_mol_pos=output_pos_gef${gef}n${nsim}.txt --define nsim=$nsim --define totgef=$gef" 

done
done
