Readme Figure 2.

Code to generate the figures is in .py files (python).
figs5molNS.py can generate Figures 2  A,B, E,F (5 molecules)
figs5000molNS.py can generate Figures 2  C,D, G,H (5 molecules)
figs5molNS_std.py can generate Figure S1  A,B, E,F (5 molecules)
figs5000molNS_std.py can generate Figure S1  C,D, G,H (5 molecules)

Code to generate spatial Gillespie raw data is in folders association5molSG (5 molecules data) and association5000molSG.
Within these folders, the code to run simulations is written in C. The .sh files contain the code to compile C files and submit jobs varying parameters to a slurm managed Linux cluster. 

Code to generate particle-based raw data is in reversiblePB. Each subfolder corresponds to a particular parameter set.
To run the simulations first the parameters are set in the Matlab file parameters_make_bash.m. Upon execution of this script, a .sh file will be generated to submit Matlab simulations to a slurm managed Linux cluster splitting the number of simulations in different jobs. Once the simulations are run, the several output files are assembled into a single file using the script assemble_and_check_files.m. This last script takes as input the folder containing the data, the metadata.mat file (one of the ouptput files of parameters_make_bash.m) and the name of the file containing the assembled data, which typically was "assembled_data.mat" (The python files to generate the figures read this last file). 