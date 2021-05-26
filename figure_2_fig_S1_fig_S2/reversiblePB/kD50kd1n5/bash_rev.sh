#!/bin/bash

sbatch -p general -N 1 -n 12 -t 144:00:00 --wrap="matlab -nodisplay -nosplash -singleCompThread -r revheteroanni_parManager\(5,5,0.0025,0.0025,0.25,1,0.005,0.5,1e+06,1.000000e-05,100,1000,0,\'tempDir0001\',\'revAt0005_Bt0005_km0-25_kr1_DA0-0025_DB0-0025_realiz0001_to_1000\',12\) -logfile 0001.log"
sbatch -p general -N 1 -n 12 -t 144:00:00 --wrap="matlab -nodisplay -nosplash -singleCompThread -r revheteroanni_parManager\(5,5,0.0025,0.0025,0.25,1,0.005,0.5,1e+06,1.000000e-05,100,1000,1000,\'tempDir0002\',\'revAt0005_Bt0005_km0-25_kr1_DA0-0025_DB0-0025_realiz1001_to_2000\',12\) -logfile 0002.log"
sbatch -p general -N 1 -n 12 -t 144:00:00 --wrap="matlab -nodisplay -nosplash -singleCompThread -r revheteroanni_parManager\(5,5,0.0025,0.0025,0.25,1,0.005,0.5,1e+06,1.000000e-05,100,1000,2000,\'tempDir0003\',\'revAt0005_Bt0005_km0-25_kr1_DA0-0025_DB0-0025_realiz2001_to_3000\',12\) -logfile 0003.log"
sbatch -p general -N 1 -n 12 -t 144:00:00 --wrap="matlab -nodisplay -nosplash -singleCompThread -r revheteroanni_parManager\(5,5,0.0025,0.0025,0.25,1,0.005,0.5,1e+06,1.000000e-05,100,1000,3000,\'tempDir0004\',\'revAt0005_Bt0005_km0-25_kr1_DA0-0025_DB0-0025_realiz3001_to_4000\',12\) -logfile 0004.log"
sbatch -p general -N 1 -n 12 -t 144:00:00 --wrap="matlab -nodisplay -nosplash -singleCompThread -r revheteroanni_parManager\(5,5,0.0025,0.0025,0.25,1,0.005,0.5,1e+06,1.000000e-05,100,1000,4000,\'tempDir0005\',\'revAt0005_Bt0005_km0-25_kr1_DA0-0025_DB0-0025_realiz4001_to_5000\',12\) -logfile 0005.log"
sbatch -p general -N 1 -n 12 -t 144:00:00 --wrap="matlab -nodisplay -nosplash -singleCompThread -r revheteroanni_parManager\(5,5,0.0025,0.0025,0.25,1,0.005,0.5,1e+06,1.000000e-05,100,1000,5000,\'tempDir0006\',\'revAt0005_Bt0005_km0-25_kr1_DA0-0025_DB0-0025_realiz5001_to_6000\',12\) -logfile 0006.log"
sbatch -p general -N 1 -n 12 -t 144:00:00 --wrap="matlab -nodisplay -nosplash -singleCompThread -r revheteroanni_parManager\(5,5,0.0025,0.0025,0.25,1,0.005,0.5,1e+06,1.000000e-05,100,1000,6000,\'tempDir0007\',\'revAt0005_Bt0005_km0-25_kr1_DA0-0025_DB0-0025_realiz6001_to_7000\',12\) -logfile 0007.log"
sbatch -p general -N 1 -n 12 -t 144:00:00 --wrap="matlab -nodisplay -nosplash -singleCompThread -r revheteroanni_parManager\(5,5,0.0025,0.0025,0.25,1,0.005,0.5,1e+06,1.000000e-05,100,1000,7000,\'tempDir0008\',\'revAt0005_Bt0005_km0-25_kr1_DA0-0025_DB0-0025_realiz7001_to_8000\',12\) -logfile 0008.log"
sbatch -p general -N 1 -n 12 -t 144:00:00 --wrap="matlab -nodisplay -nosplash -singleCompThread -r revheteroanni_parManager\(5,5,0.0025,0.0025,0.25,1,0.005,0.5,1e+06,1.000000e-05,100,1000,8000,\'tempDir0009\',\'revAt0005_Bt0005_km0-25_kr1_DA0-0025_DB0-0025_realiz8001_to_9000\',12\) -logfile 0009.log"
sbatch -p general -N 1 -n 12 -t 144:00:00 --wrap="matlab -nodisplay -nosplash -singleCompThread -r revheteroanni_parManager\(5,5,0.0025,0.0025,0.25,1,0.005,0.5,1e+06,1.000000e-05,100,1000,9000,\'tempDir0010\',\'revAt0005_Bt0005_km0-25_kr1_DA0-0025_DB0-0025_realiz9001_to_10000\',12\) -logfile 0010.log"