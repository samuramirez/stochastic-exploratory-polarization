
//  Created by Samuel Ramirez on 11/3/17.
//  Copyright (c) 2017 Samuel Ramirez. All rights reserved.
//

#include "header.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

int main(int argc, char *argv[]){
	//Paramters: N, ktype, nA=nB, ka/Deff, kd   

clock_t begin = clock();

	//first parameter is not assigned and contains the name of the executable
if(argc != 6){
	printf("Invalid number of parameters\n");exit(1);}

//printf("parameters passed %s %s %s %s\n",argv[1],argv[2],argv[3],argv[4]);
	
//int seed = time(NULL)+10*atoi(argv[2]);
int seed = time(NULL);
srand(seed);

printf("Starting main, rand seed = %d\n",seed);

int i,j,counter,reaction,subvnext, species_jump, dir,subv_end,number,event,samp,totA;
double randnum,sum_a, sum_Drates;

struct Parameters p;

initialize_parameters(argv, &p);

char filename[100];
make_filename(p,filename);
printf("dx is %f, dx/sigma = %f, ka is %.6f\n, kaDeff is %f, Dm is %f\n ", p.h, p.h/p.sigma,p.kmi[0],p.kaDeff,p.Dm );

printf("output name %s\n",filename);

//MEMORY ALLOCATION
int ** state, * positionq, * subvolumeq, ** connect; 
double ** rate, ** a, * timeq; 
//allocate(state, rate, a, timeq, positionq, subvolumeq, connect, p);
state = (int **)malloc(p.Nsubv*sizeof(int *));
for(i=0;i<p.Nsubv;i++){state[i]=(int *)malloc(p.Nsp*sizeof(int));}
//First column is the sum of reaction rates, second column is the sum of diffusion rates, third column is the sum of diffusion and reaction rates
rate = (double **)malloc(p.Nsubv*sizeof(double *));
for(i=0;i<p.Nsubv;i++){rate[i]=(double *)malloc(3*sizeof(double));}
//propensities for each reactions for each subvolume
a = (double **)malloc(p.Nsubv*sizeof(double *));
for(i=0;i<p.Nsubv;i++){a[i]=(double *)malloc(p.Nreact*sizeof(double));}
//time queue    
positionq=(int *)malloc(p.Nsubv*sizeof(int));
subvolumeq=(int *)malloc(p.Nsubv*sizeof(int));
timeq=(double *)malloc(p.Nsubv*sizeof(double));
//connectivity matrix
connect = (int **)malloc(p.Nsubv*sizeof(int *));
for(i=0;i<p.Nsubv;i++){connect[i]=(int *)malloc((2*dim)*sizeof(int));}

//Create connectivity matrix. 
create_connect(connect,p);

double t;
int Nsims=100;
double interval = 5;
int nsim;

//SAMPLING
double time_samp = 0.001; 
int Nsamp = interval/time_samp + 1; // +1 because t=0 is recorded
printf("sampling every %f seconds, %d samples\n",time_samp,Nsamp);
double t_counter;
int totstate[p.Nsp];
int counts=4000; // 4seg/0.001
int icounts=0;
int histogram[Nsims*counts];
double thistogram[Nsims*counts];
double eqtime = 1;
double  meanstate[Nsamp][p.Nsp];
double meanstate2[Nsamp][p.Nsp];
int  sampstate[Nsamp][p.Nsp];
for(i=0;i<Nsamp;i++){
	for(j=0;j<p.Nsp;j++){meanstate[i][j]=0; meanstate2[i][j]=0;}
}
for(i=0;i<Nsims*counts;i++){histogram[i]=-1; thistogram[i]=-1;}
int itime;
//Loop many simulations

for(nsim=0;nsim<Nsims;nsim++){

	//printf("sims %d\n",nsim);
	

//Initialize variables    
initialize(state,positionq,subvolumeq, timeq, rate,  &t,a,p);
//Initialize queue;
build_queue(positionq, timeq, subvolumeq,p); 
checkqueue(positionq,timeq,subvolumeq, p);    

//sample t=0
itime=0;
t_counter=0;
for(j=0;j<p.Nsp;j++){sampstate[itime][j] = p.n0[j];	}	
t_counter +=  time_samp;
itime++;

////////////
//TIME ITERATIONS
////////////

while(t < interval){

    //Update time when next event occurs
    t = timeq[0];
	//get subvolume where next event occurs
	subvnext=subvolumeq[0];
   
	randnum= (double)rand()/RAND_MAX;
	//avoid randnum==0 to avoid trying a reaction when randnum==0 and rate[subvnext][0] == 0 
	while(randnum==0.0){randnum= (double)rand()/RAND_MAX;}	
	//reaction event
	
	if(randnum <= rate[subvnext][0]/rate[subvnext][2] ){
		event=0;
		react(randnum, subvnext, state, rate, a, timeq, positionq, subvolumeq, t, p);	
//		checkqueue(positionq,timeq,subvolumeq, p)				
	}// end reaction event

	//diffusion event
	else{ //  rate[subvnext][0]/rate[subvnext][2] < randnum < 1		
		event=1;
		diffuse(randnum, subvnext, state, rate, a, timeq,  positionq, subvolumeq, t, connect, p);
//		checkqueue(positionq,timeq,subvolumeq, p);   		
	}//end diffusion event
	
	//sampling
	////event at time slightly larger or equal than t_counter is stored at t_counter
	if (t >= t_counter){ 

		//if  next t_counter does not surpass t
		if (t_counter + time_samp < t){
			printf("next tcounter less than t");
			exit(1);
			// while t_counter <= t we copy the previous sample 
			while(t_counter <= t && itime < Nsamp){

				for(j=0;j<p.Nsp;j++){sampstate[itime][j] = sampstate[itime-1][j];}
				itime++;
				t_counter += time_samp;
			}

			//now t_counter > t and we sample on the previous itime
			for(j=0;j<p.Nsp;j++){
				totstate[j]=0;
				for(i=0; i<p.Nsubv; i++){totstate[j] += state[i][j];}
				sampstate[itime-1][j] = totstate[j];
			}
		
		} //end small sampling time compared with time increment
		
		else{ //next t_counter greater than t
			for(j=0;j<p.Nsp;j++){
				totstate[j]=0;
				for(i=0; i<p.Nsubv; i++){totstate[j] += state[i][j];}
				sampstate[itime][j] = totstate[j];
			}
			t_counter += time_samp;
			itime++;
			//get counts for histrogram
			if (icounts < Nsims*counts && t > eqtime){
				histogram[icounts]=totstate[0];
				thistogram[icounts]=t;
				icounts++;}

		}

	//printf("time sampled = %f, tot A = %d\n",t, totstate[0]);		
	}//end sampling
		
}//end simulation 

//store states and states squared in mean containers
for(i=0;i<Nsamp;i++){
	for(j=0;j<p.Nsp;j++){
		meanstate[i][j]+=sampstate[i][j];
		meanstate2[i][j]+=pow(sampstate[i][j],2);
	 } 
}

checkqueue(positionq,timeq,subvolumeq, p);    

} //end many simulations

//compute means dividing by the number of simulations
for(i=0;i<Nsamp;i++){
	for(j=0;j<p.Nsp;j++){
		meanstate[i][j]/=Nsims;
		meanstate2[i][j]/=Nsims;
	 } 
}

FILE * fileTseries;

char mean_file_name[100];
strcpy(mean_file_name,filename);
strcat(mean_file_name,"hist.dat");
fileTseries=fopen(mean_file_name,"w");

for(i=0;i<icounts;i++){
	fprintf(fileTseries,"%d,%f", (int) histogram[i],(double) thistogram[i]);
	fprintf(fileTseries,"\n");
}
fclose(fileTseries);


//Free memory
for (i=0; i<p.Nsubv; i++) {free(state[i]);}  
free(state);
for (i=0; i<p.Nsubv; i++) {free(rate[i]);}  
free(rate);
for (i=0; i<p.Nsubv; i++) {free(connect[i]);}  
free(connect);
free(timeq);free(positionq);free(subvolumeq);

clock_t end = clock();
double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

printf("End main CPU time = %f\n",time_spent); 
return 0;
}//END MAIN


