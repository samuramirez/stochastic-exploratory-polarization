
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
	//see submit.sh for paramters

clock_t begin = clock();

	//first parameter is not assigned and contains the name of the executable
if(argc != 24){
	printf("Invalid number of parameters\n");exit(1);}

//printf("parameters passed %s %s %s %s\n",argv[1],argv[2],argv[3],argv[4]);
	
int seed = time(NULL)+10*atoi(argv[1]);
srand(seed);

printf("Starting main, rand seed = %d\n",seed);

int i,j,counter,reaction,subvnext, species_jump, dir,subv_end,number,event,totA;
double randnum,sum_a, sum_Drates;

struct Parameters p;

initialize_parameters(argv,argc, &p);

char filename[100];
make_filename(p,filename,argv);
printf("file name %s\n",filename);
printf("fast pf 	kmi[10]/Dtot is %f\n",p.kmi[10]/(p.D[4]+p.D[1]));
printf("fast complex 	kmi[9]/Dtot is %f\n",p.kmi[9]/(p.D[5]+p.D[0]));
printf("slow pf		kmi[4]/Dtot is %f\n",p.kmi[4]/(p.D[1]+p.D[3]));
printf("slow complex 	kmi[5]/Dtot is %f\n",p.kmi[5]/(p.D[2]+p.D[0]));
printf("basal gef 	kmi[2]/Dtot is %f\n",p.kmi[2]/(p.D[2]+p.D[3]));
printf("dx is %f, dx/sigma = %f\n", p.h, p.h/p.sigma  );

printf(" %f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n",p.kmi[0],p.kmi[1],p.kmi[2]/pi/p.sigma/p.sigma,p.kmi[3],p.kmi[4]/pi/p.sigma/p.sigma,p.kmi[5]/pi/p.sigma/p.sigma,p.kmi[6],p.kmi[7],p.kmi[8],p.kmi[9]/pi/p.sigma/p.sigma,p.kmi[10]/pi/p.sigma/p.sigma);

double t;

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

//Initialize variables    
initialize(state,positionq,subvolumeq, timeq, rate,  &t,a,p);
//Initialize queue;
build_queue(positionq, timeq, subvolumeq,p); 
checkqueue(positionq,timeq,subvolumeq, p);    
//printf("first check queue\n");
//SAMPLING
double t_counter = p.time_samp;

FILE * fileConc;
fileConc=fopen(filename,"w");

////////////
//TIME ITERATIONS
////////////
for(i=0;i<p.Nsubv;i++){get_rate_timeq(i,positionq, timeq, state,rate,t,a,p);}
build_queue(positionq, timeq, subvolumeq,p); 
checkqueue(positionq,timeq,subvolumeq, p);    

while(t < p.interval){
	
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

	//SAMPLING
	if (t >= t_counter){
		t_counter = t_counter + p.time_samp;
		if (t_counter < t){printf("Small sampling time compared with time increcments\n");}
/*	
		printf("t is %f\n",t);	
		system("clear");

		for(i=0; i<p.Nsubv; i++){
			if(i%p.N==0 && i != 0){printf("\n");}
			printf("%d,",(state[i][0]+state[i][1]));
			//printf("%2.f,",(state[i][3]));
		}
*/
		//Save data on file
		for(j=0;j<=1; j++){
			for(i=0; i<p.Nsubv; i++){fprintf(fileConc,"%d\t",state[i][j]);}
			fprintf(fileConc,"\n");
		}
	}
}//end simulation 

fclose(fileConc);	
checkqueue(positionq,timeq,subvolumeq, p);    
//printf("last check queue\n");

//save final concentration
get_data(state,p);

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

printf("Simulated time = %f, CPU time = %f\n",t,time_spent); 
return 0;
}//END MAIN


