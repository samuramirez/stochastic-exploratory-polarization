
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
///////////
double pert_time=60;
int pert=0; 
int xi,yi,subv;
double frapL=3;
int tot[p.Nsp];

while(t < p.interval){
	

 //Bleach Cdc42 at the patch
	if(t > pert_time && pert==0){
	//convert all Cdc42 at the patch into the dark versions
	//xi and yi are the grid coordinates where bleaching takes place
		for(xi=floor(p.N/2 - frapL/2.0/p.h); xi <  floor(p.N/2 + frapL/2.0/p.h); xi++ ){
			for(yi=floor(p.N/2 - frapL/2.0/p.h); yi <  floor(p.N/2 + frapL/2.0/p.h); yi++ ){
			subv = xi + p.N*yi;
	//
	// Bleach all Cdc42 at the membrane
	//	for(subv=0;subv<p.Nsubv;subv++){
			//convert BemGEF42 into BemGEF42d
			state[subv][6]= state[subv][1];
			state[subv][1]= 0;
			//convert Cdc42T into Cdc42Td
			state[subv][7]=state[subv][0];
			state[subv][0]=0;
			//convert Cdc42D into Cdc42Dd
			state[subv][8]= state[subv][3];
			state[subv][3]= 0;

			get_rate_timeq(subv,positionq,timeq,state,rate,t,a,p);
			sort(subv,positionq,timeq,subvolumeq,p);
	//		}

			}
		}
	pert=1;
 	}
 


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
		
		for(j=0; j<p.Nsp; j++){
			tot[j]=0;
			//look at all the domain
			for(i=0; i<p.Nsubv; i++){
	   		tot[j]+=state[i][j];
			}		
			//look only at the bleached region
/*			for(xi=floor(p.N/2 - frapL/2.0/p.h); xi <  floor(p.N/2 + frapL/2.0/p.h); xi++ ){
				for(yi=floor(p.N/2 - frapL/2.0/p.h); yi <  floor(p.N/2 + frapL/2.0/p.h); yi++ ){
					i = xi + p.N*yi;
   				tot[j]+=state[i][j];

			}
		}
*/
		}
	//	system("clear");

	//	for(i=0; i<p.Nsubv; i++){
	//		if(i%p.N==0 && i != 0){printf("\n");}
	//		printf("%d,",(state[i][2]+state[i][1]));
			//printf("%2.f,",(state[i][3]));
//		}

//		printf("t is %f tot42 is %d tot gef42 is %d\n",t,tot42,totgef);
	
		//Save data on file
		for(j=0; j<p.Nsp; j++){ fprintf(fileConc,"%d\t",tot[j]); }
		fprintf(fileConc,"%f",t);
		fprintf(fileConc,"\n");
		//
//		for(i=0; i<p.Nsubv; i++){fprintf(fileConc,"%d\t",state[i][0]);	}
//		fprintf(fileConc,"\n");
//		for(i=0; i<p.Nsubv; i++){fprintf(fileConc,"%d\t",state[i][1]);	}
//		fprintf(fileConc,"\n");
//		for(i=0; i<p.Nsubv; i++){fprintf(fileConc,"%d\t",state[i][2]);	}
//		fprintf(fileConc,"\n");
//		for(i=0; i<p.Nsubv; i++){fprintf(fileConc,"%d\t",state[i][3]);	}
//		fprintf(fileConc,"\n");
//		for(i=0; i<p.Nsubv; i++){fprintf(fileConc,"%d\t",state[i][6]);	}
//		fprintf(fileConc,"\n");
//		for(i=0; i<p.Nsubv; i++){fprintf(fileConc,"%d\t",state[i][7]);	}
//		fprintf(fileConc,"\n");
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


