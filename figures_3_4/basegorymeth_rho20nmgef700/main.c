
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
	//Paramters: N, nsim, ktype

clock_t begin = clock();

	//first parameter is not assigned and contains the name of the executable
if(argc != 4){
	printf("Invalid number of parameters\n");exit(1);}

	
int seed = time(NULL)+10*atoi(argv[2]);
srand(seed);

printf("Starting main, rand seed = %d\n",seed);

int i,j,counter,reaction,subvnext, species_jump, dir,subv_end,number,event,samp,totA;
double randnum,sum_a, sum_Drates;

struct Parameters p;

initialize_parameters(argv, &p);

char filename[100];
make_filename(p,filename);
strcat(filename,"sim");
strcat(filename,argv[2]);
strcat(filename,".dat");
printf("file name %s\n",filename);
printf("fast pf 	kmi[10]/Dtot is %f\n",p.kmi[10]/(p.D[4]+p.D[1]));
printf("fast complex 	kmi[9]/Dtot is %f\n",p.kmi[9]/(p.D[5]+p.D[0]));
printf("slow pf		kmi[4]/Dtot is %f\n",p.kmi[4]/(p.D[1]+p.D[3]));
printf("slow complex 	kmi[5]/Dtot is %f\n",p.kmi[5]/(p.D[2]+p.D[0]));
printf("basal gef 	kmi[2]/Dtot is %f\n",p.kmi[2]/(p.D[2]+p.D[3]));
printf("dx is %f, dx/sigma = %f\n", p.h, p.h/p.sigma  );


double t;
double interval = 300;

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

//SAMPLING
int Nsamp=1;
double ts[Nsamp]; //times
for(i=0;i<Nsamp;i++){ts[i]=-1;}
samp=0;     
int nsim;
double time_samp = 10;
double t_counter = time_samp;
int tot42t;
//Loop many simulations

for(nsim=0;nsim<Nsamp;nsim++){

//Initialize variables    
initialize(state,positionq,subvolumeq, timeq, rate,  &t,a,p);
//Initialize queue;
build_queue(positionq, timeq, subvolumeq,p); 
checkqueue(positionq,timeq,subvolumeq, p);    
//printf("first check queue\n");

FILE * fileConc;
fileConc=fopen(filename,"w");


////////////
//TIME ITERATIONS
////////////

while(t < interval){

	//if(rate[subvolumeq[0]][2]<=0){
	//printf("Last event %d\n",event);
	//printf("time =  %f\n",t);
	//printf("timeq[0] =  %f\n",timeq[0]);
	//checkqueue(positionq,timeq,subvolumeq, p);    
	//exit(1);
	//}

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
	
	if (t >= t_counter){
		t_counter = t_counter + time_samp;
		if (t_counter < t){printf("Small sampling time compared with time increcments\n");}
/*	
		tot42t=0;
		for(i=0; i<p.Nsubv; i++){
			tot42t+= state[i][0];
		}
		system("clear");

		for(i=0; i<p.Nsubv; i++){
			if(i%p.N==0 && i != 0){printf("\n");}
			printf("%d,",(state[i][0]+state[i][1]));
			//printf("%2.f,",(state[i][3]));
		}

		printf("t = %f, tot42t = %d\n",t,tot42t);
*/
		for(i=0; i<p.Nsubv; i++){fprintf(fileConc,"%d\t",state[i][0]+state[i][1]);}
		fprintf(fileConc,"\n");
/*
		for(i=0; i<p.Nsubv; i++){fprintf(fileConc,"%d\t",state[i][1]);}
		fprintf(fileConc,"\n");
		for(i=0; i<p.Nsubv; i++){fprintf(fileConc,"%d\t",state[i][2]);}
		fprintf(fileConc,"\n");
		for(i=0; i<p.Nsubv; i++){fprintf(fileConc,"%d\t",state[i][3]);}
		fprintf(fileConc,"\n");
		for(i=0; i<p.Nsubv; i++){fprintf(fileConc,"%d\t",state[i][4]);}
		fprintf(fileConc,"\n");
		for(i=0; i<p.Nsubv; i++){fprintf(fileConc,"%d\t",state[i][5]);}
		fprintf(fileConc,"\n");
*/
		
	}
		
	
}//end simulation 

fclose(fileConc);	

get_data(state,p);

//checkqueue(positionq,timeq,subvolumeq, p);    
//printf("last check queue\n");

ts[nsim]=t;
samp++;
//printf("t is %f\n",t);
} //end many simulations


printf(" t = %f interval = %f \n",t,interval);	


/*
//Print timeseries to file
FILE * fileA;
fileA=fopen(p.filename,"w");
for(i=0; i<samp; i++){
	fprintf(fileA,"%.8f\t%f\n",seriesA[0][i],seriesA[1][i]);
}
fclose(fileA);	
printf("Number of sampled points = %d\n",samp);
*/

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


