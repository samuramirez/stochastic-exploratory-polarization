#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "header.h"


void diffuse(double randnum, int subvnext, int ** state, double ** rate, double ** a, double * timeq, int * positionq, int * subvolumeq, double t, int ** connect, struct Parameters p)
{	
	int j,i,species_jump,dir, subv_end;
	double sum_Drates; 
	//rescale the random number to  [0,1]

	//at this point  rate[subvnext][0]/rate[subvnext][2] <=   randnum <= 1 
//	randnum = (randnum - rate[subvnext][0]/rate[subvnext][2])/(1 - rate[subvnext][0]/rate[subvnext][2]);
    
	//determine which type of molecule diffuse
	randnum= (double)rand()/RAND_MAX;
	while(randnum==0.0){randnum= (double)rand()/RAND_MAX;}
	sum_Drates=0;
	for(j=p.Nsp-1; j>=0; j--){//counting backwards to get cytosolic species first
		sum_Drates = sum_Drates + 4*p.D[j]*state[subvnext][j]/(p.h*p.h)  ;
		if(randnum <= sum_Drates/rate[subvnext][1]){ species_jump=j; break;}
	}
    if(j==-1){
		printf("species to jump not assigned\n"); 
		printf("normalized random number %f\n",randnum);
		printf("total propensity to jump at subvnext %f\n",rate[subvnext][1]);
		printf("total propensity to react %f\n",rate[subvnext][0]);
		printf("sum of propensities to jump %f\n",sum_Drates);
		printf("Number of molecules in subvnext %d %d\n",state[subvnext][0],state[subvnext][1]);
		exit(1);
	}

	randnum= (double)rand()/RAND_MAX;
	//determine direction of movement
	for(j=1; j<=4; j++){
		if(randnum <= j/4.0){ dir = j; break;}
	}
	
    if(j==5){printf("direction of movement not assigned\n"); exit(1);}
	//subvolume where particle jumped (connect direction index goes from 0 to 3)
	subv_end = connect[subvnext][dir-1];
	state[subvnext][species_jump]--;
	state[subv_end][species_jump]++;
		
	if(state[subvnext][species_jump]<0){printf("species less than zero after diffusing!\n"); exit(1);}
    
	//Update rate and time to next reaction for subvolum subvend
    get_rate_timeq(subv_end,positionq,timeq,state,rate,t,a,p);
	sort(subv_end,positionq,timeq,subvolumeq,p);
    //Update rate and time to next reaction for subvolume next
    get_rate_timeq(subvnext,positionq,timeq,state,rate,t,a,p);
    //update queue
    sort(subvnext,positionq,timeq,subvolumeq,p);
	
}


void react(double randnum, int subvnext, int ** state, double ** rate, double ** a, double * timeq, int * positionq, int * subvolumeq, double t, struct Parameters p){

	int j,reaction;
	double sum_a;

//rescale the random number to  [0,1]
// At this point 0 < randnum <= rate[subvnext][0]/rate[subvnext][2]
//randnum = randnum/(rate[subvnext][0]/rate[subvnext][2]);

randnum= (double)rand()/RAND_MAX;
//avoid randnum==0 to avoid selecting reaction j==0  when randnum==0  
while(randnum==0.0){randnum= (double)rand()/RAND_MAX;}

//Choose which reaction occurs
sum_a=0;
for(j=0; j<p.Nreact; j++){
	
	sum_a=sum_a + a[subvnext][j] ;
	if(randnum <= sum_a/rate[subvnext][0]){ reaction=j; break;}
}
if(j==p.Nreact){
	printf("reaction not assigned\n");
	printf("normalized random number %f\n",randnum);
	printf("total propensity %f\n",rate[subvnext][0]);
	printf("sum of propensities %f\n",sum_a);
	printf("difrate + reacrate  %f\n",rate[subvnext][2]);
	printf("Number of molecules in subvnext %d %d\n",state[subvnext][0],state[subvnext][1]);
	exit(1);
}


//if (reaction==9){printf(" species1 = %d, species2 = %d, kmeso=%f\n",state[subvnext][0],state[subvnext][5],kMesoHetero(9,5,0,subvnext,state,p));}

update_states_reaction(reaction, state, subvnext);
//Update rate and time to next reaction for subvolume next
get_rate_timeq(subvnext,positionq,timeq,state,rate,t,a,p);
//update queue
sort(subvnext,positionq,timeq,subvolumeq,p);

}



void allocate(int ** state, double ** rate, double ** a, double * timeq, int * positionq, int * subvolumeq, int ** connect, struct Parameters p){

int i;
	
//ALLOCATION
//Creat state matrix. Each row corresponds to a subvolume. The columns indicate the number of each particular species.
state = (int **)malloc(p.Nsubv*sizeof(int *));
for(i=0;i<p.Nsubv;i++){state[i]=(int *)malloc(p.Nsp*sizeof(int));}

//First column is the sum of reaction rates, second column is the sum of diffusion rates, third column is the sum of diffusion and reaction rates
rate = (double **)malloc(p.Nsubv*sizeof(double *));
for(i=0;i<p.Nsubv;i++){rate[i]=(double *)malloc(3*sizeof(double));}

//propensities for each reactions for each subvolume
a = (double **)malloc(p.Nsubv*sizeof(double *));
for(i=0;i<p.Nsubv;i++){a[i]=(double *)malloc(p.Nsubv*sizeof(double));}

//time queue    
positionq=(int *)malloc(p.Nsubv*sizeof(int));
subvolumeq=(int *)malloc(p.Nsubv*sizeof(int));
timeq=(double *)malloc(p.Nsubv*sizeof(double));

//connectivity matrix
connect = (int **)malloc(p.Nsubv*sizeof(int *));
for(i=0;i<p.Nsubv;i++){connect[i]=(int *)malloc((2*dim)*sizeof(int));}
	
}

void swap(int positionq1, int positionq2,int *  positionq,double *  timeq, int * subvolumeq){

    int subvolumeq1,subvolumeq2,subvolumeq_temp, positionq_temp;
    double time_temp;
    
//get subvolume indices to later swap positions (not necessary but the code is clearer) 	
subvolumeq1 = subvolumeq[positionq1];
subvolumeq2 = subvolumeq[positionq2];

//swap values in the queue
time_temp	 		= timeq[positionq1];
timeq[positionq1]	= timeq[positionq2];
timeq[positionq2]	= time_temp;

subvolumeq_temp			= subvolumeq[positionq1];
subvolumeq[positionq1]	= subvolumeq[positionq2];
subvolumeq[positionq2]	= subvolumeq_temp;
		
//swap values in the "positions in queue" array
positionq_temp			= positionq[subvolumeq1];
positionq[subvolumeq1]	= positionq[subvolumeq2];
positionq[subvolumeq2]	= positionq_temp;
}//end swap function

void checkqueue(int * positionq, double * timeq, int * subvolumeq, struct Parameters p){
	int i,parent,left,right;
	for(i=0;i<p.Nsubv; i++){
		parent = ceil(i/2.0) - 1;
		left = i*2+1;
		right = i*2+2;
		
		if(parent >= 0 && timeq[parent] > timeq[i]){
			printf("parent %f greater than child %f\n",timeq[parent],timeq[i]);
			exit(1);
		}
		if(left < p.Nsubv && timeq[i] > timeq[left]){
			printf("parent %f %d greater than left %f %d\n",timeq[i],i,timeq[left],left);
			exit(1);
		}
		if(right < p.Nsubv && timeq[i] > timeq[right]){
			printf("parent %f %d greater than right %f %d\n",timeq[i],i,timeq[right],right);
			
			exit(1);
		}		
	}
	
}


void sort (int subv, int * positionq, double * timeq, int * subvolumeq, struct Parameters p){
	int parent_pos, daughter1, daughter2, done_bub_down;
	
//get position of parent in queue		
parent_pos = (positionq[subv]+1)/2 - 1;
	
if(parent_pos >=0 && timeq[parent_pos] > timeq[positionq[subv]] ){

	//While time < parent time,  bubble up
	while(positionq[subv]>0 && timeq[positionq[subv]] < timeq[parent_pos]  ){
	//swap
		swap(positionq[subv],parent_pos,positionq,timeq,subvolumeq);
		//update parent_pos	
		parent_pos = (positionq[subv]+1)/2 - 1;		
	} //while bubbling up
}
	
else{
	heapify(positionq[subv],positionq,timeq,subvolumeq,p);
}	
	
	
}//end sort function


void sort_no_heapify (int subv, int * positionq, double * timeq, int * subvolumeq, struct Parameters p){
	int parent_pos, daughter1, daughter2, done_bub_down;
	
//get position of parent in queue		
parent_pos = (positionq[subv]+1)/2 - 1;
	
if(parent_pos >=0 && timeq[parent_pos] > timeq[positionq[subv]] ){

	//While time < parent time,  bubble up
	while(positionq[subv]>0 && timeq[positionq[subv]] < timeq[parent_pos]  ){
	//swap
		swap(positionq[subv],parent_pos,positionq,timeq,subvolumeq);
		//update parent_pos	
		parent_pos = (positionq[subv]+1)/2 - 1;		
	} //while bubbling up
}
	
else{
	
	//BUBBLE DOWN
	
	done_bub_down=0;

	while (done_bub_down == 0){

	//update daughters
	daughter1 = positionq[subv]*2+1;
	daughter2 = positionq[subv]*2+2;

	//if has both daughters
	if  (daughter2 < p.Nsubv ){
		//if time is larger than any of the daughters swap
		if (timeq[positionq[subv]] >  timeq[daughter1] ||  timeq[positionq[subv]] > timeq[daughter2]){

			//	swap with daughter with smaller time
			if(timeq[daughter1] < timeq[daughter2]){
				swap(positionq[subv], daughter1,positionq,timeq,subvolumeq);	}
			else { 
				swap(positionq[subv], daughter2,positionq,timeq,subvolumeq);}

		}
		//time is smaller than both of the daughters
		else { done_bub_down =1; }  
	}

	//doesn't have daughter2 but has daughter1
	else if(daughter1 < p.Nsubv){

		if (timeq[subv] > timeq[daughter1] ){ 
			swap(positionq[subv], daughter1,positionq,timeq,subvolumeq);
		}
		//done bubbling  either case
		done_bub_down =1;
	}

	//no daughters
	else{done_bub_down =1;}
		 
	}//while bubbling down
		
	
	
}	
	
	
}//end sort function


//heapify takes a node that may not root its sub-heap-tree and updates the heap-tree assuming that the daughter nodes are roots of their respective sub-heap-trees 
void heapify(int i, int *positionq,double * timeq, int * subvolumeq, struct Parameters p){

int left,right,smallest;

	left = i*2+1;	
	right = i*2+2;

	if(left < p.Nsubv && timeq[left] < timeq[i] ){
	smallest = left;
	}
	else {smallest = i; }
	

	if(right < p.Nsubv && timeq[right] < timeq[smallest] ){
	smallest = right;
	}

	if(smallest != i){
        
        swap(i,smallest, positionq, timeq, subvolumeq);
		heapify(smallest,positionq,timeq,subvolumeq,p);
	}

}

void build_queue(int *positionq,double * timeq,int * subvolumeq, struct Parameters p){
	int i;

    for(i=p.Nsubv/2; i >=0; i--){
        heapify(i,positionq, timeq, subvolumeq,p);
	}
}


void create_connect(int ** connect, struct Parameters p){

//Create connectivity matrix. Each row corresponds to the index of a subvolume. The indices of the neighbors are stored in the columns. Indices increase according to Matlab matrices indexing: [0][0] is the upper left corner of the matrix, increasing index goes down in the y direction.

		
int x,y,xtemp,ytemp,i;
    
for(i=0; i < p.Nsubv ; i++){

    //subv = y + x*N
    
	y = i%p.N;
	x = i/p.N;

	connect[i][0] = p.N*x + (y+1)%p.N ;

	connect[i][1] = p.N*((x+1)%p.N) + y ;

		ytemp = y - 1 < 0 ? p.N -1 : y - 1 ;
	connect[i][2] = p.N*x + ytemp  ;

		xtemp = x - 1 < 0 ? p.N - 1 : x - 1;
	connect[i][3] = p.N*xtemp + y ;
    
    
}
    
}


