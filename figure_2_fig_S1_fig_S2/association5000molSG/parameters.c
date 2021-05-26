#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "header.h"
#include <string.h>

double kMesoHetero(int reaction,int species1,int species2, int subv,int ** state, struct Parameters p) {
	//Make sure "reaction" is an association reaction
	//kmi units are um^2/s
	//this function returns kmeso/h^2 so units are 1/s

int j,n;
double kmeso,byj,bsr,lambda;

double Deff = p.D[species1] + p.D[species2];


//npart = largest between species1  and species2
n = state[subv][species1] >  state[subv][species2] ? state[subv][species1] : state[subv][species2] ;

	
//free area equality	
//pi*bsr^2 - pi*sigma^2 = h*h/n - pi*sigma^2 

if(p.ktype<0 ||p.ktype > 4){printf("wrong ktype\n"); exit(1);}

switch(p.ktype){
	
	case 0:
		kmeso=p.kmi[reaction];	
		break;
	case 1:
		//Hellander-Petzold
		kmeso = p.kmi[reaction]/(1 + p.kmi[reaction]/(2*pi*Deff)*( log(p.h/sqrt(pi)/p.sigma) - (3 + 0.1951*2*pi)/4 ) );
		break;
	case 2:
		//Yogurtcu-Johnson
		byj = 2*sqrt(p.h*p.h/pi/n + p.sigma*p.sigma);
		lambda = p.sigma/byj;		
		kmeso = p.sigma < byj ? p.kmi[reaction]/(1 + p.kmi[reaction]/(2*pi*Deff)*( log(1/lambda)/pow(1-lambda*lambda,2) - (3-lambda*lambda)/(1-lambda*lambda)/4 ) ) : p.kmi[reaction];
		break;
	case 3:	
		//Ours homo
		bsr = sqrt(2*p.h*p.h/pi/n);
		lambda = p.sigma/bsr;		
		kmeso = p.sigma < bsr ? p.kmi[reaction]/(1 + p.kmi[reaction]/(2*pi*Deff)*( log(1/lambda)/pow(1-lambda*lambda,2) - (3-lambda*lambda)/(1-lambda*lambda)/4 ) ) : p.kmi[reaction];	
		break;
	case 4:	
		//Ours hetero
		if(n==0){kmeso=p.kmi[reaction];}
		else{
		bsr = sqrt(p.h*p.h/pi/n);
		lambda =  p.sigma/bsr ;
		kmeso = p.sigma < bsr ? p.kmi[reaction]/(1 + p.kmi[reaction]/(2*pi*Deff)*( log(1/lambda)/pow(1-lambda*lambda,2) - (3-lambda*lambda)/(1-lambda*lambda)/4 ) ) : p.kmi[reaction];}
		break;
	case 5:	
		//Ours hetero
		if(n==0){n=1;}//to compute a reasonalbe kdmeso in low concentration situations
		bsr = sqrt(p.h*p.h/pi/n);
		lambda =  p.sigma/bsr ;
		kmeso = p.sigma < bsr ? p.kmi[reaction]/(1 + p.kmi[reaction]/(2*pi*Deff)*( log(1/lambda)/pow(1-lambda*lambda,2) - (3-lambda*lambda)/(1-lambda*lambda)/4 ) ) : p.kmi[reaction];	
		break;
}


if (p.kmi[reaction]==0.0){kmeso=0;}
if (kmeso<0){printf("kmeso = %f, lambda = %f\n", kmeso, lambda); exit(1);}
if (p.ktype ==4 && kmeso > p.kmi[reaction] ){printf("kmeso = %f greater than  kmi = %f , lambda = %f, reaction %d\n", kmeso, p.kmi[reaction], lambda, reaction); exit(1);}
return kmeso/(p.h*p.h);

}


void get_rate_timeq(int subv, int * positionq, double * timeq, int ** state, double ** rate, double t, double ** a, struct Parameters p){

int j;

//get rates for subvolume subv

//Get propensities for each reaction
//
//A + B --> C
//C --> A + B
//0 A
//1 B
//2 C
a[subv][0] = kMesoHetero(0,0,1,subv,state,p)*state[subv][0]*state[subv][1];

//REVERSIBLE REACTION kmicro/dmicro=kmeso/dmeso  
//dmeso = dmicro*kmeso*h*h/kmicro

a[subv][1] = p.kmi[0] > 0 ? p.kmi[1]*kMesoHetero(0,0,1,subv,state,p)*p.h*p.h/p.kmi[0]*state[subv][2] : p.kmi[1]*state[subv][2];
	
//Get total rate as sum of propensities 
rate[subv][0]=0;
for(j=0; j< p.Nreact ; j++){rate[subv][0]+=a[subv][j];}

//sum of diffusion rates. 
rate[subv][1]=0;
for(j=0; j<p.Nsp; j++){rate[subv][1]+=2*dim*p.D[j]/(p.h*p.h)*state[subv][j];}

//total rate (reaction+diffusion)
rate[subv][2] =  rate[subv][0]+rate[subv][1];

double randnumtime =  (double)rand()/RAND_MAX ;

//-log(1.0)/0.0 == nan

while( randnumtime == 1.0){randnumtime =  (double)rand()/RAND_MAX;}

timeq[positionq[subv]]= t + -log(randnumtime)/(rate[subv][2]);

}

void update_states_reaction(int reaction, int ** state, int subvnext){

	//Update according to reaction
switch(reaction){
	case 0:
	//A + B --> C
		state[subvnext][0]--;
		state[subvnext][1]--;
		state[subvnext][2]++;
		break;
	case 1:
	//C --> A + B
		state[subvnext][0]++;
		state[subvnext][1]++;
		state[subvnext][2]--;
		break;
	default:
		printf("invalid reaction type");
		exit(1);
}//end update according to  reaction
	
}

void initialize(int ** state, int * positionq, int * subvolumeq, double * timeq, double ** rate, double * t, double ** a, struct Parameters p)
{

int subvi,i,j;

*t=0.0;

//Set species randomly on the grid

for(j=0;j<p.Nsp;j++){

	for(i=0; i<p.Nsubv ; i++ ){state[i][j]=0;}
	for(i=0; i < p.n0[j] ; i++ ){subvi=rand()%p.Nsubv; state[subvi][j]++;}

}


//INITIALIZE TRIVIALLY subvolumeq[positionq], timeq, positioniq[subvolume]
for(i=0;i<p.Nsubv;i++){positionq[i]=i; subvolumeq[positionq[i]]=i; }

//CALCULATE SUM OF REACTION AND DIFFUSION RATES AND TIME FOR NEXT REACTION
for(i=0;i<p.Nsubv;i++){get_rate_timeq(i,positionq, timeq, state,rate,*t,a,p);}


}

void initialize_parameters(char * argv[],struct Parameters * p){
	//Future implementation should read parameters and system from an input text file
strcpy(p->filename,"rev");

	
p->N = atoi(argv[1]); //Number of voxels on the side
p->Nsubv = pow(p->N,dim); //number of voxels
p->L = 1; //box length um
p->h = p->L/p->N; //grid size

p->sigma = 0.005; //um

p->ktype = atoi(argv[2]); //0=ka, 1=kh, 2=kyj, 3=ksrhomo, 4=ksrhet

//number of species
p->Nsp = 3; 
//Number of reactions
p->Nreact = 2;

p->Dm = 0.0025; //um2s-1
p->Dc = 1.0; //um2s-1

p->D=(double *)malloc(p->Nsp*sizeof(double));
p->D[0]=p->Dm;
p->D[1]=p->Dm;
p->D[2]=p->Dm;
	

//Initial number of proteins
p->n0=(int *)malloc(p->Nsp*sizeof(int));
p->n0[0]= atoi(argv[3]); 
p->n0[1]= atoi(argv[3]);
p->n0[2]= 0;
				

p->kaDeff=atof(argv[4]);		
p->kmi=(double *)malloc(p->Nreact*sizeof(double));
//for association reactions kmi = mike's lambda * pi*sigma^2
p->kmi[0]=p->kaDeff*2*p->Dm ;
p->kmi[1]=atof(argv[5]);		
	
						
}

void make_filename(struct Parameters p, char filename[]){

//Making filename
char  gridN[10];
sprintf(gridN, "%d", p.N);

if(p.ktype<0 ||p.ktype > 4){printf("wrong ktype\n"); exit(1);}
if(p.ktype==0){strcat(p.filename, "ka");}
if(p.ktype==1){strcat(p.filename, "kh");}
if(p.ktype==2){strcat(p.filename, "kyj");}
if(p.ktype==3){strcat(p.filename, "ksr");}
if(p.ktype==4){strcat(p.filename, "ksrhet");}

char factorD[10];
if(p.kaDeff<1){sprintf(factorD, "%.2f", p.kaDeff);}
else{sprintf(factorD, "%.0f", p.kaDeff);}

char kd[10];
if(p.kmi[1]<1){sprintf(kd, "%.1f", p.kmi[1]);}
else{sprintf(kd, "%.0f", p.kmi[1]);}

char  nini[p.Nsp][10];
sprintf(nini[0], "%d", p.n0[0]);
sprintf(nini[1], "%d", p.n0[1]);

strcat(p.filename, "N");
strcat(p.filename, gridN);
strcat(p.filename, "ka");
strcat(p.filename, factorD);
strcat(p.filename, "D");
strcat(p.filename, "kd");
strcat(p.filename, kd);
strcat(p.filename, "nA0");
strcat(p.filename, nini[0]);
strcat(p.filename, "nB0");
strcat(p.filename, nini[1]);
//strcat(p.filename, ".dat");


//printf("file name %s\n",p.filename);
strcpy(filename,p.filename);

}


	



