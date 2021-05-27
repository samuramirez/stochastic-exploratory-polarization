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

if(p.ktype<0 ||p.ktype > 6){printf("wrong ktype\n"); exit(1);}

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
		if(n==0){kmeso=p.kmi[reaction];}
		else{
		//Ours hetero
		bsr = sqrt(p.h*p.h/pi/n);
		lambda =  p.sigma/bsr ;
		kmeso = p.sigma < bsr ? p.kmi[reaction]/(1 + p.kmi[reaction]/(2*pi*Deff)*( log(1/lambda)/pow(1-lambda*lambda,2) - (3-lambda*lambda)/(1-lambda*lambda)/4 ) ) : p.kmi[reaction];	
		}
		break;
	case 5:
      //hybrid
		if(n<2){
			//Hellander-Petzold
			kmeso = p.kmi[reaction]/(1 + p.kmi[reaction]/(2*pi*Deff)*( log(p.h/sqrt(pi)/p.sigma) - (3 + 0.1951*2*pi)/4 ) );
        }
		else{
			bsr = sqrt(p.h*p.h/pi/n);
			lambda =  p.sigma/bsr ;
			kmeso = p.sigma < bsr ? p.kmi[reaction]/(1 + p.kmi[reaction]/(2*pi*Deff)*( log(1/lambda)/pow(1-lambda*lambda,2) - (3-lambda*lambda)/(1-lambda*lambda)/4 ) ) : p.kmi[reaction];
		}
		break;
	case 6:
      //hybrid
		if(n<1){
			//Hellander-Petzold
			kmeso = p.kmi[reaction]/(1 + p.kmi[reaction]/(2*pi*Deff)*( log(p.h/sqrt(pi)/p.sigma) - (3 + 0.1951*2*pi)/4 ) );
        }
		else{
			bsr = sqrt(p.h*p.h/pi/n);
			lambda =  p.sigma/bsr ;
			kmeso = p.sigma < bsr ? p.kmi[reaction]/(1 + p.kmi[reaction]/(2*pi*Deff)*( log(1/lambda)/pow(1-lambda*lambda,2) - (3-lambda*lambda)/(1-lambda*lambda)/4 ) ) : p.kmi[reaction];
		}
		break;

}


if (p.kmi[reaction]==0.0){kmeso=0;}
if (kmeso<0){printf("kmeso = %f, lambda = %f\n", kmeso, lambda); exit(1);}
if (p.ktype ==4 && kmeso > p.kmi[reaction] ){printf("kmeso = %f, kmi = %f , lambda = %f, reaction %d\n", kmeso, p.kmi[reaction], lambda, reaction); exit(1);}
return kmeso/(p.h*p.h);

}


void get_rate_timeq(int subv, int * positionq, double * timeq, int ** state, double ** rate, double t, double ** a, struct Parameters p){

int j;

//get rates for subvolume subv

//0 Cdc42T
//1 BemGEF42
//2 BemGEF
//3 Cdc42D		 
//4 Cdc42Dc
//5 BemGEFc

//Get propensities for each reaction
//kmeso has dimensions kmicro/h^2 --> 1/s 

//0 k1a_uM = 10;       %s-1,       BemGEFc -> BemGEF
a[subv][0] = p.kmi[0]*state[subv][5];	
	
//1 k1b_uM = 10;       %s-1,       BemGEF -> BemGEFc
a[subv][1] = p.kmi[1]*state[subv][2];	
	
//2 k2a_uM = 0.16;     %uM-1.s-1,  BemGEF + Cdc42D -> Cdc42T + BemGEF
a[subv][2] = kMesoHetero(2,2,3,subv,state,p)*state[subv][2]*state[subv][3];

//3 k2b_uM = 1.0;    %s-1,       Cdc42T -> Cdc42D 
a[subv][3] = p.kmi[3]*state[subv][0];
	
//4 k3_uM  = 0.35;     %uM-1.s-1,  BemGEF42 + Cdc42D -> Cdc42T + BemGEF42
a[subv][4] = kMesoHetero(4,1,3,subv,state,p)*state[subv][1]*state[subv][3];	
	
//5 k4a_uM = 10;       %uM-1.s-1,  BemGEF + Cdc42T -> BemGEF42
a[subv][5] = kMesoHetero(5,2,0,subv,state,p)*state[subv][2]*state[subv][0];		
	
//6 k4b_uM = 10;       %s-1,       BemGEF42 -> BemGEF + Cdc42T
//REVERSIBLE REACTION kmicro/dmicro=kmeso*h*h/dmeso  
//dmeso = dmicro*kmeso/kmicro
a[subv][6] = p.kmi[5] > 0 ? p.kmi[6]*kMesoHetero(5,2,0,subv,state,p)*p.h*p.h/p.kmi[5]*state[subv][1] : p.kmi[6]*state[subv][1];

	
//7 k5a_uM = 36;		%s-1		Cdc42c -> Cdc42
a[subv][7] = p.kmi[7]*state[subv][4];
	
//8 k5b_uM = 0.65;		%s-1		Cdc42 -> Cdc42c
a[subv][8] = p.kmi[8]*state[subv][3];

//9 k7_uM  = 10;       %uM-1.s-1,  BemGEFc + Cdc42T -> BemGEF42
a[subv][9] = kMesoHetero(9,5,0,subv,state,p)*state[subv][0]*state[subv][5];		

//10 k8_uM  = ;       %uM-1.s-1,  Cdc42c + BemGEF42 -> Cdc42T + BemGEF42 
a[subv][10] = kMesoHetero(10,4,1,subv,state,p)*state[subv][4]*state[subv][1];		

//Get total rate as sum of propensities 
rate[subv][0]=0;
for(j=0; j< p.Nreact ; j++){rate[subv][0]+=a[subv][j];}

//sum of diffusion rates. 
rate[subv][1]=0;
for(j=0; j<p.Nsp; j++){rate[subv][1]+=2*dim*p.D[j]/(p.h*p.h)*state[subv][j];}

//total rate (reaction+diffusion)
rate[subv][2] =  rate[subv][0]+rate[subv][1];

//if(rate[subv][2]==0.0){
//	printf("Total reaction rate = 0.0\n");
//	printf("propensity react %f\n",rate[subv][0]);
//	printf("total diffuse  %f\n",rate[subv][1]);
//	printf("difrate + reacrate  %f\n",rate[subv][2]);
//	exit(1);
//}

double randnumtime =  (double)rand()/RAND_MAX ;

//-log(1.0)/0.0 == nan

while( randnumtime == 1.0){randnumtime =  (double)rand()/RAND_MAX;}

timeq[positionq[subv]]= t + -log(randnumtime)/(rate[subv][2]);

}

void update_states_reaction(int reaction, int ** state, int subvnext){

//0 Cdc42T
//1 BemGEF42
//2 BemGEF
//3 Cdc42D		 
//4 Cdc42Dc
//5 BemGEFc

	//Update according to reaction
switch(reaction){
	case 0:
	//BemGEFc -> BemGEF
		state[subvnext][5]--;
		state[subvnext][2]++;
		break;
	case 1:
	//BemGEF -> BemGEFc	
		state[subvnext][2]--;
		state[subvnext][5]++;
		break;
	case 2:
	//BemGEF + Cdc42D -> Cdc42T + BemGEF	
		state[subvnext][3]--;
		state[subvnext][0]++;	
		break;
	case 3:
	//Cdc42T -> Cdc42D	
		state[subvnext][0]--;
		state[subvnext][3]++;	
		break;
	case 4:
	//BemGEF42 + Cdc42D -> Cdc42T + BemGEF42	
		state[subvnext][3]--;
		state[subvnext][0]++;	
		break;
	case 5:
	//BemGEF + Cdc42T -> BemGEF42	
		state[subvnext][2]--;
		state[subvnext][0]--;
		state[subvnext][1]++;	
		break;
	case 6:
	//BemGEF42 -> BemGEF + Cdc42T	
		state[subvnext][2]++;
		state[subvnext][0]++;
		state[subvnext][1]--;	
		break;
	case 7:
	//Cdc42Dc -> Cdc42D	
		state[subvnext][4]--;
		state[subvnext][3]++;
		break;
	case 8:
	//Cdc42D -> Cdc42Dc	
		state[subvnext][3]--;
		state[subvnext][4]++;
		break;	
	case 9:
	//BemGEFc + Cdc42T -> BemGEF42	
		state[subvnext][5]--;
		state[subvnext][0]--;
		state[subvnext][1]++;	
		break;	
	case 10:
	//BemGEF42 + Cdc42c -> Cdc42T + BemGEF42
		state[subvnext][4]--;
		state[subvnext][0]++;
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
strcpy(p->filename,"basegory");

	
p->N = atoi(argv[1]); //Number of voxels on the side
p->Nsubv = pow(p->N,dim); //number of voxels
p->L = 8; //box length um
p->h = p->L/p->N; //grid size

p->sigma = 0.02; //um

p->ktype = atoi(argv[3]); //0=ka, 1=kh, 2=kyj, 3=ksrhomo, 4=ksrhet

//number of species
p->Nsp = 6; 
//Number of reactions
p->Nreact = 11;

p->Dm = 0.0045; //um2s-1
p->Dc = 10.0; //um2s-1

p->D=(double *)malloc(p->Nsp*sizeof(double));
p->D[0]=p->Dm;
p->D[1]=p->Dm;
p->D[2]=p->Dm;
p->D[3]=p->Dm;
p->D[4]=p->Dc;
p->D[5]=p->Dc;
	
//0 Cdc42T
//1 BemGEF42
//2 BemGEF
//3 Cdc42D		 
//4 Cdc42Dc
//5 BemGEFc

//Initial number of proteins
p->n0=(int *)malloc(p->Nsp*sizeof(int));
p->n0[0]= 0; 
p->n0[1]= 0;
p->n0[2]= 0;
p->n0[3]= 0;
p->n0[4]= 5000;
p->n0[5]= 700;
				
//Macroscopic reactionn rates
	//0 k1a_uM = 10;       %s-1,       BemGEFc -> BemGEF	
	//1 k1b_uM = 10;       %s-1,       BemGEF -> BemGEFc	
	//2 k2a_uM = 0.16;     %uM-1.s-1,  BemGEF + Cdc42D -> Cdc42T + BemGEF
	//3 k2b_uM = 1.0;    %s-1,       Cdc42T -> Cdc42D 	
	//4 k3_uM  = 0.35;     %uM-1.s-1,  BemGEF42 + Cdc42D -> Cdc42T + BemGEF42	
	//5 k4a_uM = 10;       %uM-1.s-1,  BemGEF + Cdc42T -> BemGEF42	
	//6 k4b_uM = 10;       %s-1,       BemGEF42 -> BemGEF + Cdc42T	
	//7 k5a_uM = 36;		%s-1		Cdc42c -> Cdc42	
	//8 k5b_uM = 0.65;		%s-1		Cdc42 -> Cdc42c
	//9 k7_uM  = 10;       %uM-1.s-1,  BemGEFc + Cdc42T -> BemGEF42
	//10 k8_uM  = ;       %uM-1.s-1,  Cdc42c + BemGEF42 -> Cdc42T + BemGEF42 
		
p->kmi=(double *)malloc(p->Nreact*sizeof(double));

p->kmi[0]=0.1; //k1a
p->kmi[1]=10; //k1b
p->kmi[2]=0.032; //k2a  Dtot=2*0.0025 , kmi/Dtot=8.3
p->kmi[3]=0.63; //k2b
p->kmi[4]=0.07; //k3 Dtot=2*0.0025, kmi/Dtot=282
p->kmi[5]=2.0;  //k4a Dtot=2*0.0025, kmi/Dtot=15
p->kmi[6]=10; //k4b
p->kmi[7]=4; //k5a
p->kmi[8]=6.5; //k5b
p->kmi[9]=0.2; //k7  Dtot=10+0.0025,  kmi/Dtot=0.13
p->kmi[10]=0; // Dtot=10+0.0025,  kmi/Dtot=0.13

						
}

void make_filename(struct Parameters p, char filename[]){

//Making filename
char  gridN[10];
sprintf(gridN, "%d", p.N);

if(p.ktype<0 ||p.ktype > 6){printf("wrong ktype\n"); exit(1);}
if(p.ktype==0){strcat(p.filename, "ka");}
if(p.ktype==1){strcat(p.filename, "kh");}
if(p.ktype==2){strcat(p.filename, "kyj");}
if(p.ktype==3){strcat(p.filename, "ksr");}
if(p.ktype==4){strcat(p.filename, "ksrhet");}
if(p.ktype==5){strcat(p.filename, "khyb");}
if(p.ktype==6){strcat(p.filename, "khyb0");}


//char  nini[p.Nsp][10];
//sprintf(nini[0], "%d", p.n0[0]);
//sprintf(nini[1], "%d", p.n0[1]);

strcat(p.filename, "N");
strcat(p.filename, gridN);

//strcat(p.filename, ".dat");
//printf("file name %s\n",p.filename);
strcpy(filename,p.filename);

}

void get_data(int ** state, struct Parameters p){
	
FILE * fileConc;
int i,j;

fileConc=fopen("finalconc.dat","w");

for(i=0; i<p.Nsubv; i++){
	
	if(i%p.N==0 && i != 0){fprintf(fileConc,"\n");}
	
	fprintf(fileConc,"%d\t",state[i][0]);
}

fclose(fileConc);	
	
	
}
	
	



