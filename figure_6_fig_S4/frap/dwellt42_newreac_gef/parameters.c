#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "header.h"
#include <string.h>

void initialize_parameters(char * argv[], int argc,struct Parameters * p){
	//Future implementation should read parameters and system from an input text file
strcpy(p->filename,"dynfrap");

	
p->N = atoi(argv[2]); //Number of voxels on the side
p->Nsubv = pow(p->N,dim); //number of voxels
p->L = atof(argv[3]); //box length um
p->h = p->L/p->N; //grid size

p->interval=atof(argv[4]);
p->time_samp=atof(argv[5]);

p->sigma = atof(argv[6]); //um

p->ktype = 4; //0=ka, 1=kh, 2=kyj, 3=ksrhomo, 4=ksrhet

//number of species
p->Nsp = 10; 
//Number of reactions
p->Nreact = 24;

p->Dm = atof(argv[7]); //um2s-1
p->Dc = atof(argv[8]); //um2s-1

p->D=(double *)malloc(p->Nsp*sizeof(double));
p->D[0]=p->Dm;
p->D[1]=p->Dm;
p->D[2]=p->Dm;
p->D[3]=p->Dm;
p->D[4]=p->Dc;
p->D[5]=p->Dc;
p->D[6]=p->Dm;
p->D[7]=p->Dm;
p->D[8]=p->Dm;
p->D[9]=p->Dc;
	
//0 Cdc42T
//1 BemGEF42
//2 BemGEF
//3 Cdc42D		 
//4 Cdc42Dc
//5 BemGEFc
//6 BemGEF42d
//7 Cdc42Td
//8 Cdc42Dd
//9 Cdc42Ddc

//Initial number of proteins
p->n0=(int *)malloc(p->Nsp*sizeof(int));
p->n0[0]= 0; 
p->n0[1]= 0;
p->n0[2]= 0;
p->n0[3]= 0;
p->n0[4]= atoi(argv[9]);
p->n0[5]= atoi(argv[10]);
p->n0[6]= 0; 
p->n0[7]= 0;
p->n0[8]= 0;
p->n0[9]= 0;
				
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
//There are indeed only 11 kmis, more kmis are declared here but are not used
//for association reactions kmi = mike's lambda * pi*sigma^2
p->kmi[0]=atof(argv[11]);
p->kmi[1]=atof(argv[12]);
p->kmi[2]=atof(argv[13]);
p->kmi[3]=atof(argv[14]);
p->kmi[4]=atof(argv[15]); 
p->kmi[5]=atof(argv[16]); 
p->kmi[6]=atof(argv[17]);		
p->kmi[7]=atof(argv[18]);
p->kmi[8]=atof(argv[19]);
p->kmi[9]=atof(argv[20]); 
p->kmi[10]=atof(argv[21]);

	FILE * fileparams;
	fileparams=fopen("params.dat","w");
	int i;
	for (i=2 ; i<argc; i++){
		fprintf(fileparams,"%s\n",argv[i]);
	}
	fclose(fileparams);
}

void initialize(int ** state, int * positionq, int * subvolumeq, double * timeq, double ** rate, double * t, double ** a, struct Parameters p)
{

int subvi,i,j,temp;
int	totstate[p.Nsp];
*t=0.0;

//Set species randomly on the grid
for(j=0;j<p.Nsp;j++){
	for(i=0; i<p.Nsubv ; i++ ){state[i][j]=0;}
	for(i=0; i < p.n0[j] ; i++ ){subvi=rand()%p.Nsubv; state[subvi][j]++;}
}


//Set cytosolic species in the middle of the grid
//for(j=0;j<p.Nsp;j++){
//   for(i=0; i<p.Nsubv ; i++ ){state[i][j]=0;}
//         for(i=0; i < p.n0[j] ; i++ ){subvi = p.Nsubv/2 + p.N/2; state[subvi][j]++;}
//}
//set some Cdc42T in the middle of the grid

subvi = p.Nsubv/2 + floor(p.N/2);
state[subvi][0]+=500;


/*
//read initial concentrations from file
FILE * fileConc;
fileConc=fopen("finalconcHSS.dat","r");
if (fileConc==NULL) { printf("no such file."); } 

for(j=0;j<p.Nsp; j++){
	for(i=0; i<p.Nsubv; i++){fscanf(fileConc,"%d",&state[i][j]);}
	//printf("\n");	
}
fclose(fileConc);	

//Check total amount of proteins
for(j=0;j<p.Nsp;j++){
	totstate[j]=0;
	for(i=0; i<p.Nsubv; i++){
       totstate[j] += state[i][j];
	}
}
printf("tot 42  = %d\n",totstate[0]+totstate[1]+totstate[3]+totstate[4]);
printf("tot GEF  = %d\n",totstate[1]+totstate[2]+totstate[5]);
*/
//INITIALIZE TRIVIALLY subvolumeq[positionq], timeq, positioniq[subvolume]
for(i=0;i<p.Nsubv;i++){positionq[i]=i; subvolumeq[positionq[i]]=i; }

//CALCULATE SUM OF REACTION AND DIFFUSION RATES AND TIME FOR NEXT REACTION
for(i=0;i<p.Nsubv;i++){get_rate_timeq(i,positionq, timeq, state,rate,*t,a,p);}

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



void get_rate_timeq(int subv, int * positionq, double * timeq, int ** state, double ** rate, double t, double ** a, struct Parameters p){

int j;

//get rates for subvolume subv

//0 Cdc42T
//1 BemGEF42
//2 BemGEF
//3 Cdc42D		 
//4 Cdc42Dc
//5 BemGEFc
//6 BemGEF42d
//7 Cdc42Td
//8 Cdc42Dd
//9 Cdc42Ddc

//Get propensities for each reaction
//kmeso has dimensions kmicro/h^2 --> 1/s 

//0 k1a_uM = 10;       %s-1,       BemGEFc -> BemGEF
a[subv][0] = p.kmi[0]*state[subv][5];	

//1 k1b_uM = 10;       %s-1,       BemGEF -> BemGEFc
a[subv][1] = p.kmi[1]*state[subv][2];	
	
//2 k2a_uM = 0.16;     %uM-1.s-1,  BemGEF + Cdc42D -> Cdc42T + BemGEF
a[subv][2] = kMesoHetero(2,2,3,subv,state,p)*state[subv][2]*state[subv][3];
//2 k2a_uM = 0.16;     %uM-1.s-1,  BemGEF + Cdc42Dd -> Cdc42Td + BemGEF
a[subv][11] = kMesoHetero(2,2,8,subv,state,p)*state[subv][2]*state[subv][8];

//3 k2b_uM = 1.0;    %s-1,       Cdc42T -> Cdc42D 
a[subv][3] = p.kmi[3]*state[subv][0];
//3 k2b_uM = 1.0;    %s-1,       Cdc42Td -> Cdc42Dd 
a[subv][12] = p.kmi[3]*state[subv][7];
	
//4 k3_uM  = 0.35;     %uM-1.s-1,  BemGEF42 + Cdc42D -> Cdc42T + BemGEF42
a[subv][4] = kMesoHetero(4,1,3,subv,state,p)*state[subv][1]*state[subv][3];	
//4 k3_uM  = 0.35;     %uM-1.s-1,  BemGEF42d + Cdc42D -> Cdc42T + BemGEF42d
a[subv][13] = kMesoHetero(4,6,3,subv,state,p)*state[subv][6]*state[subv][3];	
//4 k3_uM  = 0.35;     %uM-1.s-1,  BemGEF42 + Cdc42Dd -> Cdc42Td + BemGEF42
a[subv][14] = kMesoHetero(4,1,8,subv,state,p)*state[subv][1]*state[subv][8];	
//4 k3_uM  = 0.35;     %uM-1.s-1,  BemGEF42d + Cdc42Dd -> Cdc42Td + BemGEF42d
a[subv][15] = kMesoHetero(4,6,8,subv,state,p)*state[subv][6]*state[subv][8];	
	
//5 k4a_uM = 10;       %uM-1.s-1,  BemGEF + Cdc42T -> BemGEF42
a[subv][5] = kMesoHetero(5,2,0,subv,state,p)*state[subv][2]*state[subv][0];		
//5 k4a_uM = 10;       %uM-1.s-1,  BemGEF + Cdc42Td -> BemGEF42d
a[subv][16] = kMesoHetero(5,2,7,subv,state,p)*state[subv][2]*state[subv][7];		
	
//6 k4b_uM = 10;       %s-1,       BemGEF42 -> BemGEF + Cdc42T
//REVERSIBLE REACTION kmicro/dmicro=kmeso*h*h/dmeso  
//dmeso = dmicro*kmeso/kmicro
a[subv][6] = p.kmi[5] > 0 ? p.kmi[6]*kMesoHetero(5,2,0,subv,state,p)*p.h*p.h/p.kmi[5]*state[subv][1] : p.kmi[6]*state[subv][1];
//6 k4b_uM = 10;       %s-1,       BemGEF42d -> BemGEF + Cdc42Td
a[subv][17] = p.kmi[5] > 0 ? p.kmi[6]*kMesoHetero(5,2,7,subv,state,p)*p.h*p.h/p.kmi[5]*state[subv][6] : p.kmi[6]*state[subv][6];

	
//7 k5a_uM = 36;		%s-1		Cdc42c -> Cdc42
a[subv][7] = p.kmi[7]*state[subv][4];
//7 k5a_uM = 36;		%s-1		Cdc42Ddc -> Cdc42Dd
a[subv][18] = p.kmi[7]*state[subv][9];
	
//8 k5b_uM = 0.65;		%s-1		Cdc42 -> Cdc42c
a[subv][8] = p.kmi[8]*state[subv][3];
//8 k5b_uM = 0.65;		%s-1		Cdc42Dd -> Cdc42Ddc
a[subv][19] = p.kmi[8]*state[subv][8];

//9 k7_uM  = 10;       %uM-1.s-1,  BemGEFc + Cdc42T -> BemGEF42
a[subv][9] = kMesoHetero(9,5,0,subv,state,p)*state[subv][0]*state[subv][5];		
//9 k7_uM  = 10;       %uM-1.s-1,  BemGEFc + Cdc42Td -> BemGEF42d
a[subv][20] = kMesoHetero(9,5,7,subv,state,p)*state[subv][5]*state[subv][7];		

//10 k8_uM  = ;       %uM-1.s-1,  Cdc42c + BemGEF42 -> Cdc42T + BemGEF42 
a[subv][10] = kMesoHetero(10,4,1,subv,state,p)*state[subv][4]*state[subv][1];		
//10 k8_uM  = ;       %uM-1.s-1,  Cdc42Ddc + BemGEF42 -> Cdc42Td + BemGEF42 
a[subv][21] = kMesoHetero(10,9,1,subv,state,p)*state[subv][9]*state[subv][1];		
//10 k8_uM  = ;       %uM-1.s-1,  Cdc42c + BemGEF42d -> Cdc42T + BemGEF42d 
a[subv][22] = kMesoHetero(10,4,6,subv,state,p)*state[subv][4]*state[subv][6];		
//10 k8_uM  = ;       %uM-1.s-1,  Cdc42Ddc + BemGEF42d -> Cdc42Td + BemGEF42d 
a[subv][23] = kMesoHetero(10,9,6,subv,state,p)*state[subv][9]*state[subv][6];		

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
//6 BemGEFd42
//7 Cdc42Td
//8 Cdc42Dd
//9 Cdc42Ddc
//
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
	
	case 11:
	//BemGEF + Cdc42Dd -> Cdc42Td + BemGEF	
		state[subvnext][8]--;
		state[subvnext][7]++;	
		break;

	case 3:
	//Cdc42T -> Cdc42D	
		state[subvnext][0]--;
		state[subvnext][3]++;	
		break;
	case 12:
	//Cdc42Td -> Cdc42Dd	
		state[subvnext][7]--;
		state[subvnext][8]++;	
		break;

	case 4:
	//BemGEF42 + Cdc42D -> Cdc42T + BemGEF42	
		state[subvnext][3]--;
		state[subvnext][0]++;	
		break;
	case 13:
	//BemGEF42d + Cdc42D -> Cdc42T + BemGEF42d	
		state[subvnext][3]--;
		state[subvnext][0]++;	
		break;
	case 14:
	//BemGEF42 + Cdc42Dd -> Cdc42Td + BemGEF42	
		state[subvnext][8]--;
		state[subvnext][7]++;	
		break;
	case 15:
	//BemGEF42d + Cdc42Dd -> Cdc42Td + BemGEF42d	
		state[subvnext][8]--;
		state[subvnext][7]++;	
		break;

	case 5:
	//BemGEF + Cdc42T -> BemGEF42	
		state[subvnext][2]--;
		state[subvnext][0]--;
		state[subvnext][1]++;	
		break;
	case 16:
	//BemGEF + Cdc42Td -> BemGEF42d	
		state[subvnext][2]--;
		state[subvnext][7]--;
		state[subvnext][6]++;	
		break;

	case 6:
	//BemGEF42 -> BemGEF + Cdc42T	
		state[subvnext][2]++;
		state[subvnext][0]++;
		state[subvnext][1]--;	
		break;
	case 17:
	//BemGEF42d -> BemGEF + Cdc42Td	
		state[subvnext][2]++;
		state[subvnext][7]++;
		state[subvnext][6]--;	
		break;

	case 7:
	//Cdc42Dc -> Cdc42D	
		state[subvnext][4]--;
		state[subvnext][3]++;
		break;
	case 18:
	//Cdc42Ddc -> Cdc42Dd should not occur	
		state[subvnext][9]--;
		state[subvnext][8]++;
		break;

	case 8:
	//Cdc42D -> Cdc42Dc	
		state[subvnext][3]--;
		state[subvnext][4]++;
		break;	
	case 19:
	//Cdc42Dd -> Cdc42Ddc	
		state[subvnext][8]--;
		state[subvnext][4]++; //this is the only change with respect to frap
		break;	

	case 9:
	//BemGEFc + Cdc42T -> BemGEF42	
		state[subvnext][5]--;
		state[subvnext][0]--;
		state[subvnext][1]++;	
		break;
	case 20:
	//BemGEFc + Cdc42Td -> BemGEF42d	
		state[subvnext][5]--;
		state[subvnext][7]--;
		state[subvnext][6]++;	
		break;

	case 10:
	//BemGEF42 + Cdc42c -> Cdc42T + BemGEF42
		state[subvnext][4]--;
		state[subvnext][0]++;
		break;	
	case 21:
	//BemGEF42 + Cdc42Ddc -> Cdc42Td + BemGEF42
		state[subvnext][9]--;
		state[subvnext][7]++;
		break;	
	case 22:
	//BemGEF42d + Cdc42c -> Cdc42T + BemGEF42d
		state[subvnext][4]--;
		state[subvnext][0]++;
		break;	
	case 23:
	//BemGEF42d + Cdc42Ddc -> Cdc42Td + BemGEF42d
		state[subvnext][9]--;
		state[subvnext][7]++;
		break;	

	default:
		printf("invalid reaction type");
		exit(1);
}//end update according to  reaction
	
		
}

void make_filename(struct Parameters p, char filename[], char * argv[] ){

//Making filename
char  gridN[10];
sprintf(gridN, "%d", p.N);

//char param1[10];
//sprintf(param1,"%.2f",p.kmi[0]);

/*
if(p.ktype<0 ||p.ktype > 4){printf("wrong ktype\n"); exit(1);}
if(p.ktype==0){strcat(p.filename, "ka");}
if(p.ktype==1){strcat(p.filename, "kh");}
if(p.ktype==2){strcat(p.filename, "kyj");}
if(p.ktype==3){strcat(p.filename, "ksr");}
if(p.ktype==4){strcat(p.filename, "ksrhet");}
*/

strcat(p.filename, "N");
strcat(p.filename, gridN);

strcat(p.filename, argv[22]);
strcat(p.filename, argv[23] );
strcat(p.filename, "sim" );
strcat(p.filename, argv[1] );
strcat(p.filename, ".dat");

//printf("file name %s\n",p.filename);
strcpy(filename,p.filename);

}


	
	



