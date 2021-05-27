#include <math.h>
#define dim 2
#define pi acos(-1.0)

struct Parameters {
	
    double interval;
	double time_samp;
	
	double * kmi;
	
	double Dm;
	double Dc;
			
    double * D;
	int * n0;
	
    double sigma;
    double L;
    int N;
    int Nsubv;
    double h;
    int Nsp;
    int Nreact;
	
	int ktype;
	
	char filename[100];
    
};

void allocate(int ** state, double ** rate, double ** a, double * timeq, int * positionq, int * subvolumeq, int ** connect, struct Parameters p);
void swap(int position1, int position2,int *  positionq,double *  timeq, int * subvolumeq);
void heapify(int i, int * positionq, double * timeq, int * subvolumeq, struct Parameters p);
void build_queue(int *positionq,double * timeq,int * subvolumeq, struct Parameters p);
void update_states_reaction(int reaction, int ** state, int subvnext);
void get_rate_timeq(int i, int * positionq , double * timeq, int ** state , double ** rate, double t, double ** a, struct Parameters p);
void create_connect(int ** connect, struct Parameters p);
void initialize(int ** state, int * positionq, int * subvolumeq, double * timeq, double ** rate, double * t, double ** a, struct Parameters p);
void sort (int subv, int * positionq, double * timeq, int * subvolumeq, struct Parameters p );
void sort_no_heapify (int subv, int * positionq, double * timeq, int * subvolumeq, struct Parameters p);
void checkqueue(int * positionq, double * timeq, int * subvolumeq, struct Parameters p);
void react(double randnum, int subvnext, int ** state, double ** rate, double ** a, double * timeq, int * positionq, int * subvolumeq, double t, struct Parameters p);
void diffuse(double randnum, int subvnext, int ** state, double ** rate, double ** a, double * timeq, int * positionq, int * subvolumeq, double t, int ** connect, struct Parameters p);
void initialize_parameters(char * argv[], int argc,struct Parameters * p);
void make_filename(struct Parameters p, char filename[], char * argv[]);
double kMesoHetero(int reaction,int species1,int species2, int subv,int ** state, struct Parameters p);
void get_data(int ** state, struct Parameters p);
