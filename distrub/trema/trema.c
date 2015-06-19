#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <time.h>


int N = 100000, N1 = 100, d=3, n=10000;
double R = 1.0, D=2.7, alpha;
double **FOO;

//#include "cube_2d.c"
#include "sphere.c"

int 
select_trema( /*char* inp_name,*/  char* out_name ) { 
	double x, y, z;
	FILE *out /*, *inp*/;
	double *point;
	int i, count;
	double r;


	if( (out = fopen(out_name, "w")) == NULL ) {
		printf("Can't create file %s \n", out_name);
		return 1; 
	}
	
	count = 0;
	void_centers();
	point = (double *) calloc( d, sizeof( double ) );
	while ( count < n ) {
		//count++;
		r = 0;
		for(i=0; i<d; i++) {
			point[i] = ((double) rand() / RAND_MAX) ;	
			r += point[i]*point[i]; 
		}
		if(r>1) continue; 
		if(is_it_good(point)) { 
			for(i=0; i<d; i++) {
				fprintf(out, " %lf ", point[i] * R);
			}	
			fprintf(out, "\n");
			count++;
			printf("\r%d", count);
		}
	}

	printf("\n");
	fclose(out);
	/*
	out = fopen("pewpew.dat", "w");
	double C = d-D;

	for(i=0; i<N; i++) {
		double S = C/(2.0*(i+1));
		fprintf(out, "%lf %lf %lf\n", FOO[i][0], FOO[i][1], S/M_PI);
	}
	fclose(out);
	*/
	//fclose(inp);
	return N;
}

int main(int argc, char** argv) {
	char* out_name = "trema.dat";
	char optString[] = {"s:p:N:n:R:D:d:r:"};
    double r = -1;
	int opt; 
	opt = getopt( argc, argv, optString );
	while( opt != -1 ) {
        switch( opt ) {
		case 's':
	                out_name = optarg; 
        	        break;
		case 'n':
			if(!sscanf(optarg, "%d", &n)) {
				printf("n should be int\n");
			}
                	break;	
		case 'N':
			if(!sscanf(optarg, "%d", &N)) {
				printf("N should be int\n");
			}
	                break;	
		case 'D':
			if(!sscanf(optarg, "%lf", &D)) {
				printf("D should be real\n");
			}
			break;
		case 'r':
			if(!sscanf(optarg, "%lf", &r)) {
				printf("r is wrong\n");
			}
			break;

		case 'd':
			if(!sscanf(optarg, "%d", &d)) {
				printf("d should be integer\n");
			}
                	break;
		case 'R':
			if(!sscanf(optarg, "%lf", &R)) {
				printf("R should be real\n");
			}
                	break;
            default:
			printf("What?");
			break;
        }
        opt = getopt( argc, argv, optString );
	}
	if(d<D) {
		printf("d<D, I can't do that.\n");
	}
	switch( d ) {
		case 1:
			alpha = 1; // Srsly o_O
			break;
		case 2: 
			alpha = M_PI; // S = \pi r^2	
			break;
		case 3: 
			alpha = 4.0/3.0*M_PI;
			break; 
		case 4: 
			alpha = M_PI*M_PI/2; // lol?
			break; 
		default:
			printf("Sorry, what exactly dimention?\n");
			return 42;
	}
	srand(time(0));
    N1 = limit_N(0.3);
    if( r != -1 ) { 
        N = limit_N( r ); 
    } else {
        r = limit_r( N );
    }
	printf("N = %d, N1 = %d,  d=%d, D=%e, r = %e\n", N, N1, d, D, r );
	select_trema(/*inp_name,*/ out_name);
	return 0;
}

