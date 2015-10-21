#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <time.h>


int N = 100000, N1 = 10, d=3, n=10000;
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
	
    FILE * out2; 
    out2 = fopen("mad.dat", "w");
	count = 0;
	void_centers();
	point = (double *) calloc( d, sizeof( double ) );
	while ( count < n ) {
		//count++;
		r = 0;
		for(i=0; i<d; i++) {
			point[i] = (((double) rand() / RAND_MAX)) ;	
			r += point[i]*point[i]; 
		}
		if(r>1) continue; 
		if(is_it_good(point)) { 
			for(i=0; i<d; i++) {
				fprintf(out, " %3.10lf ", point[i] * R);
			}	
			fprintf(out, "\n");
			count++;
			printf("\r%d%%", (int) (1.*count/n*100.));
		}
		for(i=0; i<d; i++) {
			fprintf(out2, " %lf ", point[i] * R);
		}	
		fprintf(out2, "\n");
        
	}

	printf("\n");
	fclose(out);
	fclose(out2);

	out = fopen("pewpew.dat", "w");
	double C = d-D;

	for(i=0; i<N; i++) {
		double S = C/(2.0*(i+1));
        if( sqrt(FOO[i][0] * FOO[i][0] +
                 FOO[i][1] * FOO[i][1] +
                 FOO[i][2] * FOO[i][2] ) > 1 )
            continue;
		fprintf(out, "%lf %lf %lf %lf\n", FOO[i][0]*R, FOO[i][1]*R, FOO[i][2]*R, FOO[i][3]);
	}
	fclose(out);
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
    // N1 = limit_N(0.3);
    if( r != -1 ) { 
        N = limit_N( r ); 
    } else {
        r = limit_r( N );
    }
	printf("voids N = %d\t points n = %d\nskip N1 = %d\tspace  d=%d\tdim D=%1.3e\t???? r = %1.3e\talpha = %1.3e\n", N, n, N1, d, D, r, alpha );
	select_trema(/*inp_name,*/ out_name);
	return 0;
}

