#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#define MAX(a, b)       ((a) > (b) ? (a) : (b))
#define MIN(a, b)       ((a) < (b) ? (a) : (b))

double r2d = 180/M_PI, d2r = M_PI / 180;
int SELECT_COORDINATES = 0; // CARTESIAN
char   *input_suff, *output_suff; 
double A[2] = {0,0}, B[2] = {0, 0}, R[2] = {0, 300};
double q =1;
int N;
double n0, solid_angle;

void
spherical_coordinats(int n, double *x, double *y, double *z, double *a, double *b, double *r) {
	int i; 
	for( i=0; i<n; i++ ) {
		r[i] = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
        a[i] = atan2(y[i], x[i]) ;  
		if( a[i] < 0 ) 
			a[i] += 2*M_PI; 
		if( x[i] == 0 & y[i] == 0 ) {
			if( z[i] > 0 ) b[i] = +M_PI/2.;
			if( z[i] < 0 ) b[i] = -M_PI/2.;
		}
        b[i] = asin(z[i]/r[i]) ; 
	}
}	

void
cartesian_coordinats(int n, double *x, double *y, double *z, double *a, double *b, double *r) {
	int i; 
	for( i=0; i<n; i++ ) {
		x[i] = r[i] * cos(b[i]) * cos(a[i]); 
		y[i] = r[i] * cos(b[i]) * sin(a[i]); 
		z[i] = r[i] * sin(b[i]); 
	}
}	

int
inside_sp(double a,double b,double r) {
	if(r < R[0]) { return 1; }
	if(r > R[1]) { return 2; }
    if( A[0] != A[1] ) {
        if(a < A[0]) { return 3; }
        if(a > A[1]) { return 4; }
    }      
    if( B[0] != B[1] ) { 
	    if(b < B[0]) { return 5; }
	    if(b > B[1]) { return 6; }
    }
	return 0;
}


void
read_data() {
	FILE 	*finput, *foutput; 
	char	 input_name[80];
	char	output_name[80];
	double 	x, y, z, a, b, r;
	int i, oor = 0, hack; // out of range;	
	int NMAX = N;
	double alim[] = {360*d2r, 0}, blim[] = {90*d2r, -90*d2r}, rlim[] = {1e10, 0};
	
	sprintf(input_name, "%s.dat", input_suff);
	finput = fopen(input_name, "r");
	if( finput == NULL ) {
		printf("Can't find file %s.\n", input_name);
		exit(10);
	}

	sprintf(output_name, "%s", output_suff);
	foutput = fopen(output_name, "w");
	if( foutput == NULL ) {
		printf("Can't open file %s.\n", output_name);
		exit(10);
	}
    printf("%s -> %s\n", input_name, output_name);
	double toss;
    int nfields = 3;
    size_t lineno = 0;

    // Read a file containing a list of 3-D vectors as floats.
    N = 0;
    for (;fscanf(finput, "%lf %lf %lf", &x, &y, &z) == 3; ) {
        if( SELECT_COORDINATES == 0 ) {
            spherical_coordinats(1, &(x), &(y), &(z), &(a), &(b), &(r));
        } else { 
            a = x*d2r; b = y*d2r; r = z;
            cartesian_coordinats(1, &(x), &(y), &(z), &(a), &(b), &(r));
        } 
		if( a < alim[0] ) alim[0] = a; 
		if( a > alim[1] ) alim[1] = a; 
		if( b < blim[0] ) blim[0] = b; 
		if( b > blim[1] ) blim[1] = b; 
		if( r < rlim[0] ) rlim[0] = r; 
		if( r > rlim[1] ) rlim[1] = r; 
		toss=((double) rand() / RAND_MAX);
		if( (inside_sp(a, b, r) == 0) && ( q >= toss))  {
			N++;
            fprintf(foutput, "%lf %lf %lf\n", x, y, z);
		}
    }

	fclose(finput);
	fclose(foutput);

    // Parameters of selected objects
    if( (A[1] == A[0]) & (B[0] == B[1]) ) {
        solid_angle = 4 * M_PI;
    } else {
        if( A[1] == A[0] ) 
            solid_angle = 2*M_PI *  (sin(B[1]) - sin(B[0]));
        else 
	        solid_angle = (A[1] - A[0]) * (sin(B[1]) - sin(B[0]));
    }
	n0 = (double) 3. * N / (solid_angle * (pow(R[1],3) - pow(R[0],3)) );
	printf("%d selected.\n", N);
	printf("A : %+3.3lf %+3.3lf -> %+3.3lf %+3.3lf \n", alim[0]*r2d, alim[1]*r2d, A[0]*r2d, A[1]*r2d );
	printf("B : %+3.3lf %+3.3lf -> %+3.3lf %+3.3lf \n", blim[0]*r2d, blim[1]*r2d, B[0]*r2d, B[1]*r2d );
	printf("R : %+3.3lf %+3.3lf -> %+3.3lf %+3.3lf \n", rlim[0]    , rlim[1]    , R[0]    , R[1]     );
	printf("Solid angle : %1.3le\nDensety : %1.3le\n", solid_angle, n0);
}



int main( int argc, char *argv[] ) {
	char optString[] = {"d:a:b:r:q:o:c"};
	int opt, got_out = 0; 
	opt = getopt( argc, argv, optString );
	while( opt != -1 ) {
		switch( opt ) {
			case 'd':
                input_suff = optarg; 
				break;
			case 'a':
				sscanf(optarg, "%lf", &A[0]); 
				sscanf(argv[optind], "%lf", &A[1]);
				optind++;
				break ;
			case 'b':
				sscanf(optarg, "%lf", &B[0]); 
				sscanf(argv[optind], "%lf", &B[1]);
				optind++;
				break ;
			case 'r':
				sscanf(optarg, "%lf", &(R[0])); 
				sscanf(argv[optind], "%lf", &(R[1]));
				optind++;
				break ;
			case 'q': 
				sscanf(optarg, "%lf", &q); 
				//sscanf(optarg, "%lf", %q);
				break;
			case 'o':
				got_out = 1;
				output_suff = optarg; 
				break;
            case 'c': 
                SELECT_COORDINATES = 1; // spherical
                printf("I'll work with spherical coordinates.\n");
                break;

			default:
				printf("What? %d\n", opt);
				break;
    	}
	    opt = getopt( argc, argv, optString );
	}

	A[0]=A[0]*M_PI/180;  
	A[1]=A[1]*M_PI/180;  
	B[0]=B[0]*M_PI/180;  
	B[1]=B[1]*M_PI/180;  
    
	if( got_out != 1 ) {
        char tmp[80];
        sprintf(tmp, "%s_sample.dat", input_suff);
        output_suff = tmp;
    }

    read_data();

	return(0);	
}

