#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <time.h>

int N = 10000;
double A[2] = {0, 360}, B[2] = {-90, 90}, R[2] = {0, 300};
const double r2d = 180/M_PI, d2r = M_PI / 180;
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
	if(a < A[0]) { return 3; }
	if(a > A[1]) { return 4; }
	if(b < B[0]) { return 5; }
	if(b > B[1]) { return 6; }
	return 0;
}

int
inside_ca(double x,double y, double z) {
	double a,b,r; 
	spherical_coordinats(1,&x,&y,&z,&a,&b,&r);
	return inside_sp(a,b,r);
}

int
get_poisson_from_square( char * output ) {
	double x, y, z;
    double a, b, r;
	FILE *out;
	int i;
	
	if( (out = fopen(output, "w")) == NULL ) {
		printf("Can't create file %s \n", output);
		return 1; 
	}
	i = 0;
	while ( i<N ) {
		x = ((double) rand() / RAND_MAX) * 2 * R[1] - R[1] ;
		y = ((double) rand() / RAND_MAX) * 2 * R[1] - R[1] ;
		z = ((double) rand() / RAND_MAX) * 2 * R[1] - R[1] ;
	    spherical_coordinats(1,&x,&y,&z,&a,&b,&r);
        if( inside_ca(x, y, z) == 0 )    {
		    fprintf(out, "%le %le %le\n", x, y, z);
		    //fprintf(out, "%le %le %le %le %le %le\n", x, y, z, a, b, r);
			i++;
		}
	}
	fclose(out);
	return N;
}

/*
int
get_poisson_sector2( char * output ) {
	double x, y, z;
	double ri, ra, dec; 
	FILE *out;
	int i;
	
	if( (out = fopen(output, "w")) == NULL ) {
		printf("Can't create file %s \n", output);
		return 1; 
	}
	
	r[0] = pow(r[0], 3);
	r[1] = pow(r[1], 3);

	Dec[0] = sin(Dec[0]);
	Dec[1] = sin(Dec[1]);

	for(i=0; i<N; i++) {
		ri  = pow(  ((double) rand() / RAND_MAX) * ( r  [1] - r  [0] ) + r  [0] , 1./3. );
		ra  =       ((double) rand() / RAND_MAX) * ( Ra [1] - Ra [0] ) + Ra [0] ; 
		dec = asin( ((double) rand() / RAND_MAX) * ( Dec[1] - Dec[0] ) + Dec[0]);
		
		x =  ri * cos(dec) * cos(ra);
		y =  ri * cos(dec) * sin(ra);
		z =  ri * sin(dec)          ;
		//fprintf(out, "%lf %lf %lf\n", ri, ra, dec);
		fprintf(out, "%le %le %le\n", x, y, z);
	}
	fclose(out);
	return N;
}

int
get_poisson_sphere2( char * output ) {
	double r, x, y, z;
	FILE *out;
	int i;
	
	if( (out = fopen(output, "w")) == NULL ) {
		printf("Can't create file %s \n", output);
		return 1; 
	}
	i = 0;
	while ( i<N ) {
		x = ((double) rand() / RAND_MAX) * 2 * R - R ;
		y = ((double) rand() / RAND_MAX) * 2 * R - R ;
		z = ((double) rand() / RAND_MAX) * 2 * R - R ;
		if( (r = sqrt(x*x + y*y + z*z) ) < R )    {
			fprintf(out, "%lf %lf %lf\n", x, y, z);
			i++;
		}
	}
	fclose(out);
	return N;
}
*/

int main(int argc, char** argv) {
	char optString[] = {"d:N:R:r:a:b:"};
	char* data = "poisson.dat";
	FILE *LOG; 

	srand(time(NULL));

	int opt; 
	opt = getopt( argc, argv, optString );
	while( opt != -1 ) {
        switch( opt ) {
		case 'd':
			data = optarg; 
			break;
		case 'R':  
			sscanf(optarg, "%lf", &R[1]); 
           	break;
		case 'N': 
			sscanf(optarg, "%d", &N); 
           	break;
			
		case 'r': 
			sscanf(optarg, "%lf", &R[0]); 
			sscanf(argv[optind], "%lf", &R[1]);
            optind++;
			break;
			
		case 'a':
			sscanf(optarg, "%lf", &A[0]); 
			sscanf(argv[optind], "%lf", &A[1]);
            optind++;
            break;

		case 'b':
			sscanf(optarg, "%lf", &B[0]); 
			sscanf(argv[optind], "%lf", &B[1]);
            optind++;
        	break;
		
		default:
			printf("What?\n");
			break;
        }
        opt = getopt( argc, argv, optString );
	}
    A[0] *= d2r;  A[1] *= d2r; 
    B[0] *= d2r;  B[1] *= d2r; 
    get_poisson_from_square( data );
    /*
	LOG = fopen("distrub.log", "a");
	fprintf(LOG, "******** Poisson *******\n" );
	fprintf(LOG, "File name = %s\n", data);
	fprintf(LOG, "N : %d\n", N);
	if( (Ra[0]==Ra[1]) & (Dec[0] == Dec[1]) ) {
		fprintf(LOG, "R : %lf\n", R)  ;
		fprintf(LOG, "Whole sphere\n");
		get_poisson_sphere(data);	
	} else {
		fprintf(LOG, "R   : %lf %lf\n", r[0],   r[1])  ;
		fprintf(LOG, "RA  : %lf %lf\n", Ra[0]*r2d,  Ra[1]*r2d) ;
		fprintf(LOG, "DEC : %lf %lf\n", Dec[0]*r2d, Dec[1]*r2d);	
		get_poisson_sector(data);
	}
	fclose(LOG);
    */
	return 0;
}

