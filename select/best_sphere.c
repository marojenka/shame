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

int SELECT_COORDINATES = 0; // CARTESIAN

char *input_suff, *name, *output; 
int N, n; // how many points
double **p;    // Array of points. p[i][N]
double R = 1; // radii of desired sphere 
double A = 2; // size of used box
int INDEX = 0; // INDEX of new center;

double q = 1;

void
find_INDEX() {
	int i, j, k;
	double rd;
    double x,   y,  z;
    double dx, dy, dz;
    int number_of_points, max_number_of_points = 0; 

	for(i=0; i<n; i++) {
        printf("\r%d%%", (int) 100 * i / n );
        x = p[i][0]; 
        y = p[i][1]; 
        z = p[i][2]; 
        if( ( fabs(x) > A - R ) | 
            ( fabs(y) > A - R ) | 
            ( fabs(z) > A - R )   ) {
            continue; 
        }
        number_of_points = 0;
		for(j=0; j<n; j++) {
            if( i == j ) continue;

            dx = fabs(x - p[j][0]);
            dy = fabs(y - p[j][1]);
            dz = fabs(z - p[j][2]);
            if( ( dx > R ) | 
                ( dy > R ) | 
                ( dz > R )   ) {
                // if points are out of box with side 2*R
                // then distance between them clearly is bigger than R.
                continue;
            } else {
                rd = sqrt( dx*dx + dy*dy + dz*dz ); 
                if( rd < R ) {
                    number_of_points++;
                }
            }
		}
        if( number_of_points > max_number_of_points ) {
            max_number_of_points = number_of_points; 
            INDEX = i; 
        }
	}
	printf("\n");
}

void
read_data() {
	FILE 	*file; 
	char	input_name[80];
	double 	x, y, z, r;
	int hack; // out of range;	
	int NMAX = N;
	
	sprintf(input_name, "%s.dat", input_suff);
	file = fopen(input_name, "r");
	if( file == NULL ) {
		printf("Can't read file %s.\n", input_name);
		exit(10);
	}

	double foo;
    int i = -1, k;
	for(k=0; k<N; k++) {
		foo=((double) rand() / RAND_MAX);
        hack = fscanf(file, "%lf %lf %lf", &x, &y, &z);

		if( q > foo )  {
            i++; 
            p[i][0] = x; 
            p[i][1] = y; 
            p[i][2] = z;
		} 	
    }
	fclose(file);
    if( i == -1 ) { 
        printf("empty\n");
        exit(30);
    }
    n = i+1;
    printf("N centers = %d\n", n);
}

void
init_begining() {
	int i, k;
	p = (double **) calloc(N, sizeof(double *));
	for(i=0; i<N; i++)  {
		p[i] = (double *) calloc(3, sizeof(double));
	}
}

void
write_sphere() {
	int i, j, k, hack; 
    double  X,  Y,  Z;
    double  x,  y,  z;
    double dx, dy, dz;
    double rd; 
	FILE    *foutput; 
	FILE 	*file; 
	char	input_name[80];


    printf("write to %s\n", output);
	foutput = fopen(output, "w");
	sprintf(input_name, "%s.dat", input_suff);
	file = fopen(input_name, "r");

    i = INDEX;
    X = p[i][0]; 
    Y = p[i][1];
    Z = p[i][2]; 
    for( j=0; j<N; j++ ) {
        hack = fscanf(file, "%lf %lf %lf", &x, &y, &z);
        dx = fabs(x - X);
        dy = fabs(y - Y);
        dz = fabs(z - Z);
        if( ( dx > R ) | 
            ( dy > R ) | 
            ( dz > R )   ) {
            continue;
        } else {
            rd = sqrt( dx*dx + dy*dy + dz*dz ); 
            if( rd < R ) {
                fprintf(foutput, "%le %le %le \n", x - X, y - Y, z - Z);
            }
        }
    }
    fclose(foutput);
    fclose(file);
}

void
_gogogo() {
    int INDEX; 
	init_begining();
	read_data();
    find_INDEX();
    if( INDEX == -1 ) {
        printf("Something wrong\n");
        exit(30);
    }
    write_sphere(); 
}

int main( int argc, char *argv[] ) {
	char optString[] = {"d:N:q:R:o:A:"};
	int opt, got_out = 0; 
	opt = getopt( argc, argv, optString );
	while( opt != -1 ) {
		switch( opt ) {
			case 'd':
				name = optarg; 
                input_suff = optarg; 
				break;
			case 'N':
				sscanf(optarg, "%d", &N); 
				break;
			case 'q': 
				sscanf(optarg, "%lf", &q); 
				break;
			case 'R': 
				sscanf(optarg, "%lf", &R); 
				break;
			case 'A': 
				sscanf(optarg, "%lf", &A); 
				break;
			case 'o':
				got_out = 1;
				output = optarg; 
				break;

			default:
				printf("What? %d\n", opt);
                exit(42);
				break;
    	}
	    opt = getopt( argc, argv, optString );
	}
	if( got_out == 1 ) {
        name = output; 
    }
	char tmp[80]; 
	sprintf(tmp, "%s_sphere.dat", name);
	output = tmp;
	_gogogo();
	return(0);	
}

