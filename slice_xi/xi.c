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

char *input_suff, *out_suff, *name, *output; 
int N, NS, NP; // how many points
double q = 1, qran = 10;
double **p, **P;  // Array of points. p[i][N] and P[i][N]
double **s, **S;  // Array of points. slice part 
double **POISS;  // Array of points. RANDOM
int    *slice_index; // connection of indexes slice -- main set.
double R0 = 50 , dr = 10, r[2]; 

int grid_n = 1000;
double *grid, grid_h1, grid_R[2] = {0.001, M_PI}; // grid itself, h+1 and limits. 

typedef unsigned long long int myint; 
myint  **DD, **DR, **RR; 
myint  max = 1;

// double A[2] = {0, 90}, B[2] = {0, 90}, R[2] = {1, 300};
double A[2] = {0, 0}, B[2] = {0, 0}, R[2] = {0, 300};
double r2d = 180/M_PI, d2r = M_PI / 180;
double n0, solid_angle, density;

char *limits_file = "limits.txt";
double slice_limits[100][2];
int slice_limits_N;
//double step[];; 

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
    if( ( A[0] == A[1] ) & (B[0] == B[1]) ) return 0; 
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

void
count(  int same, 
        double **S1, int N1, 
        double **S2, int N2, 
        myint ** pairs ) {
	
    int i, j, k;
    double rd;
    
    int start; 
    
	for(i=0; i<N1; i++) {
        printf("\r%d%%", (int) 100 * i / N1 );
		for(j=(i+1)*same; j < N2; j++) {
            rd = acos(
                    sin( S1[i][1] ) * sin( S2[j][1] ) + 
                    cos( S1[i][1] ) * cos( S2[j][1] ) * cos( S1[i][0] - S2[j][0] ) ); 
            k = index_grid( rd ); 
            if( same == 1 ) 
                pairs[i][k] += 2; 
            else 
                pairs[i][k] ++; 

		}
	}
    printf("\n");
}

/*
void
write_gamma( FILE * foutput ) {
    int i,j,k; 
    double tmp; 
    double * mean  = (double *) calloc( grid_n, sizeof(double) );
    myint  * count = (myint *) calloc( grid_n, sizeof(myint) );
    double    * nc    = (double *) calloc( grid_n, sizeof(double) );
    double    * zer0  = (double *) calloc( grid_n, sizeof(double) );
    double    * zer1  = (double *) calloc( grid_n, sizeof(double) );
    double    * zer2  = (double *) calloc( grid_n, sizeof(double) );
    double    * delta  = (double *) calloc( grid_n, sizeof(double) );

    double closest = 0;
    for( i=0; i<NS; i++ ) {
        closest += neighbor[i]; 
    }
    closest /= NS; 
    

    for( k=0; k<grid_n; k++ ) {
        nc[k] = 0; count[k] = 0;
        for(i=0; i<NS; i++) {
            // if( s_R_max(i) < 0.3 ) continue; 
            if( s_R_max_ind(i)-1 >= k) {
                count[k] += m[i][k]; 	
                if( m[i][k] < 0 ) printf("M < 0. OPPA \n");
                if( count[k] < 0 ) printf(" count too small HAHA \n");
                nc[k]++; 
                if( m[i][k] == 0 ) zer0[k]++;
                if( m[i][k] <= 1 ) zer1[k]++;
                if( m[i][k] <= 2 ) zer2[k]++;
            } 	
        }
        mean[k] = (nc[k] != 0) ? (double) count[k] / nc[k]  :  -1;
        if( mean[k] < 0 & grid[k] < 0.38 ) printf("k = %d, r = %lf mean = %lf\n", k, grid[k], mean[k] ); 
        if( mean[k] < 0 ) continue; //printf("k = %d, r = %lf mean = %lf\n", k, grid[k], mean[k] );
        for( i=0; i<NS; i++ ) {
            tmp = 0;
            if( s_R_max_ind(i)-1 >= k ) {
                tmp += pow( m[i][k] - mean[k] , 2);
            }
            delta[k] = sqrt( (double) 1. / (nc[k]-1) * tmp );
        }
        nc[k] = (double) nc[k] / NS * 100;
        zer0[k] = (double) zer0[k] / NS * 100; 
        zer1[k] = (double) zer1[k] / NS * 100; 
        zer2[k] = (double) zer2[k] / NS * 100; 
    }

	for(k=0; k<grid_n; k++) {
		if( zer0[k] >= 100 ) continue; 
        if( nc[k] <= 0 ) continue; 
        if( mean[k] == 0 ) continue; 
		fprintf(foutput, "%le ", grid[k]  ); 
		fprintf(foutput, "%le ", mean[k]  / ( 2*M_PI * ( 1 - cos(grid[k]) ) ) / density );
		fprintf(foutput, "%le ", delta[k] / ( 2*M_PI * ( 1 - cos(grid[k]) ) ) );
		fprintf(foutput, "%le ", grid[k] / closest  ); 
		fprintf(foutput, "%le ", nc   [k]  ); 
		fprintf(foutput, "%le ", zer0 [k]  ); 
		fprintf(foutput, "\n");
	}


    free(mean );
    free(count);
    free(nc   );
    free(zer0 );
    free(zer1 );
    free(zer2 );
}

void
write_foo( double scale ) {
    int i,j,k;
    int scale_ind; 
    scale_ind = index_grid( scale ); 

    char foo_name[80];
    sprintf( foo_name, "foo_%s_%03.f-%03.f_%1.3f.dat", out_suff, r[0], r[1], scale); 
    FILE *foo_file = fopen( foo_name, "w"); 

    for( i=0; i<NS; i++ ) {
        if( m[i][scale_ind]==0) continue;
        fprintf(foo_file, "%f %f %ld\n", S[i][0], S[i][1], (long int)  m[i][ scale_ind ]);    
    }
    fclose( foo_file );
}
*/

void
write_xi( FILE * foutput ) {
	int k, i; 
    myint Dd = 0, Rr = 0, Dr = 0;
	for(k=0; k<grid_n; k++) {
        Dd = 0; Rr = 0; Dr = 0;
            for( i = 0; i < NS; i++ ) { 
                Dd += DD[i][k]; 
                Dr += DR[i][k];
            }
            for( i = 0; i < NP; i++ ) {
                Rr += RR[i][k];
            }
        /// :if(Dd == 0 || Dr == 0 || Rr == 0) continue; 
		fprintf(foutput, " %le", grid[k] ); 
		fprintf(foutput, " %le", (double) Dd / Rr * ((double) NP / NS)  );  
		fprintf(foutput, " %le", (double) ( NP / NS * ( (double) NP / NS * (double) Dd / Rr - 2. * (double) Dr / Rr) + 1 ) + 1 );
        fprintf(foutput, " %le", (double) ( NP-1 ) / NS / Rr * ( (double) NP / ( NS-1 ) * Dd - (double) 2 * Dr ) + 1 + 1);
		fprintf(foutput, " %le", (double) Dd / NS);  
		fprintf(foutput, " %le", (double) Dr / NS);  
		fprintf(foutput, " %le", (double) Rr / NP);  
		fprintf(foutput, "\n");
    }
    
}

void
write_results() {
	FILE *foutput; 
	int k, i; 

    char output[80]; 
	sprintf(output, "%s_xi_%03.f-%03.f.dat", out_suff, r[0], r[1]);
    printf("write to %s\n", output); 

	foutput = fopen(output, "w");
    fprintf(foutput, "# %s %5f %5f %5f %5f %5f %5f\n", name, A[0], A[1], B[0], B[1], R[0], R[1]);
    fprintf(foutput, "# %10s %10s %10s %10s %10s\n", "r1", "r2", "N", "dens", "neig");
    fprintf(foutput, "# %10f %10f %10d %10f\n", r[0], r[1], NS, density);
    fprintf(foutput, "#\n");
    //foutput = fopen("gamma_2d.dat", "w");

    write_xi ( foutput ); 

	fclose(foutput);
    
}


void
make_grid() {
    grid[0] == grid_R[0]; 
	if(grid[0] == 0) grid[0] = 1e-10;
	grid_h1 = pow(grid_R[1]/grid[0], (double) 1./((double) grid_n - 1));
	int k;
    //printf("%lf %lf\n", grid[0], grid_h1);
	for( k=1; k<grid_n; k++)	{
		grid[k] = grid[k-1] * grid_h1;
	}	
    //printf("pewpew\n");
}

int
index_grid(double r) {
	int i;
	i = (int) (log(r / grid[0] ) / log(grid_h1)); 
	if(i<=0) {
        i = 0;
    }
	return i;
}

void
init_begining() {
	A[0]=A[0]*M_PI/180;  
	A[1]=A[1]*M_PI/180;  
	B[0]=B[0]*M_PI/180;  
	B[1]=B[1]*M_PI/180;  
	int i, k;
	// points. N of them. 
	
	p = (double **) calloc(N, sizeof(double *));
	for(i=0; i<N; i++)  {
		p[i] = (double *) calloc(3, sizeof(double));
	}
	P = (double **) calloc(N, sizeof(double *));
	for(i=0; i<N; i++)  {
		P[i] = (double *) calloc(3, sizeof(double));
	}
}

void
read_data() {
	FILE 	*file; 
	char	input_name[80];
	double 	x, y, z, a, b, r;
	int i, oor = 0, hack; // out of range;	
	int NMAX = N;
	double alim[] = {360*d2r, 0}, blim[] = {90*d2r, -90*d2r}, rlim[] = {1e10, 0};
	
	sprintf(input_name, "%s.dat", input_suff);
	file = fopen(input_name, "r");
	if( file == NULL ) {
		printf("Can't find file %s.\n", input_name);
		exit(10);
	}
	double foo;
	for(i=0; i<N; i++) {
		foo=((double) rand() / RAND_MAX);
        if( SELECT_COORDINATES == 0 ) {
            //printf("!!! read cartesian coordinates !!!\n");
            hack = fscanf(file, "%lf %lf %lf", &(p[i][0]), &(p[i][1]), &(p[i][2]));
            spherical_coordinats(1, &(p[i][0]), &(p[i][1]), &(p[i][2]), &(P[i][0]), &(P[i][1]), &(P[i][2]));
        } else { 
            //printf("!!! read spherical coordinates !!!\n");
            hack = fscanf(file, "%lf %lf %lf", &(P[i][0]), &(P[i][1]), &(P[i][2]));
            P[i][0] *= d2r; P[i][1] *= d2r; 
            cartesian_coordinats(1, &(p[i][0]), &(p[i][1]), &(p[i][2]), &(P[i][0]), &(P[i][1]), &(P[i][2]));
        } 

		if( P[i][0] < alim[0] ) alim[0] = P[i][0]; 
		if( P[i][0] > alim[1] ) alim[1] = P[i][0]; 
		if( P[i][1] < blim[0] ) blim[0] = P[i][1]; 
		if( P[i][1] > blim[1] ) blim[1] = P[i][1]; 
		if( P[i][2] < rlim[0] ) rlim[0] = P[i][2]; 
		if( P[i][2] > rlim[1] ) rlim[1] = P[i][2]; 
		if( (inside_sp(P[i][0], P[i][1], P[i][2]) != 0) || ( q < foo))  {
			oor++;
			i--; N--;
		} 	//else
			//printf("%lf\n", P[i][1]);
	}
	fclose(file);

    // Parameters of selected objects
    if( (A[1] == A[0]) & (B[0] == B[1]) ) {
        solid_angle = 4 * M_PI;
    } else {
	    solid_angle = (A[1] - A[0]) * (sin(B[1]) - sin(B[0]));
    }
	// n0 = (double) 3. * N / (solid_angle * (pow(R[1],3) - pow(R[0],3)) );
    n0 = ( double ) N / solid_angle;
	printf("%d selected.\n", N);
	printf("A : %+5le %+5le -> %+5le %+5le \n", alim[0]*r2d, alim[1]*r2d, A[0]*r2d, A[1]*r2d );
	printf("B : %+5le %+5le -> %+5le %+5le \n", blim[0]*r2d, blim[1]*r2d, B[0]*r2d, B[1]*r2d );
	printf("R : %+5le %+5le -> %+5le %+5le \n", rlim[0]    , rlim[1]    , R[0]    , R[1]     );
	printf("Solid angle : %le\nDensety : %le\n", solid_angle, n0);
	
	sprintf(input_name, "%s_c2.dat", name);
	file = fopen(input_name, "w"); 
	for(i=0; i<N; i++) {
		//fprintf(file, "%lf %lf %lf %lf %lf %lf\n", p[i][0],  p[i][1], p[i][2], P[i][0],  P[i][1], P[i][2]);
		fprintf(file, "%lf %lf %lf\n", p[i][0],  p[i][1], p[i][2]);
	}
	fclose(file);
}

void
init_second() {
	int i,k;
    grid = (double *) calloc(grid_n, sizeof(double));

	DD = (unsigned long long int **) calloc(N, sizeof(unsigned long long int *));
	for(i=0; i<N; i++) {
		DD[i] = (unsigned long long int *) calloc(grid_n, sizeof(unsigned long long int));
	}

	DR = (unsigned long long int **) calloc(N, sizeof(unsigned long long int *));
	for(i=0; i<N; i++) {
		DR[i] = (unsigned long long int *) calloc(grid_n, sizeof(unsigned long long int));
	}

	RR = (unsigned long long int **) calloc(N * qran, sizeof(unsigned long long int *));
	for(i=0; i<N * qran; i++) {
		RR[i] = (unsigned long long int *) calloc(grid_n, sizeof(unsigned long long int));
	}
	
	POISS = (double **) calloc(N*qran, sizeof(double *));
	for(i=0; i<N*qran; i++)  {
		POISS[i] = (double *) calloc(3, sizeof(double));
	}
	
	s = (double **) calloc(N, sizeof(double *));
	for(i=0; i<N; i++)  {
		s[i] = (double *) calloc(3, sizeof(double));
	}

	S = (double **) calloc(N, sizeof(double *));
	for(i=0; i<N; i++)  {
		S[i] = (double *) calloc(3, sizeof(double));
	}
    slice_index = (int *) calloc( N, sizeof( int ) );
}

void
select_slice( ) {
    // We need to count 2d slice. 
    // 1st, we need to select propper points. 
    int i, j=0; 
    for( i=0; i<N; i++ ) {
        if( (P[i][2] > r[0]) & (P[i][2] < r[1]) ) {
            s[j] = p[i];
            S[j] = P[i]; 
            j++;
        }
    }
    j--;
    density = (double) j / (solid_angle);
    NS = j;
    printf("slice %e %e -- %d\n", r[0], r[1], NS);
    /*
    char slice_file[80];
    sprintf(slice_file, "slice_%3e_%3e.dat", r[0], r[1]);
	FILE* file = fopen(slice_file, "w"); 
	for(i=0; i<NS; i++) {
		fprintf(file, "%lf %lf %lf %lf %lf %lf %lf\n", s[i][0],  s[i][1], s[i][2], S[i][0],  S[i][1], S[i][2], s_R_max(i));
	}
	fclose(file);
    */
}

void
whipe() {
    int i, j, k;
    for( i=0; i<N; i++ ) {
        for( k=0; k<grid_n; k++ ) {
            DD[i][k] = 0;
        }
    }
    for( i=0; i<N; i++ ) {
        for( k=0; k<grid_n; k++ ) {
            DR[i][k] = 0;
        }
    }
    for( i=0; i<NP; i++ ) {
        for( k=0; k<grid_n; k++ ) {
            RR[i][k] = 0;
        }
    }
}

void
read_limits() {
    limits_file;
    FILE * F_LIMITS;
    char line[80];
    int readed, i = 0; 

    F_LIMITS = fopen( limits_file, "r" ) ; 
    if( F_LIMITS == NULL ) {
        printf(" can't open file %s \n", limits_file);
        exit(20); 
    }
    
    i=-1;
    while(fgets(line, 80, F_LIMITS) != NULL) {
        i++; 
        if( i >= 100 ) { 
            printf("slice_limits is too small.\n");
            exit(42);
        }
        readed = sscanf (line, "%lf %lf", &(slice_limits[i][0]), &(slice_limits[i][1]));
        if( readed != 2 ) {
            printf("I need lines of 2 float in limits_file, is it that hard?\n");
            exit(42);
        }
    }
    fclose(F_LIMITS);
    slice_limits_N = i; 
}

void
generate_poisson() { 
    NP = NS * qran; 
	int i; 
    double sinb[2];

    sinb[0] = sin( B[0] ); 
    sinb[1] = sin( B[1] ); 
    
    printf("A: %3.1f %3.1f\n", A[0] * 180 / M_PI, A[1] * 180 / M_PI);
    printf("B: %3.1f %3.1f\n", B[0] * 180 / M_PI, B[1] * 180 / M_PI);

    for( i=0 ; i < NP ; i ++ ) { 
        POISS[i][0] = ((double) rand() / RAND_MAX) * ( A[1] - A[0] ) + A[0] ;
        POISS[i][1] = asin ( ((double) rand() / RAND_MAX) * ( sinb[1] - sinb[0] ) + sinb[0] );
    }
    /*
    FILE  * file = fopen("pos.dat", "w") ;
    for( i=0 ; i < NP ; i ++ ) { 
        fprintf( file, "%lf %lf 1\n", POISS[i][0], POISS[i][1] ) ;
    }
    fclose( file );
    */
}

void
gogogo() {
    double old_N; 

	init_begining();
	read_data();
    if( N == 0 ) exit(10);
	init_second();
    make_grid();
    read_limits(); 
    
    //double step[] = {1,2,3,10} ;
    int i;
    for( i=0; i <= slice_limits_N; i++ ) {
        old_N = N;
        r[0] = MIN( slice_limits[i][0], slice_limits[i][1] ); 
        r[1] = MAX( slice_limits[i][0], slice_limits[i][1] ); 
        // r[0] = slice_limits[i][0];
        // r[1] = slice_limits[i][1]; 
        printf("%e %e\n", slice_limits[i][0], slice_limits[i][1]);
        
        
        dr = r[1] - r[0];
        R0 = r[0] + dr/2;
        select_slice(); 
        generate_poisson(); 
        count( 1, S, NS, S, NS, DD );
        count( 0, S, NS, POISS, NP, DR );
        count( 1, POISS, NP, POISS, NP, RR );
        write_results();
        whipe();
        N = old_N;
    }
}

int main( int argc, char *argv[] ) {
	char optString[] = {"d:N:q:a:b:r:n:R:co:p:"};
	int opt, got_out = 0; 
	opt = getopt( argc, argv, optString );
	while( opt != -1 ) {
		switch( opt ) {
            // READ FILE
			case 'd':
				name = optarg; 
				input_suff = optarg; 
				break;
			case 'N':
				sscanf(optarg, "%d", &N); 
				break;
            // Read just some portion of a file
			case 'q': 
				sscanf(optarg, "%lf", &q); 
				break;
			case 'p': 
				sscanf(optarg, "%lf", &qran); 
				break;
            // 2d AREA 
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
			case 'l':
				limits_file = optarg; 
				break;
			case 'R':
				sscanf(optarg, "%lf", &(grid_R[0])); 
				sscanf(argv[optind], "%lf", &(grid_R[1]));
				optind++;
				break ;
            case 'n':
                sscanf(optarg, "%d", &grid_n); 
                break;
            case 'c': 
                SELECT_COORDINATES = 1; // spherical
                printf("I'll work with spherical coordinates.\n");
                break;
			case 'o':
				got_out = 1;
				output = optarg; 
                out_suff = optarg;
                break;
			default:
				printf("What? %s\n", optarg);
				break;
    	}
	    opt = getopt( argc, argv, optString );
	}
	if( got_out != 1 ) {
	    out_suff = input_suff;
    }
	gogogo();
	return(0);	
}

    //int i, j, scale;
    //double rd; 
    //FILE * file; 
    //char filename[80] ;
    //scale = 10;
    //for( i=0; i<10; i++ ) {
    //    sprintf(filename, "circle_%d_%d.dat", i, scale);
    //    file = fopen( filename, "w" ); 
    //    for( j=0; j<N; j++ ) {
    //        if( i == j ) continue; 
	//		rd = R0 * acos(
    //                sin( S[i][1] ) * sin( S[j][1] ) + 
    //                cos( S[i][1] ) * cos( S[j][1] ) * cos( S[i][0] - S[j][0] ) ); 
    //        if( rd < scale ) {
    //            fprintf( file, "%le %le %le %le %le %le\n", s[j][0], s[j][1], s[j][2], S[j][0], S[j][1], S[j][2]); 
    //        }
    //    }
    //    fclose( file );
    //}
