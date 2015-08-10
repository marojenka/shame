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
int N; // how many points
double **p;  // Array of points. p[i][N]
double **P;  // Array of points. p[i][N]
double *R_max; // Rmax ofc
int    *R_max_ind; // Rmax compared to grid. grid[R_max_ind] < R_max[i] < grid[R_max_ind + 1] ???

int grid_n = 800;
double *grid, grid_h1, grid_R[2] = {1e-3, 1e3}; // grid itseld, h+1 and limits. 
double *max_dist; 

typedef unsigned long long int myint; 
myint  **m; 
double *neighbor; 

double A[2] = {0,0}, B[2] = {0, 0}, R[2] = {0, 300};
double q =1;
double r2d = 180/M_PI, d2r = M_PI / 180;
double n0, solid_angle;

int index_grid(double r);
typedef struct statK {
    int Nc; 
    myint c_sh, c_sp; 
} statK;

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

int
inside_ca(double x,double y, double z) {
	double a,b,r; 
	spherical_coordinats(1,&x,&y,&z,&a,&b,&r);
	return inside_sp(a,b,r);
}

double
update_max_dist(int i, int j, int K, double rd) {
	int k;
    for(k=0; k<=MIN(R_max_ind[i]+1, R_max_ind[j]+1); k++) {
		if( max_dist[k] < rd ) max_dist[k] = rd; 
    }
}

// void
// add(int i, int j, int K) {
// 	int k;
//     //for(k=K+1; k<=MIN(R_max_ind[i]+1, grid_n); k++) {
// 	//for(k=K+1; k<=(grid_n); k++) {
// 		m[i][k] ++;
// 	}
// }

void
count() {
	int i, j, k;
	double rd;
	for(i=0; i<N; i++) {
	    //printf("\ri=%d\t%d%%", i, (int) 100 * i / N );
        printf("\r%d%%", (int) 100 * i / N );
		for(j=i+1; j<N; j++) {
            rd = sqrt(pow(p[i][0]-p[j][0],2)+pow(p[i][1]-p[j][1],2)+pow(p[i][2]-p[j][2],2)); 
			if( rd < neighbor[i] ) neighbor[i] = rd;
			if( rd < neighbor[j] ) neighbor[j] = rd;
            k = index_grid(rd); // grid[k] < rd < grid[k+1]

            if( k < R_max_ind[i]) {
                m[i][k] ++ ;
				//add(i, j, k ); 
            }
			if( k < R_max_ind[j]) {
                m[j][k] ++ ;
				//add(j, i, k ); 
            }
            // if( (rd < R_max[i]) & (rd < R_max[j]) ) {
            // update_max_dist(i, j, k, rd);
            // }
		}
	}
	printf("\n");
}

// double
// ave_N(int K) {
// 	int i, j, k;
// 	int nc; 
// 	myint count;
// 	count = 0; nc = 0;
// 	for(i=0; i<N; i++) {
//         if( R_max_ind[i] <= K ) 
//             continue; 
// 
//         nc++; 
// 
//         for( k=0; k<=K; k++ ) {
// 			count += m[i][k]; 	
//         }
// 		//if(R_max_ind[i]-1 >= k) {
// 		//	count += m[i][k]; 	
// 		//	nc++; 
// 		//} 	
// 	}
// 	//if( nc == 0 | count == 0 ) return -1; 
// 	if( nc == 0 ) return -1; 
// 	return (double) count/nc;
// }
// 
// double
// ave_N_shell(int K) {
// 	int i, j;
// 	int nc; 
// 	myint count;
// 	count = 0; nc = 0;
// 	for(i=0; i<N; i++) {
//         if( R_max_ind[i] <= K ) 
//             continue; 
//         nc++; 
// 		count += m[i][k]; 	
//         
// 		//if(R_max_ind[i]-1 >= k) {
// 		//	count += m[i][k]; 	
// 		//	nc++; 
// 		//} 	
// 	}
// 	//for(i=0; i<N; i++) {
// 	//	if(R_max_ind[i]-1 >= k) {
// 	//		count += m[i][k] - m[i][k-1]; 	
//     //        if( m[i][k] - m[i][k-1] < 0 ) 
//     //            printf("PEW\n");
// 	//		nc++; 
// 	//	} 	
// 	//}
// 	//if( nc == 0 | count == 0 ) return -1; 
// 	if( nc == 0 ) return -1; 
// 	return (double) count/nc;
// }
// 

// double
// calc_nc(int K) {
// 	int i;
// 	double nc; 
// 	nc = 0.0;
// 	for(i=0; i<N; i++) {
// 		if(R_max_ind[i] <= K)
//         // if(m[i][k] > 0) 
// 			nc++; 
// 	}
// 	return (double) nc;
// }
// 
// double
// sd(int k) {
// 	int i;
// 	int nc; 
// 	int count;
// 	count = 0; nc = 0;
// 	for(i=0; i<N; i++) {
// 		if(R_max_ind[i] >= k) {
// 			count += m[i][k]; 	
// 			nc++; 
// 		} 	
// 	}
// 
// 	double mean = (double) count / nc; 
// 	double delta = 0;
// 
// 	for(i=0; i<N; i++) {
// 		if(R_max_ind[i] >= k) {
// 			delta += pow( (double) m[i][k] - mean , 2);
// 		} 	
// 	}
// 	return (double) sqrt( (double) 1. / (nc-1) * sqrt(delta) );
// }
// 
// double
// sigma(int k) {
// 	int i;
// 	int nc; 
// 	int count;
// 	count = 0; nc = 0;
// 	for(i=0; i<N; i++) {
// 		if(R_max_ind[i] >= k) {
// 			count += m[i][k]; 	
// 			nc++; 
// 		} 	
// 	}
// 
// 	double mean = (double) count / nc; 
// 	double delta = 0;
// 
// 	for(i=0; i<N; i++) {
// 		if(R_max_ind[i] >= k) {
// 			delta += pow( (double) m[i][k] - mean , 2);
// 		} 	
// 	}
// 	return (double) ( (double) 1. / (nc-1) * delta / pow(mean,2) );
// }

statK
foo(int K) {
	int i, j, k;
    statK F; 
	F.c_sh = 0;
    F.c_sp = 0;
    F.Nc = 0; 
	for(i=0; i<N; i++) {
        if( R_max_ind[i] <= K ) 
            continue; 
        F.Nc++; 
        for( k=0; k<=K; k++ ) {
			F.c_sp += m[i][k]; 	
        }
        F.c_sh = m[i][K];
	}
	//if( nc == 0 | count == 0 ) return -1; 
	return F;
}

double volume_sp(double r) { return( 4./3. * M_PI * pow(r,3)); }
double volume_sh(double r1, double r2) { return( 4./3. * M_PI * pow(MAX(r1,r2),3) - pow(MIN(r1,r2),3)); }

void
write_result() {
	FILE *foutput; 
	int k, i; 
    
    double closest = 0;
    for( i=0; i<N; i++ ) {
        closest += neighbor[i]; 
        //if( neighbor[i] == 0 ) printf(" OPPA \n");
    }
    closest /= N;

    printf("write to %s\n", output); 
	foutput = fopen(output, "w");
    fprintf(foutput, "# %s %5f %5f %5f %5f %5f %5f\n", name, A[0], A[1], B[0], B[1], R[0], R[1]);
    fprintf(foutput, "# %10s %10s %10s\n","N",  "dens", "neig");
    fprintf(foutput, "# %10d %10f %2.2f\n", N, n0,  closest);
    fprintf(foutput, "#\n");
    //foutput = fopen("gamma_2d.dat", "w");

    int HAI = index_grid(10.);
    //double HAX = ave_N(HAI) / (4./3.*M_PI*pow(grid[HAI],3));
    statK F;
	for(k=1; k<grid_n-1; k++) {
        F = foo(k); 
        if( F.Nc == 0 ) continue; 
		fprintf(foutput, "%le ", grid[k] ); 
		fprintf(foutput, "%le ", (double) F.c_sp / F.Nc / volume_sp(grid[k]) ); 
		fprintf(foutput, "%le ", (double) F.c_sh / F.Nc / volume_sh(grid[k], grid[k+1]) ); 
		fprintf(foutput, "%d  ", F.c_sp ); 
		fprintf(foutput, "%d  ", F.c_sh ); 
		fprintf(foutput, "%d  ", F.Nc ); 
		//fprintf(foutput, "%le ", sigma(k) ); 
		//fprintf(foutput, "%le ", calc_nc(k) / N * 100 ); 
		//fprintf(foutput, "%le ", ave_N(k) ); 
		//fprintf(foutput, "%le ", (double) sqrt( 1. / ave_N(k) ) ); 
		//fprintf(foutput, "%le ", grid[k] / closest ); 
        //fprintf(foutput, "%le ", max_dist[k]);
		fprintf(foutput, "\n");
	}
	fclose(foutput);
    /*
    FILE * fpdf;
    double scales [] = {1, 5, 10, 20, 30, 50, 80, 100}; 
    char pdfname[80];
    int j;
    for( j=0; j<8; j++ ) {
        sprintf( pdfname, "pdf_%s_%03.0f.dat", name, scales[j] );
        fpdf = fopen(pdfname, "w");
        k  = index_grid( scales[j] ) ;
        for( i=0; i<N; i++ ) {
            if(R_max_ind[i]-1 >= k) {
                fprintf( fpdf, "%ld \n", (long int) m[i][k] );
            } 	
        }
        fclose(fpdf); 
    }
    */
}

void
write_mad() {
	FILE *foutput; 
    char output[80];
    double max, min, h; 
	int i, j, k; 
    const int  length = 100; 
    int hist[length]; 

	sprintf(output, "%s_mad.dat", name);
	foutput = fopen(output, "w");

    /*
    max = m[0][k]; min = m[0][k];
    for( k = 0; k < grid_n; k++ ) {
        for( i = 0; i < N; i++ ) {
            if( m[i][k] > max ) max = m[i][k]; 
            if( m[i][k] < min ) min = m[i][k];
        }
    }
    //if( max == min ) return 0;
    h = ( max - min ) / length; 
    for( k = 0; k < grid_n; k++ ) {
        for( j = 0; j < length; j++ ) 
            hist[ j ] = 0;

        for( i = 0; i < N; i++ ) {
            if( R_max_ind[i] >= k )  
                hist[ (int) floor((m[i][k] - min ) /  h) ] ++; 
        }

        for( j = 0; j < length; j++ )
            fprintf(foutput, "%d ", hist[j]);
        fprintf(foutput, "\n");
    }
    */
    /*
    for( k = 0; k < grid_n; k++ ) {
        max = m[0][k]; min = m[0][k];
        for( i = 0; i < N; i++ ) {
            if( m[i][k] > max ) max = m[i][k]; 
            if( m[i][k] < min ) min = m[i][k];
        }
        if( max == min ) continue;
        h = ( max - min ) / length; 
        for( j = 0; j < length; j++ ) 
            hist[ j ] = 0;

        for( i = 0; i < N; i++ ) {
            hist[ (int) floor((m[i][k] - min ) /  h) ] ++; 
        }
        printf("\n");

        for( j = 0; j < length; j++ )
            fprintf(foutput, "%d ", hist[j]);
        fprintf(foutput, "\n");
    }
    */
    fprintf(foutput, "%d %d %d %d %d %d ", 0, 0, 0, 0, 0, 0);
    for(k=0; k<grid_n; k++) { fprintf(foutput, "%le ", grid[k]); } fprintf(foutput, "\n"); 
    fprintf(foutput, "%d %d %d %d %d %d ", 0, 0, 0, 0, 0, 0);
    for(k=0; k<grid_n; k++) { fprintf(foutput, "%le ", (4./3.*M_PI*pow(grid[k],3))); } fprintf(foutput, "\n"); 
    fprintf(foutput, "%d %d %d %d %d %d ", 0, 0, 0, 0, 0, 0);
    for(k=0; k<grid_n; k++) { fprintf(foutput, "%le ", (4./3.*M_PI*pow(grid[k],3) - 4./3.*M_PI*pow(grid[k-1],3))); } fprintf(foutput, "\n"); 

    for( i=0; i<N; i++ ) {
        fprintf(foutput, "%d %le %le ", R_max_ind[i], R_max[i], P[i][2]);
        fprintf(foutput, "%le %le %le ", p[i][0], p[i][1], p[i][2]);
        for( k=0; k<grid_n; k++ ) {
            fprintf(foutput, "%llu ", m[i][k]);
            //if( R_max_ind[i] < k ) {
            //    if( m[i][k] != 0 ) printf("PewPew? %1.3f\t%1.3f\t%1.3f\n", grid[k], R_max[i], grid[k+1] );
            //}
        } fprintf(foutput, "\n");
    }
    fclose(foutput);
}

double
max( double *x, int n ) {
	int i; 	double m = x[0]; 
	for( i = 0; i < n; i++ ) {
		if( m < x[i] ) m = x[i];
	}
	return m; 
}

double
min( double *x, int n ) {
	int i; 	double m = x[0]; 
	for( i = 0; i < n; i++ ) {
		if( m > x[i] ) m = x[i];
	}
	return m; 
}

void
range( double *x, int n, double res[2] ) {
	int i;
	double max = x[0], min = x[0]; 
	for( i = 0; i < n; i++ ) {
		if( max < x[i] ) max = x[i];
		if( min > x[i] ) min = x[i];
	}
	res[0] = min; 
	res[1] = max; 
}

/*
void
pdf( double *x, int n ) { 
	int i; 
	i = 
}
*/

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
        if( A[1] == A[0] ) 
            solid_angle = 2*M_PI *  (sin(B[1]) - sin(B[0]));
        else 
	        solid_angle = (A[1] - A[0]) * (sin(B[1]) - sin(B[0]));
    }
	n0 = (double) 3. * N / (solid_angle * (pow(R[1],3) - pow(R[0],3)) );
	printf("%d selected.\n", N);
	printf("A : %+le %+le -> %+le %+le \n", alim[0]*r2d, alim[1]*r2d, A[0]*r2d, A[1]*r2d );
	printf("B : %+le %+le -> %+le %+le \n", blim[0]*r2d, blim[1]*r2d, B[0]*r2d, B[1]*r2d );
	printf("R : %+le %+le -> %+le %+le \n", rlim[0]    , rlim[1]    , R[0]    , R[1]     );
	printf("Solid angle : %le\nDensety : %le\n", solid_angle, n0);
	
    /*
	sprintf(input_name, "%s_sel.dat", name);
	file = fopen(input_name, "w"); 
	for(i=0; i<N; i++) {
		fprintf(file, "%lf %lf %lf %lf %lf %lf\n", p[i][0],  p[i][1], p[i][2], P[i][0],  P[i][1], P[i][2]);
		//fprintf(file, "%lf %lf %lf\n", P[i][0],  P[i][1], P[i][2]);
	}
	fclose(file);
    */
	/*
	sprintf(input_name, "%s_c2.dat", name);
	file = fopen(input_name, "w"); 
	for(i=0; i<N; i++) {
		//fprintf(file, "%lf %lf %lf %lf %lf %lf\n", p[i][0],  p[i][1], p[i][2], P[i][0],  P[i][1], P[i][2]);
		fprintf(file, "%lf %lf %lf\n", p[i][0],  p[i][1], p[i][2]);
	}
	fclose(file);
	*/	
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
init_second() {
	int i,k;
	R_max     = (double *) calloc(N, sizeof(double));
	R_max_ind = (int    *) calloc(N, sizeof(int   ));
	grid = (double *) calloc(grid_n, sizeof(double));
	max_dist = (double *) calloc(grid_n, sizeof(double));
	//mass = (double ***) calloc(N, sizeof(double));
	//for(i=0; i<N; i++) {
	//	mass[i] = (double **) calloc(grid_n, sizeof(double));
	//	for(k=0; k<grid_n; k++) {
	//		mass[i][k] = (double *) calloc(3, sizeof(double)); 
	//	}
	//}
	neighbor = (double *) calloc(N, sizeof(double *));
    for( i = 0; i < N; i++ ) {
        neighbor[i] = R[1]; 
    }
	m = (myint **) calloc(N, sizeof(myint *));
	for(i=0; i<N; i++) {
		m[i] = (myint *) calloc(grid_n, sizeof(myint));
        for( k = 0; k<grid_n; k++ ) {
            m[i][k] = 0;
        }
	}
	/*	
	pos = (int ***) calloc(N, sizeof(int **));
	for(i=0; i<N; i++) {
		pos[i] = (int **) calloc(grid_n, sizeof(int *));
		if(pos[i] == NULL) { printf("I CANT SORRY\n"); exit(42);}
		//for(k=0; k<grid_n; k++) {
		//	pos[i][k] = (int *) calloc(2, sizeof(int));
		//}
	}
	*/
}

void
make_grid() {
	grid_h1 = pow(grid_R[1]/grid_R[0], (double) 1./((double) grid_n - 1));
	grid[0] = grid_R[0]; 
	if(grid[0] == 0) grid[0] = 1e-10;
	int k;
	for( k=1; k<grid_n; k++)	{
		grid[k] = grid[k-1] * grid_h1;
	}	
}

int
index_grid(double r) {
	int i;
	i = (int) (log(r / grid[0] ) / log(grid_h1)); 
	if(i<0) {
        printf("index_grid ALARM %lf\n", r);
        i = 0;
    }
	return i;
}

void
print_summ() {
	printf("%d points\n.", N);
	printf("A:%1.3lf-%1.3lf B:%1.3lf-%1.3lf R:%1.3lf-%1.3lf\n", A[0]*180/M_PI, A[1]*180/M_PI, B[0]*180/M_PI, B[1]*180/M_PI, R[0], R[1]);
	printf("SA: %lf   n0 = %lf\n", solid_angle, n0);
	printf("GRID :: %lf - %lf length = %d\n", grid[0], grid[grid_n - 1], grid_n);
	printf("Output : %s\n", output );
}

void
R_max_find() {
	double r1, r2, r3, r4; 
	double a,b,r;
	double da = 2*M_PI, db = 2*M_PI; 
	int i; 
	
	//printf("%d\n", N);
    if( (A[0] == A[1]) & (B[0] == B[1])) {
        printf("!! Whole sphere !!\n");
        for(i=0; i<N; i++) {
	    	spherical_coordinats(1, &p[i][0], &p[i][1], &p[i][2], &a, &b, &r);
            if( R[0] == 0 ) 
                R_max[i] = (R[1] - r);
            else 
	    	    R_max[i] = MIN(r - R[0], R[1] - r); 
	    	R_max_ind[i] = index_grid(R_max[i]);
        }
    } else {
	    for(i=0; i<N; i++) { 
	    	//i = 4; {
	    	spherical_coordinats(1, &p[i][0], &p[i][1], &p[i][2], &a, &b, &r);
	    	//printf("%lf %lf %lf\n", a*180/M_PI, b*180/M_PI, r);
	    	da = 2*M_PI, db = 2*M_PI;
	    	if( (A[0] != A[1] - 2*M_PI) &
                (A[0] != A[1])  ) {
	    		da = MIN(a - A[0], A[1] - a); 
	    	}
	    	//if( B[0] != -M_PI/2 ) { 
	    		db = b - B[0];
	    	//}	
	    	//if( B[1] != M_PI/2 ) {
	    		db = MIN(db, B[1] - b);
	    	//}	
	    	//printf("da = %lf db = %lf\n", da*180/M_PI, db*180/M_PI);
	    	R_max[i] = MIN(r - R[0], R[1] - r); 
	    	//R_max_ind[i] = 1; 
	    	if( da <= M_PI/2 ) {
	    		R_max[i] = MIN(R_max[i], r * cos(b) * sin( da )) ; 
	    		//if(R_max[i] == r * cos(b) * sin(da)) R_max_ind[i] = 2; 
	    	}
	    	if( db <= M_PI/2 ) {
	    		R_max[i] = MIN(R_max[i], r * sin(db) ); 
	    		//if(R_max[i] == r * sin(db)) R_max_ind[i] = 3; 
	    	}

	    	R_max_ind[i] = index_grid(R_max[i]);
	    	if(R_max[i] == 0) {
	    		printf("Rmax == 0; %lf %lf %lf\n", a, b, r);
	    	}
	    
	    	if(R_max[i] < 0) {
	    		printf("Rmax <0; %lf %lf %lf\n", a, b, r);
	    	}
	    	//printf("%lf %d\n\n", R_max[i], R_max_ind[i]);
	    }	
    }
}

void
R_max_check() {
	int i = 0; 
	double a,b,r,a2,b2,r2;
	double x[100000], y[100000],  z[100000];
	int MISS; int n;
	//FILE *ol, *ol2;
	//ol =fopen("ol.dat", "w");
	///ol2=fopen("ol2.dat", "w");
	for(i=0; i<100; i++) {
		//printf("\ri=%d", i);
	//i = 4; {
		MISS=0; n=0;
		for(b=-M_PI/2+0.001;b<=M_PI/2;b+=M_PI/1000 ) {
    		for(a=0.001;a<=2*M_PI;a+=(2*M_PI)/1000) {
				r = R_max[i]; 
				//printf("\n");
				//printf("%lf %lf %lf\n", a, b, r );
				cartesian_coordinats(1, &(x[n]), &(y[n]), &(z[n]), &a, &b, &r);
				//printf("%lf %lf %lf\n", p[i][0], p[i][1], p[i][2]);
				//printf("%lf %lf %lf\n", x[n], y[n], z[n]);
				x[n] += p[i][0];
				y[n] += p[i][1];
				z[n] += p[i][2];
				//spherical_coordinats(1, &(x[n]), &(y[n]), &(z[n]), &a, &b, &r);
				//fprintf(ol, "%lf %lf %lf\n", x[n], y[n], z[n]);
				//printf("%lf %lf %lf\n", a, b, r );
				//printf("MISS = %d; i = %d %lf %lf %lf %lf %lf %lf %d\n", MISS, i, x[n], y[n], z[n], a*180/M_PI, b*180/M_PI, r, inside_ca(x[n], y[n], z[n]));
				
				if(inside_ca(x[n], y[n], z[n]) != 0) {
					MISS++;
					//spherical_coordinats(1, &(x[n]), &(y[n]), &(z[n]), &a2, &b2, &r2);
					//fprintf(ol2, "%lf %lf %lf\n", x[n], y[n], z[n]);
					//printf("%lf %lf %lf %lf\n", a*180/M_PI, b*180/M_PI, r, R_max[i]);
					//printf("MISS = %d; i = %d %lf %lf %lf  ||  %lf %lf %lf || %d %d\n", MISS, i, a2*180/M_PI, b2*180/M_PI, r2, a*180/M_PI, b*180/M_PI, r, inside_ca(x[n], y[n], z[n]),  inside_sp(a, b, r));
				}
				//n++; 
			}
		}
		if(MISS!=0) {printf("MISS=%d %lf\n", MISS, R_max[i]);}
	}
	//fclose(ol);
	//fclose(ol2);
}

void
hist( double *x, int n, double dr, char *out ) {
	double rg[2]; 
	range(x, n, rg);
	int n_bin = ( rg[1] - rg[0] ) / dr;
	double *h = (double *) calloc( n_bin+1, sizeof(double) );	
	
	int i, j; 
	for( i=0; i<n; i++ ) {
		h[ ( int ) ( (rg[1] - x[i] ) / dr) ] ++ ; 
 	}

	FILE *fout = fopen( out, "w" ); 
	for( i=0; i<n_bin; i++ ) {
		fprintf(fout, "%lf %lf\n", rg[0] + dr*i, h[i]);
	}
	fclose(fout);
}

void
SL() {
    int i, k;
    int max_rmax; 
    for( i=0; i<N; i++ ) { 
        if( max_rmax < R_max_ind[i] ) 
            max_rmax = R_max_ind[i]; 
    }
    printf( "\nSL\n%d\n", max_rmax );

    double scale[4] = { 1., 5., 10., 20. }; 
    int scale_index, scale_n = 4; 
    char filename_sl[80]; 
    for( scale_index = 0; scale_index < scale_n; scale_index++ ) {
	    sprintf(filename_sl, "%s_SL_%03f.dat", name, scale[scale_index]);
        FILE * file = fopen(filename_sl, "w");
        k = index_grid( scale[scale_index] ); 
        for( i=0; i<N; i++ ) {
            if( R_max_ind[i]-1 >= k ) {
                fprintf( file, "%le %d\n", P[i][2], (int) m[i][k] ); 
            }
        }
        fclose(file); 
    }
}

void
ololo() {
    int i,j,k;
    for( i=0; i<N; i++ ) {
        if( grid[R_max_ind[i]] > R_max[i] ) 
            printf("FUCK #1 %d\t%e\t%e\t%e\n", i, grid[R_max_ind[i]],  R_max[i], grid[R_max_ind[i]+1]);
        if( grid[R_max_ind[i]+1] < R_max[i] ) 
            printf("FUCK #2 %d\n", i);
    }
}

void
gogogo() {
	init_begining();
	read_data();
	init_second();
	make_grid();
		//int k = index_grid(10.); 
		//printf("%lf %lf %lf\n", grid[k-1], grid[k], grid[k+1] );
		//exit(0);
	print_summ();
	R_max_find();
	//write_geom();
	//R_max_check();
	count(); 
	write_result(); 
	write_mad(); 
    int k=index_grid(1); 
    printf("%1.3f\t%1.3f\t%1.3f\t%d\n", grid[k-1], grid[k], grid[k+1], index_grid(grid[k]));
	//hist( R_max, N, 1, "rmax_hist.dat" );
	//double foo[10] = {-1,-1,-1,-2,-1,4,0,1,1,1};
	//hist( foo, 10, ( max(foo, 10) - min(foo, 10) ) / 5, "rmax_hist.dat" );
	//printf("\n %lf %lf %lf\n", max(foo, 10) , min(foo, 10), ( max(foo, 10) - min(foo, 10) ) / 5);
	//SL();
    //ololo();
}

int main( int argc, char *argv[] ) {
	char optString[] = {"d:N:g:n:h:a:b:r:q:o:c"};
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
			case 'g':
				sscanf(optarg, "%lf", &grid_R[0]); 
				sscanf(argv[optind], "%lf", &grid_R[1]);
				optind++;
				break;
			case 'n':
				sscanf(optarg, "%d", &grid_n); 
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
				output = optarg; 
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
	if( got_out == 1 ) {
        name = output; 
    }
	char tmp[80]; 
	sprintf(tmp, "%s_gamma.dat", name);
	output = tmp;
	gogogo();
	return(0);	
}

