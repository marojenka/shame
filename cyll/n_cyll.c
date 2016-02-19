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
int FULL_SPHERE = 0; // is it a full sphere? 

char *input_suff, *name, *prefix,  *output_par, *output_ort; 
int N; // how many points
double **p;  // Array of points. p[i][N]
double **P;  // Array of points. p[i][N]
double *R_max; // Rmax ofc
int    *R_max_ind; // Rmax compared to grid. grid[R_max_ind] < R_max[i] < grid[R_max_ind + 1] ???

int grid_n = 200;
double *grid, grid_h1, grid_R[2] = {1e-3, 1e3}; // grid itseld, h+1 and limits. 
double H = -1; 

typedef unsigned long long int myint; 
myint  **m_par, **m_ort; 

double A[2] = {0, 360}, B[2] = {-90, 90}, R[2] = {0, 300};
double q = 1;
double h = 1;
double r2d = 180/M_PI, d2r = M_PI / 180;
double n0, solid_angle;

int grid_index(double r);
typedef struct statK {
    int Nc; 
    myint c_sh, c_sp; 
    double s2_sh, s2_sp;
} statK;

typedef struct cyl {
    double *R_max;
    int    *R_max_ind; 
    myint **m; 
} cyl; 

cyl cyl_par, cyl_ort;

void
spherical_coordinats(int n, double *x, double *y, double *z, double *a, double *b, double *r) {
	int i; 
	for( i=0; i<n; i++ ) {
		r[i] = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
        a[i] = atan2(y[i], x[i]) ;  
		if( a[i] < 0 ) 
			a[i] += 2*M_PI; 
		if( x[i] == 0 & y[i] == 0 ) {
			if( z[i] >= 0 ) b[i] = +M_PI/2.;
			if( z[i] <  0 ) b[i] = -M_PI/2.;
		} else
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

void
cross_product(double *a, double *b, double *res) {
	res[0] = a[1]*b[2]-a[2]*b[1];
	res[1] = a[2]*b[0]-a[0]*b[2];
	res[2] = a[0]*b[1]-a[1]*b[0];
}

double
scalar_product(double *a, double *b) {
    return( a[0]*b[0] + a[1]*b[1] + a[2]*b[2] );
}

void
count() {
	int i, j, k;
    int index;
	double rd;
    double e[3], eo[3], v[3];
    double r1, r2, r3, r;
    double l, h2, l2;

	for(i=0; i<N; i++) {
        if( inside_sp(P[i][0], P[i][1], P[i][2]) != 0  ) continue; 
	    printf("\ri=%d\t%d%%", i, (int) 100 * i / N );
        if( cyl_par.R_max[i] == -1 ) continue; 
        if( cyl_ort.R_max[i] == -1 ) continue; 

        r = P[i][2];
        if( r == 0 ) {
            e[0] = 0; e[1] = 0; e[2] = 1;
        } else {
            e[0] = p[i][0] / r; 
            e[1] = p[i][1] / r; 
            e[2] = p[i][2] / r; 
        }
        eo[0] = - e[1] / sqrt( e[1]*e[1] + e[0]*e[0]); 
        eo[1] =   e[0] / sqrt( e[1]*e[1] + e[0]*e[0]); 
        eo[2] =   0;
        r1 = r - cyl_par.R_max[i];
        r2 = r + cyl_par.R_max[i];
        r3 = cyl_ort.R_max[i];
            
		for(j=0; j<N; j++) {
            if( i == j ) continue; 

            l  = scalar_product(p[j], e);
            h2 = sqrt( P[j][2]*P[j][2] - l*l );
            if( ( h2 < h ) & (l > r1)&(l < r2) ) {
                if( l > r ) 
                    l = l - r;
                else 
                    l = r - l;
                index = MAX(0, grid_index( l ) + 1 );
                cyl_par.m[i][ index ] ++;
            }
            
            // oort 
            v[0] = p[j][0] - p[i][0]; 
            v[1] = p[j][1] - p[i][1]; 
            v[2] = p[j][2] - p[i][2]; 

            l  = fabs( scalar_product( v, eo ) );
            if( l < r3 )  { 
                h2 = sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] - l*l );
                if( h2 < h ) { 
                    index = MAX(0, grid_index( l ) + 1 );
                    cyl_ort.m[i][ index ] ++;
                }
            }
		}
	}
	printf("\n");
}

/*
double volume_sp(double r) { 
    return( 4./3. * M_PI * pow(r,3)); 
}

double volume_sh(double r1, double r2) { 
    return( 4./3. * M_PI * pow(MAX(r1,r2),3) - pow(MIN(r1,r2),3)); 
}

statK
foo(int K) {
	int i, j, k;
    double sp, sp_mean, sh_mean; 
    statK F; 
	F.s2_sh = 0;
    F.s2_sp = 0;
	F.c_sh = 0;
    F.c_sp = 0;
    F.Nc = 0; 
    // Усы на графике Г (r)  это дсперсия (точнее ср.кв.откл.) значений Г в бинах dr:
    // с.к.в. = (сигма^2)^1/2
    // sigma^2 = (1/(N_c -1)) SUMMA(x_i - <x>)^2
    // а само значение Г - это среднее значение плотности в этих бинах (= <x>).
    // Для самых больших масштабов надо брать не с.к.в. , а относительную величину
    // dN/N = 1/(N)^1/2  где  N - это число центров (а не число точек в шаре), т.к.
    // для максимального шара N_c = 1  и  с.к.в. = бесконечности.  
    // Т.е. надо просто брать в качестве ошибки большее из величин - 
    // либо с.к.в. из дисперсии числа точек в тестовых сферах,
    // либо dN/N = 1/(N)^1/2
	for(i=0; i<N; i++) {
        if( R_max_ind[i] <= K ) 
            continue; 
        F.Nc++; 
        F.c_sh += m[i][K];
        for( k=0; k<=K; k++ ) {
			F.c_sp += m[i][k]; 	
        }
	}
    if( F.Nc == 0 ) {
        F.s2_sp = 0; F.s2_sh = 0;
    } else { 
        sp_mean = (double) F.c_sp / F.Nc;
        sh_mean = (double) F.c_sh / F.Nc;
        for(i=0; i<N; i++) {
            if( R_max_ind[i] <= K ) 
                continue; 
            sp = 0;
            for( k=0; k<=K; k++ ) {
                sp += m[i][k]; 	
            }
            F.s2_sp += pow( sp - sp_mean, 2 ) ;
            F.s2_sh += pow( m[i][K] - sh_mean , 2);
        }
        F.s2_sp = sqrt( (double) 1./(F.Nc - 1) * F.s2_sp  ); 
        F.s2_sh = sqrt( (double) 1./(F.Nc - 1) * F.s2_sh  ); 
    }
	//if( nc == 0 | count == 0 ) return -1; 
	return F;
}
*/

void
write_cyl(cyl *foo, char* name_out) {
	FILE *foutput; 
	int k, i, K; 
    
    double closest = 0;

    printf("write to %s\n", name_out); 
	foutput = fopen(name_out, "w");
    fprintf(foutput, "# %s %5f %5f %5f %5f %5f %5f\n", name, A[0], A[1], B[0], B[1], R[0], R[1]);
    fprintf(foutput, "# %10s %10s %10s\n","N",  "dens", "h");
    fprintf(foutput, "# %10d %10f %2.2f\n", N, n0, h);
    fprintf(foutput, "#\n");

        double min = 3000;
        double max = 0;
        int max_ind; 
        int boo, max_boo = 0, max_boo_ind = 0; 
        for( i=0; i<N; i++ ) {
            if( max < foo->R_max[i] ) {
                max = foo->R_max[i];
                k = i;
            }
            boo = 0;
            for( k=0; k<grid_n; k++ ) {
                boo += foo->m[i][k]; 
            }
            if( max_boo < boo ) {
                max_boo = boo; 
                max_boo_ind = i; 
            }
        }
        max_ind = k;
        printf("i = %d R_max = %lf\n", max_ind, max);
        printf("i = %d boo = %d\n", max_boo_ind, max_boo);


    double mean;
    double mean2;
    int Nc2;
    int Nc = 0; 
    // double count = 0;
    int *count;
    int dcount;
    int sum; 
    double dn;

    count = (int *) calloc(N, sizeof(int));

    // k = 0;
    // for( i=0; i<N; i++ ) {
    //     // if( foo->R_max_ind[i] < k ) continue; 
    //     count[i] = foo->m[i][0];
    // }
    for( k = 0; k<grid_n-1; k++ ) {
        Nc = 0; 
        sum   = 0;
        dcount = 0;
        Nc2 = 0;
        for( i=0; i<N; i++ ) {
            if( foo->R_max[i] <= grid[k] ) continue; 
 
            count[i] += foo->m[i][k]; 
            dcount   += foo->m[i][k];
            sum += count[i];
            Nc++;
            if( foo->m[i][k] != 0 ) Nc2++;
        }
        if( Nc != 0 ) 
            mean  =  (double) sum / Nc;
        else {
            continue;
        }
        if( sum == 0 ) { 
            continue;
        }
        
        // mean = count[0];
        
        fprintf(foutput, "%1.10le ", grid[k]);
        fprintf(foutput, "%1.10le ", mean / (M_PI *2* grid[k] * h*h));
        fprintf(foutput, "%1.10le ", (double) dcount / Nc /(M_PI *2* (grid[k]-grid[k-1]) * h*h) );
        fprintf(foutput, "%f ", 100*((double) Nc/N) );
        fprintf(foutput, "%d ", foo->m[max_boo_ind][k] );
        fprintf(foutput, "%d ", sum );
        fprintf(foutput, "%d ", dcount );
        fprintf(foutput, "%lf ", (double) dcount / Nc );
        fprintf(foutput, "%lf ", (double) Nc2 != 0 ? sum / Nc2 : -1);
        fprintf(foutput, "\n");
    }
    fclose( foutput );
    free(count);

    // foutput = fopen("rmax.dat", "w");
    // for( i=0; i<N; i++ ) {
    //     fprintf(foutput, "%le %le %le %le\n", p[i][0], p[i][1], p[i][2], foo->R_max[i]);   
    // }
    // fclose( foutput );


    //int HAI = index_grid(10.);
    //double HAX = ave_N(HAI) / (4./3.*M_PI*pow(grid[HAI],3));
    // statK F;
	// for(k=1; k<grid_n-1; k++) {
    //     F = foo(k); 
    //     if( F.c_sp == 0 ) continue; 
	// 	fprintf(foutput, "%le ", grid[k] ); 
	// 	fprintf(foutput, "%le ", (double) F.c_sp / F.Nc / volume_sp(grid[k]) ); 
	// 	fprintf(foutput, "%le ", (double) F.c_sh / F.Nc / volume_sh(grid[k], grid[k+1]) ); 
	// 	fprintf(foutput, "%le ", (double) F.s2_sp / volume_sp(grid[k]) );             
	// 	fprintf(foutput, "%le ", (double) F.s2_sh / volume_sh(grid[k], grid[k+1]) );  
	// 	fprintf(foutput, "%d  ", F.c_sp ); 
	// 	fprintf(foutput, "%d  ", F.c_sh ); 
	// 	fprintf(foutput, "%le ", F.s2_sp ); 
	// 	fprintf(foutput, "%le ", F.s2_sh ); 
	// 	fprintf(foutput, "%d  ", F.Nc ); 
	// 	fprintf(foutput, "%d  ", ((F.Nc > N/2.) & (grid[k] > closest)) ? 1 : 0 ); 
	// 	//fprintf(foutput, "%le ", sigma(k) ); 
	// 	//fprintf(foutput, "%le ", calc_nc(k) / N * 100 ); 
	// 	//fprintf(foutput, "%le ", ave_N(k) ); 
	// 	//fprintf(foutput, "%le ", (double) sqrt( 1. / ave_N(k) ) ); 
	// 	//fprintf(foutput, "%le ", grid[k] / closest ); 
    //     //fprintf(foutput, "%le ", max_dist[k]);
	// 	fprintf(foutput, "\n");
	// }
	// fclose(foutput);
}


void
write_rmax() {
    int i;
    FILE *foutput;
    char xyz_rmax[80];

    sprintf(xyz_rmax, "%s_rmax.dat", name );
    foutput = fopen(xyz_rmax, "w");

    for( i=0; i<N; i++ ) {
        if( ((double) rand() / RAND_MAX) < .3 ) continue; 
        fprintf(foutput, "%le %le %le %le %le \n", p[i][0], p[i][1], p[i][2], cyl_par.R_max[i], cyl_ort.R_max[i]);
    }
    fclose( foutput );
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

	solid_angle = (A[1] - A[0]) * (sin(B[1]) - sin(B[0]));
	n0 = (double) 3. * N / (solid_angle * (pow(R[1],3) - pow(R[0],3)) );
	printf("%d selected.\n", N);
	printf("A : % 8.1f % 8.1f -> % 8.1f % 8.1f \n", alim[0]*r2d, alim[1]*r2d, A[0]*r2d, A[1]*r2d );
	printf("B : % 8.1f % 8.1f -> % 8.1f % 8.1f \n", blim[0]*r2d, blim[1]*r2d, B[0]*r2d, B[1]*r2d );
	printf("R : % 8.1f % 8.1f -> % 8.1f % 8.1f \n", rlim[0]    , rlim[1]    , R[0]    , R[1]     );
	printf("Solid angle :\t % 8.1e (%3.2f%%)\n", solid_angle, solid_angle / (4 * M_PI) * 100);
    printf("Densety :\t %8.1le (N=%d V=%8.1e)\n", n0, N, solid_angle / 3 * (pow(R[1],3) - pow(R[0],3)) );
	
    if( q < 1 ) {
        sprintf(input_name, "%s_sel.dat", name);
        file = fopen(input_name, "w"); 
        for(i=0; i<N; i++) {
            fprintf(file, "%lf %lf %lf\n", p[i][0],  p[i][1], p[i][2]);
        }
        fclose(file);
    }
}

void
init_begining() {
    if( A[0] == 0 & A[1] == 360 & 
        B[0] == -90 & B[1] == 90 ) 
        FULL_SPHERE = 1;
	A[0]=A[0]*d2r;  
	A[1]=A[1]*d2r;  
	B[0]=B[0]*d2r;  
	B[1]=B[1]*d2r;  
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

cyl
init_cyl(cyl *foo) {
    int i,k;
	foo->R_max     = (double *) calloc(N, sizeof(double));
	foo->R_max_ind = (int    *) calloc(N, sizeof(int   ));
	foo->m = (myint **) calloc(N, sizeof(myint *));
	for(i=0; i<N; i++) {
		foo->m[i] = (myint *) calloc(grid_n, sizeof(myint));
        for( k = 0; k<grid_n; k++ ) {
            foo->m[i][k] = 0;
        }
	}
}

void
init_second() {
	int i,k;
    init_cyl(&cyl_par);
    init_cyl(&cyl_ort);
}

void
grid_make() {
    int k;

    if( H == -1 ) {
        grid = (double *) calloc(grid_n, sizeof(double));
        grid[0] = grid_R[0]; 
        if(grid[0] == 0) grid[0] = 1e-10;
        grid_h1 = pow(grid_R[1]/grid_R[0], (double) 1./((double) grid_n - 1));
        for( k=1; k<grid_n; k++)	{
            grid[k] = grid[k-1] * grid_h1;
        }	
    } else {
        grid_n = (int) ( grid_R[1] - grid_R[0] ) / H + 1; 
        grid = (double *) calloc(grid_n, sizeof(double));
        grid_h1 = H; 
        grid[0] = grid_R[0]; 
        for( k=1; k<grid_n; k++ ) {
            // grid[k] = grid[0] + grid_h1 * k; 
            grid[k] = grid[k-1] + grid_h1; 
        }
    }
}

int
grid_index(double r) {
	int i;
    if( H == -1 ) {
	    i = (int) (log(r / grid[0] ) / log(grid_h1)); 
    } else {
        // r = grid[0] + i * grid_h1; 
        i = (int) ( r - grid[0] ) / grid_h1; 
    }
	if(i<0) {
        i = 0;
    }
	return i;
}

void
grid_check(double R) {
    printf("Check GRID\n");
    printf("Check GRID\n");
    printf("Check GRID\n");
    printf("GRID: %3.3f %3.3f\n", grid[0], grid[grid_n-1]);
    printf("grid_h1 : %3.3f %d\n", grid_h1, grid_n);
    int kk; 
    kk = grid_index( R );
    printf(" R = %1.3e\n", R);
    printf(" grid[kk-1] : %lf  grid[kk] : %lf  grid[kk+1] : %lf\n", grid[kk-1], grid[kk], grid[kk+1] );
    exit(10);
}

double
rmax(double x, double y, double z) {
	double a,b,r;
	double da = 2*M_PI, db = 2*M_PI; 
    double R_max; 

	spherical_coordinats(1, &x, &y, &z, &a, &b, &r);

    if( inside_sp(a, b, r) != 0 ) {
        return -1;
    }

    if( ((A[0] == A[1]) & (B[0] == B[1])) ||
        ((A[0] == A[1] - 2*M_PI ) & (B[0] == B[1] - M_PI)) )
    {
        // printf("!! Whole sphere !!\n");
            if( R[0] == 0 ) 
                R_max = (R[1] - r);
            else 
	    	    R_max = MIN(r - R[0], R[1] - r); 
    } else {
	    da = 2*M_PI;
        db = 2*M_PI;
	    if( (A[0] != A[1] - 2*M_PI) &
            (A[0] != A[1])  ) {
	    	da = MIN(a - A[0], A[1] - a); 
	    }
	    db = MIN(b - B[0], B[1] - b);
	    R_max = MIN(r - R[0], R[1] - r); 
	    if( da <= M_PI/2 ) {
	    	R_max = MIN(R_max, r * cos(b) * sin( da )) ; 
	    }
	    if( db <= M_PI/2 ) {
	    	R_max = MIN(R_max, r * sin(db) ); 
	    }

	    if(R_max == 0) { printf("Rmax == 0; %lf %lf %lf\n", a, b, r); }
	    if(R_max < 0)  { printf("Rmax <0; %lf %lf %lf\n"  , a, b, r); }
    }
    return(R_max);
}

void
R_max_find_par() {
    double r1, r2;
    double e1, e2, e3;
    double x,y,z,a,b,r;
    double l;
    int i, k, K;
    // double min = 3000;
    // for( i=0; i<N; i++ ) {
    //     if( min > P[i][2] ) {
    //         min = P[i][2];
    //         k = i;
    //     }
    // }
    // printf("%lf\t%d\n", min, k);
    if( FULL_SPHERE ) {
        printf("boundary as in full sphere\n");
    }

	double da = 2*M_PI, db = 2*M_PI; 
    for( i=0; i<N; i++ ) {
    //    r1 = sqrt(R[1]*R[1] - h*h);
    //    r2 = R[0]; 
    //    if( FULL_SPHERE == 0 ) {
    //        // so we have limits in ra and dec
	//		da = 2*M_PI, db = 2*M_PI;
	//		if( A[0] != A[1] - 2*M_PI) {
	//			da = MIN(P[i][0] - A[0], A[1] - P[i][0]); 
	//		}
	//		db = MIN(P[i][1] - B[0], B[1] - P[i][1]); 
	//		if( da < M_PI/2 ) 
	//			r2 = MAX(r2, h/(sin(da)*cos(P[i][1]))) ; 
	//		if( db < M_PI/2 ) 
	//			r2 = MAX(r2, h/tan(db) ); 
	//	
	//		if( (r2>=P[i][2])||(P[i][2]>r1) ) {
	//			cyl_par.R_max[i] = -1; 
    //            cyl_par.R_max_ind[i] = 0;
	//		} else { 
	//			cyl_par.R_max[i]     = MIN(P[i][2] - r2, r1 - P[i][2]); 
    //            cyl_par.R_max_ind[i] = grid_index(cyl_par.R_max[i]);
	//		}
	//	} else {
    //        // if there is no limit in ra and dec
    //        if( r2 == 0 ) 
    //            cyl_par.R_max[i]     = r1 - P[i][2]; 
    //        else 
	//		    cyl_par.R_max[i]     = MIN(P[i][2] - r2, r1 - P[i][2]); 
    //        cyl_par.R_max_ind[i] = grid_index(cyl_par.R_max[i]);
    //    }
    //}
        K = -1;
        for( k=0; k<grid_n; k++ ) {
            r = P[i][2] + grid[k];
            cartesian_coordinats(1, &x, &y, &z, &P[i][0], &P[i][1], &r);
            r1 = rmax(x, y, z);
            //if( r1 < h ) break;
            r = P[i][2] - grid[k];
            a = P[i][0]; b = P[i][1];
            if( r < 0 ) {
                a = a + M_PI; 
                b = - b;
                r = - r;
            }
            cartesian_coordinats(1, &x, &y, &z, &a, &b, &r);
            r2 = rmax(x, y, z);
            if( r1 < h ) break;
            if( r2 < h ) break;
            K = k;
        }
        if( K == -1 ) {
            cyl_par.R_max[i] = -1; 
            K = 0; 
        } else {
            cyl_par.R_max[i] = grid[K];
        }
        cyl_par.R_max_ind[i] = K;
      }
}

void
R_max_find_ort() {
    double r1, r2;
    double e[3], er;
    double x,y,z,a,b,r;
    double l;
    int i, k, K;
    
    for( i=0; i<N; i++ ) {
        K = 0;
        er = sqrt( p[i][1]*p[i][1] + p[i][0]*p[i][0] ); 
        e[0] = -p[i][1] / er;
        e[1] =  p[i][0] / er; 
        e[2] = 0;
        for( k=0; k<grid_n; k++ ) {
            // r = P[i][2] + grid[k];
            // cartesian_coordinats(1, &x, &y, &z, &P[i][0], &P[i][1], &r);
            z = p[i][2];

            x = p[i][0] + e[0] * grid[k]; 
            y = p[i][1] + e[1] * grid[k]; 
            r1 = rmax(x, y, z);
            x = p[i][0] - e[0] * grid[k]; 
            y = p[i][1] - e[1] * grid[k]; 
            r2 = rmax(x, y, z);

            if( r2 < h ) break;
            if( r1 < h ) break;
            K = k;
        }
        cyl_ort.R_max_ind[i] = K;
        cyl_ort.R_max[i] = grid[K];
    }
}

/*
void
R_max_find() {
	double r1, r2, r3, r4; 
	double a,b,r;
	double da = 2*M_PI, db = 2*M_PI; 
	int i; 
	
	//printf("%d\n", N);
    if( ((A[0] == A[1]) & (B[0] == B[1])) ||
        ((A[0] == A[1] - 2*M_PI ) & (B[0] == B[1] - M_PI)) )
    {
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
*/

void
_gogogo() {
	init_begining();
	read_data();
	init_second();
	grid_make();
    // grid_check(10.); 
    // A[0] = 10*M_PI/180;
    // A[1] = 60*M_PI/180;
    // B[0] = 10*M_PI/180;
    // B[1] = 60*M_PI/180;
    // R[0] = 10;
    // R[1] = 250;
    
    ////
	R_max_find_par();
	R_max_find_ort();
	count(); 
    write_cyl(&cyl_par, output_par);
    write_cyl(&cyl_ort, output_ort);

    FILE *F1, *F2, *F3, *F4;
    F1 = fopen("count_par.dat", "w");
    F2 = fopen("count_ort.dat", "w");
    F3 = fopen("rmax.dat", "w");
    F4 = fopen("grid.dat", "w");
    int i, k; 
    for(i=0; i<N; i++) {
        for( k=0; k<grid_n; k++ ) {
            fprintf( F1, " %d ", cyl_par.m[i][k] );
            fprintf( F2, " %d ", cyl_ort.m[i][k] );
        }
        fprintf( F1, "\n" );
        fprintf( F2, "\n" );
    }
    for(i=0; i<N; i++) {
        fprintf( F3, "%lf %lf %lf ", p[i][0], p[i][1], p[i][2] ); 
        fprintf( F3, "%lf ", cyl_ort.R_max[i] );
        fprintf( F3, "%lf ", cyl_par.R_max[i] );
        fprintf( F3, "\n");
    }
    for( k=0; k<grid_n; k++ ) {
        fprintf( F4, "%lf \n", grid[k] );
    }
    fclose(F1);
    fclose(F2);
    fclose(F3);
    fclose(F4);

    // write_rmax();
}

int main( int argc, char *argv[] ) {
	char optString[] = {"d:N:g:n:h:a:b:r:q:o:cH:"};
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
				break;
			case 'o':
				got_out = 1;
				prefix = optarg; 
				break;
            case 'c': 
                SELECT_COORDINATES = 1; // spherical
                printf("I'll work with spherical coordinates.\n");
                break;
			case 'h': 
				sscanf(optarg, "%lf", &h); 
				break;
			case 'H': 
				sscanf(optarg, "%lf", &H); 
				break;


			default:
				printf("What? %d\n", opt);
				break;
    	}
	    opt = getopt( argc, argv, optString );
	}
	if( got_out == 1 ) {
        name = prefix; 
    }
	char tmp1[80]; 
	sprintf(tmp1, "%s_n_cyl_par.dat", name);
	output_par = tmp1;
	char tmp2[80]; 
	sprintf(tmp2, "%s_n_cyl_ort.dat", name);
	output_ort = tmp2;
    _gogogo();
	return(0);	
}

