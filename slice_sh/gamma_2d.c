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
int N, NS; // how many points
double *R_max; // Rmax ofc
int    *R_max_ind; // Rmax compared to grid. grid[R_max_ind] < R_max[i] < grid[R_max_ind + 1] ???
double **p, **P;  // Array of points. p[i][N] and P[i][N]
double **s, **S;  // Array of points. slice part 
int    *slice_index; // connection of indexes slice -- main set.
double R0 = 50 , dr = 10, r[2]; 

int grid_n = 100;
double *grid, grid_h1, grid_R[2] = {0.001, M_PI/2}; // grid itself, h+1 and limits. 

typedef unsigned long long int myint; 
myint  **m; 
myint  max = 1;
double *neighbor; 

// double A[2] = {0, 90}, B[2] = {0, 90}, R[2] = {1, 300};
double A[2] = {0, 360}, B[2] = {-90, 90}, R[2] = {0, 300};
int _NO_RA;
double q =1;
double r2d = 180/M_PI, d2r = M_PI / 180;
double n0, solid_angle, density;

double s_R_max( int i); 
int    s_R_max_ind( int i);

char *limits_file = "limits.txt";
double slice_limits[100][2];
int slice_limits_N;
//double step[];; 

int index_grid(double r); 

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
add(int i, int j, int K) {
	int k;
    for(k=K+1; k<=MIN( s_R_max_ind(i)+1, grid_n); k++) {
    //for(k=K+1; k<=(grid_n); k++) {
        m[i][k] ++;
    }
}

void
count() {
	int i, j;
	double rd;
    FILE *file; 
	for(i=0; i<NS; i++) {
	    printf("\r%d%%", (int) 100 * i / NS );
		for(j=i+1; j<NS; j++) {
			rd = acos(
                    sin( S[i][1] ) * sin( S[j][1] ) + 
                    cos( S[i][1] ) * cos( S[j][1] ) * cos( S[i][0] - S[j][0] ) ); 
			if( rd < neighbor[i] ) neighbor[i] = rd;
			if( rd < neighbor[j] ) neighbor[j] = rd;
			if( rd < s_R_max(i) ) 
				add(i, j, index_grid(rd)); 
			if( rd < s_R_max(j) )
				add(j, i, index_grid(rd)); 
		}
	}
    printf("\n");
}

double
ave_N(int k) {
	int i, j;
	int nc; 
	unsigned long long int count;
	count = 0; nc = 0;
	for(i=0; i<NS; i++) {
		if( s_R_max_ind(i)-1 >= k) {
			count += m[i][k]; 	
			nc++; 
            //printf("\r%llu",count);
		} 	
	}
    //printf("\n");
	if( nc == 0 | count == 0 ) return -1; 
	return (double) count/nc;
}

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
    double    * sigma  = (double *) calloc( grid_n, sizeof(double) );
    double    * variance  = (double *) calloc( grid_n, sizeof(double) );

    double closest = 0;
    for( i=0; i<NS; i++ ) {
        closest += neighbor[i]; 
    }
    closest /= NS; 
    

    for( k=0; k<grid_n; k++ ) {
        nc[k] = 0; count[k] = 0; mean[k] = 0;
        zer0[k] = 0; 
        for(i=0; i<NS; i++) {
            // if( s_R_max(i) < 0.3 ) continue; 
            if( s_R_max_ind(i)-1 >= k) {
                count[k] += m[i][k]; 	
                if( m[i][k] < 0 ) printf("M < 0. OPPA \n");
                if( count[k] < 0 ) printf(" count too small HAHA \n");
                nc[k]++; 
                if( m[i][k] > 0 ) zer0[k]++;
                if( m[i][k] <= 1 ) zer1[k]++;
                if( m[i][k] <= 2 ) zer2[k]++;
            } 	
        }
        mean[k] = (nc[k] != 0) ? (double) count[k] / nc[k]  :  -1;
        //if( mean[k] < 0 & grid[k] < 0.38 ) printf("k = %d, r = %lf mean = %lf\n", k, grid[k], mean[k] ); 
        if( mean[k] <= 0 ) continue; //printf("k = %d, r = %lf mean = %lf\n", k, grid[k], mean[k] );
        if( mean[k] != 0 & count[k] == 0 ) printf(" ddddd \n"); 
        tmp = 0;
        for( i=0; i<NS; i++ ) {
            if( s_R_max_ind(i)-1 >= k ) {
                tmp += pow( m[i][k] - mean[k] , 2);
            }
        }
        //if( tmp == 0 ) printf("oops\n"); 
        delta[k] = sqrt( (double) tmp / (nc[k]-1) / nc[k] );
        sigma[k] = sqrt( (double) tmp / (nc[k]-1) ) / mean[k] ;
        tmp = 0;
        for( i=0; i<NS; i++ ) {
            if( s_R_max_ind(i)-1 >= k ) {
                tmp += pow( m[i][k] , 2);
            }
        }
        variance[k] = (double) tmp / nc[k] - pow( mean[k], 2 );
        // nc[k] = (double) nc[k] / NS;
        zer0[k] = (double) zer0[k] / NS * 100; 
        zer1[k] = (double) zer1[k] / NS * 100; 
        zer2[k] = (double) zer2[k] / NS * 100; 
    }

	for(k=2; k<grid_n; k++) {
          if( nc[k] <= 0 ) continue; 
          if( mean[k] == 0 ) continue; 
          if( mean[k-1] == 0 ) continue; 
		fprintf(foutput, "%le ", R0 * grid[k] ); 
		fprintf(foutput, "%le ", mean[k]  / ( 2*M_PI * ( 1 - cos(grid[k]) ) * (pow(r[1],3) - pow(r[0],3) ))*3 );
		fprintf(foutput, "%le ", (mean[k] - mean[k-1])  / ( 2*M_PI * ( 1 - cos(grid[k  ]) ) - 
                                                            2*M_PI * ( 1 - cos(grid[k-1]) ) 
                                                           ) * 3 / (pow(r[1],3) - pow(r[0],3) ) );
		fprintf(foutput, "%le ", grid[k] ); 
		fprintf(foutput, "%le ", delta[k] / ( 2*M_PI * ( 1 - cos(grid[k]) ) ) );
        fprintf(foutput, "%le ", sigma[k] );
        fprintf(foutput, "%le ", variance[k] / pow( 2*M_PI * ( 1 - cos(grid[k]) )  ,2) );
		fprintf(foutput, "%le ", R0 * 2 * sin( grid[k] / 2 ) ); 
		fprintf(foutput, "%le ", nc   [k]  ); 
		fprintf(foutput, "%le ", zer0 [k]  ); 
		fprintf(foutput, "\n");
	}


    free(mean );
    free(count);
    free(nc   );
    free(delta);
    free(sigma);
    free(variance);
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
        fprintf(foo_file, "%f %f %d\n", S[i][0], S[i][1], (long int)  m[i][ scale_ind ]);    
    }
    fclose( foo_file );
}

void
write_results() {
	FILE *foutput; 
	int k, i; 

    double closest = 0;
    for( i=0; i<NS; i++ ) {
        closest += neighbor[i]; 
    }
    closest /= NS; 
    
    char output[80]; 
	sprintf(output, "%s_%03.f-%03.f.dat", out_suff, r[0], r[1]);
    printf("write to %s\n", output); 

	foutput = fopen(output, "w");
    fprintf(foutput, "# %s %5f %5f %5f %5f %5f %5f\n", name, A[0], A[1], B[0], B[1], R[0], R[1]);
    fprintf(foutput, "# %10s %10s %10s %10s %10s %10s\n", "r1", "r2", "N", "dens", "neig_ang", "neig_dsit");
    fprintf(foutput, "# %10.1f %10.1f %10d %10.2f %10.2f %10.2f\n", r[0], r[1], NS, density, closest, R0 * 2 * sin(closest/2));
    fprintf(foutput, "#\n");
    //foutput = fopen("gamma_2d.dat", "w");

    write_gamma( foutput ); 

	fclose(foutput);
    
    // write_foo( 0.01 );
    // write_foo( 0.05 );
    // write_foo( 0.1 );
    // write_foo( 0.3 );
    // double l; 
    // double stepts[] = {0.001, 0.005, 0.01};
    // for( i = 0; i < 3; i += 1 ) {
    //     write_foo( stepts[i] ) ;
    // }

    // for( l=0.02; l<0.4; l+=0.02 ) {
    //     write_foo( l ) ;
    //     //if( l == 0.01 ) l = 0.0; 
    // }

	//foutput = fopen("foo.dat", "w");
    //for( i=0; i<NS; i++ ) {
    //    fprintf( foutput, "%le %d\n", R_max[slice_index[i]], R_max_ind[slice_index[i]] );
    //}
	//fclose(foutput);

	//foutput = fopen("foo.dat", "w");
    //for( k=0; k<grid_n; k++ ) {
    //    fprintf( foutput, "%le ", grid[k] );
    //    for( i=0; i<10; i++ ) {
    //        fprintf( foutput, "%llu ", m[i][k]  );
    //    }
    //    fprintf( foutput, "\n" );
    //}
	//fclose(foutput);

	//sprintf(output, "slice_%03.f-%03.f.dat", r[0], r[1]);
	//foutput = fopen(output, "w");
    //fprintf(foutput, "0 0 0 0 0 0 0 ");
    //for( k=0; k<grid_n; k++ ) {
    //    fprintf( foutput, "%le ", grid[k] );
    //}
    //fprintf( foutput, "\n" ); 
    //for( i=0; i<N; i++ ) {
    //    fprintf( foutput, "%le %le %le %le %le %le %le", s[i][0], s[i][1],  s[i][2], S[i][0], S[i][1], S[i][2], R_max[i] );
    //    for( k=0; k<grid_n; k++ ) {
    //        if( R_max_ind[i]-1 >= k ) 
    //            fprintf(foutput, " %llu ", m[i][k]);
    //        else 
    //            fprintf(foutput, " NaN ");
    //    }
    //    fprintf(foutput, "\n");
    //}
	//fclose(foutput);
}

void
make_grid() {
    grid[0] = grid_R[0]; 
	if(grid[0] == 0) grid[0] = 1e-10;
	grid_h1 = pow(grid_R[1]/grid[0], (double) 1./((double) grid_n - 1));
	int k;
    printf("GRID: %lf %lf %lf\n", grid_R[0], grid_R[1], grid_h1);
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

double 
s_R_max( int i ) {
    return( R_max[ slice_index[ i ] ] ); 
}

int 
s_R_max_ind( int i ) {
    return( R_max_ind[ slice_index[ i ] ] );
}

int
R_max_find( int i ) {
    int j, k;
    int length = 1000;
    double alpha, delta; 
    double rd;

    double l1, l2, l3, l4; 

    if( R_max[i] != 0 ) return 0;

    l1 = asin( cos(P[i][1]) * sin( MIN(M_PI/2 , P[i][0] - A[0] ) ) );
    l2 = asin( cos(P[i][1]) * sin( MIN(M_PI/2 , A[1] - P[i][0] ) ) );
    l3 = P[i][1] - B[0];
    l4 = B[1] - P[i][1];               

    if( (A[0] == 0 ) & ( A[1] == M_PI*2 ) ){ 
        l1 = l2 = M_PI; 
        _NO_RA = 1;
    }
    if( B[0] == - M_PI / 2 ) {
        l3 = M_PI; 
    }
    if( B[1] == M_PI / 2 ) {
        l4 = M_PI; 
    }
   
    //  
    // Dumb way to estimate boundary 
    //  R_max[i] = M_PI; 
    //  for( j = 0; j < length; j ++ ) {
    //      for( k = 0; k < 4; k++ ) {
    //          switch( k ) {
    //              case 0 : alpha = A[0]; 
    //                       delta = B[0] + (double) j * ( B[1] - B[0] ) / length; 
    //                      break;
    //              case 1 : alpha = A[1]; 
    //                       delta = B[0] + (double) j * ( B[1] - B[0] ) / length;
    //                      break;
    //              case 2 : alpha = A[0] + (double) j * ( A[1] - A[0] ) / length;
    //                       delta = B[0]; 
    //                      break;
    //              case 3 : alpha = A[0] + (double) j * ( A[1] - A[0] ) / length; 
    //                       delta = B[1]; 
    //                      break;
    //          }
    //          rd = acos( sin( P[i][1] ) * sin( delta ) +
    //                     cos( P[i][1] ) * cos( delta ) * cos( P[i][0] - alpha ) );
    //          if( rd < R_max[i] ) {
    //              R_max[i] = rd;
    //          }
    //      }
    //  }
    // 
    // 
    
    R_max[i] = MIN( M_PI, l1 );
    R_max[i] = MIN( R_max[i], l2 );
    R_max[i] = MIN( R_max[i], l3 );
    R_max[i] = MIN( R_max[i], l4 );
    R_max_ind[i] = index_grid( R_max[i] ); 
}

void
R_max_find_test() {
	double r1, r2, r3, r4; 
	int i; 
	
	for(i=0; i<N; i++) { 
        r1 = R0 * acos( pow( sin(S[i][1]), 2) + pow( cos(S[i][1]), 2) * cos( S[i][0] - A[0] ) );
        r2 = R0 * acos( pow( sin(S[i][1]), 2) + pow( cos(S[i][1]), 2) * cos( S[i][0] - A[1] ) );
        r3 = R0 * (S[i][1] - B[0]); 
        r4 = R0 * (B[1] - S[i][1]);
	    
        R_max[i] = MIN( r1, r2 ); 
        R_max[i] = MIN( R_max[i], r3 );
        R_max[i] = MIN( R_max[i], r4 );
		if(R_max[i] <= 0) {
			printf("Rmax <0\n");
		}
        R_max_ind[i] = index_grid( R_max[i] ); 
	}	

    double alpha, delta, a, d; 
    FILE *file1, *file2, *file;
    char filename[80];
    int j, k, length = 10000, bankai;
    double rd;
    
    // write border
    sprintf(filename, "border.dat");
    file = fopen( filename, "w" ); 
    for( j = 0; j < length; j ++ ) {
        for( k = 0; k < 4; k++ ) {
            switch( k ) {
            case 0 : alpha = A[0]; 
                     delta = B[0] + (double) j * ( B[1] - B[0] ) / length; 
                    break;
            case 1 : alpha = A[1]; 
                     delta = B[0] + (double) j * ( B[1] - B[0] ) / length;
                    break;
            case 2 : alpha = A[0] + (double) j * ( A[1] - A[0] ) / length;
                     delta = B[0]; 
                    break;
            case 3 : alpha = A[0] + (double) j * ( A[1] - A[0] ) / length; 
                     delta = B[1]; 
                    break;
            }
            fprintf( file, "%le %le\n", alpha, delta ); 
        }
    }
    fclose( file ); 

    for( i=0; i<10; i++ ) {
        bankai = 0;
        sprintf(filename, "border_%d.dat", i);
        file = fopen( filename, "w" ); 
        //for( a = A[0]; a < A[1]; a += (A[1] - A[0]) / 100 ) {
        //    for( d = B[0]; d < B[1]; d += (B[1] - B[0]) / 100 ) {
        for( j = 0; j < length; j ++ ) {
            for( k = 0; k < 4; k++ ) {
                switch( k ) {
                    case 0 : alpha = A[0]; 
                             delta = B[0] + (double) j * ( B[1] - B[0] ) / length; 
                            break;
                    case 1 : alpha = A[1]; 
                             delta = B[0] + (double) j * ( B[1] - B[0] ) / length;
                            break;
                    case 2 : alpha = A[0] + (double) j * ( A[1] - A[0] ) / length;
                             delta = B[0]; 
                            break;
                    case 3 : alpha = A[0] + (double) j * ( A[1] - A[0] ) / length; 
                             delta = B[1]; 
                            break;
                }

                rd = R0 * acos( sin( S[i][1] ) * sin( delta ) +
                                cos( S[i][1] ) * cos( delta ) * cos( S[i][0] - alpha ) );
                if( rd < R_max[i] ) {
                    //R_max[i] = rd;
                    bankai ++; 
                    fprintf( file, "%le %le %le\n", alpha, delta, rd - R_max[i] ); 
                }
            }
        }
        //R_max_ind[i] = index_grid( R_max[i] ); 
        if( bankai != 0 ) {
            printf( "bainkai! : %d %d\n", i, bankai );
        }
        fclose( file ); 
        
        sprintf( filename, "circle_%d.dat", i); 
        file = fopen( filename, "w" );
        fprintf( file, "%le %le\n", S[i][0], S[i][1] );
        for( j = 0; j < length; j++ ) {
            delta = B[0] + (double) j * ( B[1] - B[0] ) / length; 
            alpha = ( ( cos( R_max[i] / R0 ) - sin( S[i][1] ) * sin( delta ) ) / ( cos( S[i][1] ) * cos( delta ) )  ); 
            if( fabs( alpha ) > 1 ) {
                continue; 
            } else { 
                alpha = acos( alpha );
            }
            fprintf( file, "%le %le\n%le %le\n", S[i][0] + alpha, delta, S[i][0] - alpha, delta); 
        }
        fclose( file ); 
    }
}

void
R_max_check() {
        //int i = 0; 
        //double a,b,r,a2,b2,r2;
        //double x[100000], y[100000],  z[100000];
        //int MISS; int n;
        //FILE *ol, *ol2;
        //ol =fopen("ol.dat", "w");
        /////ol2=fopen("ol2.dat", "w");
        ////for(i=0; i<100; i++) {
        //i = 4; {
        //    MISS=0; n=0;
        //	for(b=-M_PI/2+0.001;b<=M_PI/2;b+=M_PI/100 ) {
        //		for(a=0.001;a<=2*M_PI;a+=(2*M_PI)/100) {
        //			r = R_max[i]; 
        //			//printf("\n");
        //			//printf("%lf %lf %lf\n", a, b, r );
        //			cartesian_coordinats(1, &(x[n]), &(y[n]), &(z[n]), &a, &b, &r);
        //			//printf("%lf %lf %lf\n", p[i][0], p[i][1], p[i][2]);
        //			//printf("%lf %lf %lf\n", x[n], y[n], z[n]);
        //			x[n] += s[i][0];
        //			y[n] += s[i][1];
        //			z[n] += s[i][2];
        //			spherical_coordinats(1, &(x[n]), &(y[n]), &(z[n]), &a2, &b2, &r2);
        //			fprintf(ol, "%lf %lf %lf %lf %lf %lf\n", x[n], y[n], z[n], a2, b2, r2);
        //			//printf("%lf %lf %lf\n", a, b, r );
        //			//printf("MISS = %d; i = %d %lf %lf %lf %lf %lf %lf %d\n", MISS, i, x[n], y[n], z[n], a*180/M_PI, b*180/M_PI, r, inside_ca(x[n], y[n], z[n]));
        //			
        //			//if(inside_ca(x[n], y[n], z[n]) != 0) {
        //			//	MISS++;
        //			//	//spherical_coordinats(1, &(x[n]), &(y[n]), &(z[n]), &a2, &b2, &r2);
        //			//	//fprintf(ol2, "%lf %lf %lf\n", x[n], y[n], z[n]);
        //			//	//printf("%lf %lf %lf %lf\n", a*180/M_PI, b*180/M_PI, r, R_max[i]);
        //			//	//printf("MISS = %d; i = %d %lf %lf %lf  ||  %lf %lf %lf || %d %d\n", MISS, i, a2*180/M_PI, b2*180/M_PI, r2, a*180/M_PI, b*180/M_PI, r, inside_ca(x[n], y[n], z[n]),  inside_sp(a, b, r));
        //			//}
        //			//n++; 
        //		}
        //	}
        //	///if(MISS!=0) {printf("MISS=%d %lf\n", MISS, R_max[i]);}
        //}
        //fclose(ol);
        ////fclose(ol2);
    double a[2], b[2], r[2]; 
    a[0] = A[0] + (A[1] - A[0] ) / 6; 
    a[1] = A[1] - (A[1] - A[0] ) / 6; 
    b[0] = B[0] + (B[1] - B[0] ) / 6; 
    b[1] = B[1] - (B[1] - B[0] ) / 6; 
    int bingo = 0;
    int i, j, k = 0;
    double r1, r2, r3, r4;
    double rd, rmax;
    for( i=0; i<N; i++ ) {
        k = 0;
        bingo = 0;
        if( ( S[i][0] < a[0] ) || ( S[i][0] > a[1] ) ||
            ( S[i][1] < b[0] ) || ( S[i][1] > b[1] ) )
            continue;    
        r1 = R0 * acos( pow( sin(S[i][1]), 2) + pow( cos(S[i][1]), 2) * cos( S[i][0] - a[0] ) );
        r2 = R0 * acos( pow( sin(S[i][1]), 2) + pow( cos(S[i][1]), 2) * cos( S[i][0] - a[1] ) );
        r3 = R0 * (S[i][1] - b[0]); 
        r4 = R0 * (b[1] - S[i][1]);
	    
        rmax = MIN( r1, r2 ); 
        rmax = MIN( rmax, r3 );
        rmax = MIN( rmax, r4 );
        if( rmax == r1 ) k |= 1 ; 
        if( rmax == r2 ) k |= 2 ; 
        if( rmax == r3 ) k |= 4 ; 
        if( rmax == r4 ) k |= 8 ; 
        for( j = i+1; j<N; j++ ) {
            rd = R0 * acos( sin( S[i][1] ) * sin( S[j][1] ) + 
                            cos( S[i][1] ) * cos( S[j][1] ) * cos( S[i][0] - S[j][0] ) ); 
            if( rd < rmax ) {
                if( ( S[j][0] < a[0] ) || ( S[j][0] > a[1] ) ||
                    ( S[j][1] < b[0] ) || ( S[j][1] > b[1] ) ) {
                    bingo ++; 
                    printf("\t%le %le %le\n", rd-rmax, rd, rmax);
                }
            }
        }
        if( bingo != 0 ) {
            printf( "pew :: %d %d %d\n", i, bingo, k );
        }
    }
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
init_second() {
	int i,k;
    R_max     = (double *) calloc(N, sizeof(double));
    R_max_ind = (int    *) calloc(N, sizeof(int   ));
    grid = (double *) calloc(grid_n, sizeof(double));

	neighbor = (double *) calloc(N, sizeof(double *));
    for( i = 0; i < N; i++ ) {
        neighbor[i] = R[1]; 
    }
	m = (unsigned long long int **) calloc(N, sizeof(unsigned long long int *));
	for(i=0; i<N; i++) {
		m[i] = (unsigned long long int *) calloc(grid_n, sizeof(unsigned long long int));
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
            R_max_find(i); 
            slice_index[j] = i; 
            j++;
        }
    }
    j--;
    density = (double) j / (solid_angle);
    NS = j;
    printf("slice %e %e -- %d\n", r[0], r[1], NS);
    char slice_file[80];
    sprintf(slice_file, "slice_%3e_%3e.dat", r[0], r[1]);
	FILE* file = fopen(slice_file, "w"); 
	for(i=0; i<NS; i++) {
		fprintf(file, "%lf %lf %lf %lf %lf %lf %lf\n", s[i][0],  s[i][1], s[i][2], S[i][0],  S[i][1], S[i][2], s_R_max(i));
	}
	fclose(file);
}

void
whipe() {
    int i, j, k;
    for( i=0; i<N; i++ ) {
        for( k=0; k<grid_n; k++ ) {
            m[i][k] = 0;
        }
        neighbor[i] = R[1]; 
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
        count();
        write_results();
        //SL();
        whipe();
        N = old_N;
    }
}

int main( int argc, char *argv[] ) {
	char optString[] = {"d:N:q:a:b:r:n:R:co:"};
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
    // if( verbose ) { 
    //     if( _NO_RA ) prntf("No RA in bounders\n");
    // }
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
