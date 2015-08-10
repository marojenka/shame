/*
double
calc_nc(int k) {
	int i;
	int nc; 
	nc = 0;
	for(i=0; i<NS; i++) {
		if(R_max_ind[slice_index[i]] >= k - 1)
			nc++; 
	}
	return (double) nc;
}

double
sd(int k) {
	int i;
	int nc; 
	int count;
	count = 0; nc = 0;
	for(i=0; i<N; i++) {
		if(R_max_ind[i] >= k) {
			count += m[i][k]; 	
			nc++; 
		} 	
	}

	double mean = (double) count / nc; 
	double delta = 0;

	for(i=0; i<N; i++) {
		if(R_max_ind[i] >= k) {
			delta += pow( (double) m[i][k] - mean , 2);
		} 	
	}
	return (double) sqrt( (double) 1. / (nc-1) * sqrt(delta) );
}

double
sigma(int k) {
	int i;
	int nc; 
	int count;
	count = 0; nc = 0;
	for(i=0; i<N; i++) {
		if(R_max_ind[i] >= k) {
			count += m[i][k]; 	
			nc++; 
		} 	
	}

	double mean = (double) count / nc; 
	double delta = 0;

	for(i=0; i<N; i++) {
		if(R_max_ind[i] >= k) {
			delta += pow( (double) m[i][k] - mean , 2);
		} 	
	}
	return (double) ( (double) 1. / (nc-1) * delta / pow(mean,2) );
}
*/
/*
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
*/
/*
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
*/

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
    //printf( "\nSL\n%d\n", max_rmax );

    double scale[4] = { 0.1, 1, 10, 60 }; 
    int scale_index, scale_n = 4; 
    char filename_sl[80]; 
    //for( scale_index = 0; scale_index < scale_n; scale_index++ ) {
    for( k = index_grid( 0.1 ); k < index_grid(100) ; k++ ) {
	    sprintf(filename_sl, "SL_%03.f-%03.f_%03f.dat", r[0], r[1], grid[k]);
        FILE * file = fopen(filename_sl, "w");
        //k = index_grid( scale[scale_index] ); 
        for( i=0; i<N; i++ ) {
            if( R_max_ind[i]-1 >= k ) {
                fprintf( file, "%le %llu\n", P[i][2], (int) m[i][k] ); 
            }
        }
        fclose(file); 
    }
}

