
void
void_centers() {
	int i, j; 
	double h, C, S, r; 
	C = d-D;
	
	FOO = (double **) calloc( N, sizeof( double* ) );
	for(i = 0; i<N; i++) {
		FOO[i] = (double *) calloc( d+2, sizeof( double ) );
		for(j=0; j<d; j++) {
			FOO[i][j] = ((double) rand() / RAND_MAX);
		}
		S = (double)C/(2.0*(i+1));
		// radis of the cube
		// mmm. Area between r - 2r 
		FOO[i][d] = pow(S/(pow(2,d)-1) , 1.0/d) / 2;
		//printf("FOO[i][d]=%lf\n", FOO[i][d]);
	}
}

int
is_it_good(double *point) {
	int i, j;
	double C, S, l, r;

	for(i=0; i<N; i++) {
		l = 0;
		for(j=0; j<d; j++) {
			if ( fabs( point[j]-FOO[i][j] ) > 2 * FOO[i][d] ) { 
				l = 1;
			}
		}
		if( l == 0 ) {
			for(j=0; j<d; j++) {
				if ( fabs( point[j]-FOO[i][j] ) >     FOO[i][d] ) { 
					l = 1;
				}
			}
			if( l == 1 ) {
				return 0;
			}
		}
		

	}
	return 42; // YES IT'S PERFECT!
}


