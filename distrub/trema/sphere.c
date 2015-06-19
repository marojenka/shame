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
		//S = (double)C/(2.0*(i+500));
		//S = (double)C/(2.0*(i));
		/*
		if( d == 2 ) {
			FOO[i][2] = pow(S / alpha, 1.0/d);
		} else 
		if( d == 3) {
			r = (S / alpha);
			FOO[i][3] = pow( r/4 + r*sqrt(17.0/16.0), 1./3. ) - pow( -r/4 + r*sqrt(17.0/16.0), 1./3. );
		}
		*/
		// radis of the sphere
		//FOO[i][d] = pow(S / alpha / 10, 1.0/d);
	}
}

int
is_it_good(double *point) {
	int i, j;
	double C, S, l, r;

	C = d-D;

	for(i=N1; i<N; i++) {
		S = (double) C/(2.0*(i+1));
		//S = (double) C/(2.0*(i+1));
		l = 0; 
		for(j=0; j<d; j++)
			l += pow(point[j]-FOO[i][j],2); 
		l = sqrt(l);
		r = pow(S / alpha, 1.0/d);
		//printf("%lf %lf\n", l, r);
		//if( l<sqrt(M_PI / S ) )
		if( l < r ) { // l in sphere
			//printf("pew\n");
		//if( fabs(l - r) < (S/(4*M_PI*r)) ) { // l betewin 2 shperas
		//if( fabs(l - r) < FOO[i][d]) { // l betewin 2 shperas
			//printf("%lf %lf %lf %lf\n", sqrt(l), sqrt(S/M_PI), l, S/M_PI );
			return 0; // Not, we don't need it.
		}
	}
	return 42; // YES IT'S PERFECT!
}

int 
limit_N(double r) {
    printf( "%le\n", (double) ( 2. * 4./3. * M_PI * r*r*r ));
    return (int) floor(( d - D ) / ( 2. * 4./3. * M_PI * r*r*r )) ;
}

double 
limit_r(int N) {
    return (double) pow( (d-D)/(2. *4./3.*M_PI) , 1./3.) ;
}

