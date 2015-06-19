RedshiftToDistance <- function(z, h = 1, Omega.m = 0.25, Omega.L = 0.75) {
	r <- 1:length(z)*0
	for( i in 1:length(z)) {
		r[i] <- integrate(function(y) {1 / (y*(Omega.m/y+Omega.L*y^2)^0.5)}, 1/(1+z[i]), 1)$value * 3000
	}
	return(r / h)
}

VisibleToAbsolute <- function( m, r, z, K=0 ) {
	M <- 1:length(m)
	for( i in 1:length(m) ) {
		if( r[i] <= 0 ) {
			M[i] <- 0
			next
		}
		M[i] <- m[i] - 5 * log10(r[i] * (1 + z[i]) ) - 25 
	}
	if( K!=0 ) 
		M <- M - K
	return( M )
}

