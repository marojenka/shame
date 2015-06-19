source('~/work/code/R/functions/astronomy.R')
source('~/work/code/R/functions/coordinates.R')
source('~/work/code/R/functions/bounders.R')

GetRandomSphere <- function(N=100, center=c(0,0,0), R=1) {
	cartesian <- matrix(ncol=3, nrow=N, dimnames=list(1:N, c("x", "y", "z")) )
	spherical <- matrix(ncol=3, nrow=N, dimnames=list(1:N, c("a", "b", "r")) )
	
	spherical  <- cbind(
			runif(N, 0, 2*pi)     ,
			asin(runif(N, -1, 1 )),
			(runif(N, 0, R^3))^(1/3)
	)
	cartesian <- cbind(
			cos(spherical[,1]) * cos(spherical[,2]) * spherical[,3] + center[1],
			sin(spherical[,1]) * cos(spherical[,2]) * spherical[,3] + center[2],
								 sin(spherical[,2]) * spherical[,3] + center[3]
	)
	spherical[,1:2] <- spherical[,1:2] * 180 / pi

	for(i in 1:N) {
		spherical[i,] <- cbind( 
						atan2( cartesian[i,2], cartesian[i,1] ) * 180 / pi, # [-180  180] 
						90 - atan2( sqrt(cartesian[i,2]^2 + cartesian[i,1]^2), cartesian[i,3] ) * 180 / pi,  
						sqrt(sum(cartesian[i,]**2)) 
		)
	}
	index <- spherical[,1] < 0
	spherical[index, 1] <- 360 + spherical[index, 1]

	random_sphere <- list(
			spherical = spherical ,
			cartesian = cartesian ,
			description = "A sphere of random points."
	)
	return(random_sphere)
}

GetRandomSlice <- function(N=100, R=c(0,1), A=c(0,360), B=c(-90, 90)) {
	cartesian <- matrix(ncol=3, nrow=N, dimnames=list(1:N, c("x", "y", "z")) )
	spherical <- matrix(ncol=3, nrow=N, dimnames=list(1:N, c("a", "b", "r")) )
	
	A <- A*pi/180
	B <- B*pi/180

	spherical  <- cbind(
			runif(N, A[1], A[2])      ,
			asin(runif(N, sin(B[1]), sin(B[2]) )), 
			runif(N, R[1]^3, R[2]^3)^(1/3)
	)
	cartesian <- cbind(
			cos(spherical[,1]) * cos(spherical[,2]) * spherical[,3],
			sin(spherical[,1]) * cos(spherical[,2]) * spherical[,3],
								 sin(spherical[,2]) * spherical[,3]
	)
	spherical[,1:2] <- spherical[,1:2] * 180 / pi
	
	rand <- cbind(spherical, cartesian)
	colnames(rand) <- c("a", "b", "r", "x", "y", "z")
	return(rand)
}

make.VL.index <- function(DATA, VL.limits) {
	index <-  DATA[, 4] > VL.limits[1] &
			  DATA[, 4] < VL.limits[2] &
			  DATA[, 5] > VL.limits[3] &
			  DATA[, 5] < VL.limits[4] 
	return(index)
}

write.VL <- function(CATALOG, i=1, PREFIX="", make_random=1) {
	index <- make.VL.index(CATALOG$data, CATALOG$VL_lim[i,])
	tmp.abrxyz <- CoordinatesSphericalToDecart( a = CATALOG$data[index, 1], 
			                                    b = CATALOG$data[index, 2], 
												r = CATALOG$data[index, 4] )

	Rsp <- bounders( a = tmp.abrxyz[,1], 
			         b = tmp.abrxyz[,2], 
					 r = tmp.abrxyz[,3],
					 limit = c(CATALOG$VL_lim[i,1:2], CATALOG$area) )
	
	write.table(tmp.abrxyz[,4:6], file=paste0(PREFIX, "VL", i,      '.dat'), row.names=F, col.names=F)
	write.table(Rsp             , file=paste0(PREFIX, "VL", i, '_Rmax.dat'), row.names=F, col.names=F)
	
	if(make_random!=0) {
		write.mock.VL(CATALOG, i, PREFIX, k=make_random)
	}
}

write.mock.VL <- function(CATALOG, i=1, PREFIX="", N=0, k=1) {
	index <- make.VL.index(CATALOG$data, CATALOG$VL_lim[i,])
	count <- max(N, k*sum(index))
	foo <- GetRandomSlice( N = count, 
			               R = CATALOG$VL_lim[i,1:2], 
						   A = CATALOG$area[1:2], 
						   B = CATALOG$area[3:4] )

	foo.Rsp <- bounders( a = foo[,1], b = foo[,2], r = foo[,3], limit = c(CATALOG$VL_lim[i,1:2], CATALOG$area) )
	write.table(foo[,4:6], file=paste0(PREFIX, "RN", i,      '.dat'), row.names=F, col.names=F)
	write.table(foo.Rsp  , file=paste0(PREFIX, "RN", i, '_Rmax.dat'), row.names=F, col.names=F)
}

volume <- function( volume_limits ) {
	r <- volume_limits[1:2]
	a <- volume_limits[3:4]
	b <- volume_limits[5:6]
	
	V <- 1/3 * (r[2]^3 - r[1]^3) * ( a[2] - a[1] ) * (sin(b[2]) - sin(b[1]))
	return(V)
}

mock <- function(r, limits, br=20) {
	foo <- hist(r, br=br, plot=FALSE)	
	rand <- matrix(ncol=6, nrow=length(r))
	j <- 0
	for(i in 1:length(foo$counts) ) {
		if( foo$counts[i] == 0) next 
		boo <- GetRandomSlice(N=foo$counts[i], R=c(foo$breaks[c(i,i+1)]), A=limits[1:2], B=limits[3:4])
		rand[j+1:foo$counts[i],] <- boo
		j <- j + foo$counts[i] 
	}
	hist(rand[,3], br=br, add=T, col='lightblue')
	#rand.Rsp <- bounders( rand[,1], rand[,2], rand[,3], c(R, QSO$S_limits) )
	#savePlot(filename=paste0(RAND_NAME, '_', br, '.png'))
	return(rand)
}	
