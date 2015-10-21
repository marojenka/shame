# test rmax
A <- c(0, 90) * d2r
B <- c(0, 90) * d2r
R <- c(0,  100)

tmp_abr <- random.sphere( N = 1e4 ) ; tmp_abr[,3] <- 1
tmp_xyz <- spherical_to_cartesian( tmp_abr )

check <- function( p = c(0,0,0), rlim = 1 ) {
    if( not_inside(cartesian_to_spherical(p)) ) return( NA )
    test <- tmp_xyz * rlim + matrix(p, ncol=3, nrow=length(tmp_xyz[,1]), by=T)
    return( sum( not_inside( cartesian_to_spherical(test)) ) )
}
N <- 1e3
centers_abr <- cbind( runif(N, A[1], A[2]), asin(runif(N,sin(B[1]),sin(B[2]))), ( 50 ) ) 
centers_xyz <- spherical_to_cartesian( centers_abr ) 
rlim <- rmax( centers_xyz ) 
out <- sapply( 1:length(rlim), function(x) {  check( centers_xyz[x,], rlim[x]  ) } )

wt( tmp_xyz )
wt( test, '/tmp/test.dat' )
wt( as.matrix(tmp_xyz)  + matrix(c(0,0,10), ncol=3, nrow=length(tmp_xyz[,1]), by=T), file = '/tmp/test.dat' )
wt( (as.matrix(tmp_xyz * rlim[1])  + matrix(centers_xyz[1,], ncol=3, nrow=length(tmp_xyz[,1]), by=T)), file = '/tmp/test.dat' )
p <- centers_xyz[1,]
tmp <- (as.matrix(tmp_xyz * rlim[1])  + matrix(centers_xyz[1,], ncol=3, nrow=length(tmp_xyz[,1]), by=T))
sum( not_inside( cartesian_to_spherical(tmp) ) )


