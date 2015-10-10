d2r <- pi / 180
r2d <- 180 / pi

rotate <- function( v, ax=0, ay=0, az=0 ) { 
# rotate a vector for a certan angle
    # v <- c(0,0,1)
    # rotate( v, ax=pi/2) c( 0, 1, 0)
    # rotate( v, ay=pi/2) c(-1, 0, 0)
    # rotate( c(0,0,1), az=pi/2) c( 0, -1, 0)
    Mx <- matrix( c(1,      0 ,        0, 
                    0, cos(ax), -sin(ax), 
                    0, sin(ax),  cos(ax)  ), 
                 ncol=3, nrow=3 )
    My <- matrix( c( cos(ay), 0 , sin(ay), 
                           0, 1 ,       0,
                    -sin(ay), 0 , cos(ay)  ), 
                 ncol=3, nrow=3 )
    Mz <- matrix( c( cos(az), -sin(az), 0, 
                     sin(az),  cos(az), 0, 
                           0,        0, 1 ), 
                 ncol=3, nrow=3 )
    result <- v;
    if( ax != 0 ) result  <- Mx %*% result;
    if( ay != 0 ) result  <- My %*% result;
    if( az != 0 ) result  <- Mz %*% result;

    return( result ) 
}

radius <- function(xyz) {
    xyz <- matrix(xyz, ncol=3)
    return( sqrt(xyz[,1]**2 + xyz[,2]**2 + xyz[,3]**2) )
}

spherical_to_cartesian <- function(abr) {
    abr <- matrix( abr, ncol=3 )
    
	xyz <- cbind( abr[,3] * cos(abr[,2]) * cos(abr[,1]), 
                  abr[,3] * cos(abr[,2]) * sin(abr[,1]),
                  abr[,3] * sin(abr[,2])                 )
    return( xyz )
}

cartesian_to_spherical <- function(xyz) {
    xyz <- matrix(xyz, ncol=3)
    abr <- matrix(nrow=length(xyz[,1]), ncol=3)
    abr[,3] <- radius(xyz)
    abr[,1] <- atan2(xyz[,2], xyz[,1]) 
    ind <- abr[,1]<0 
    abr[ind, 1] <- abr[ind,1] + 2 * pi
    ind <- ( xyz[,1] == 0 & xyz[,2] == 0 )
    abr[  ind, 2] <- sign(xyz[ ind,3]) * pi/2
    abr[ !ind, 2] <- asin(xyz[!ind,3] / abr[!ind, 3] )

    return( abr )
}

row.min <- function(x, y) {
    return( apply(cbind(x,y), 1, min) )
}

rmax <- function( xyz) {
    xyz <- matrix( xyz, ncol = 3 ) 
    N <- length(xyz[,1])
    rlim <- array( dim = N )
    abr <- cartesian_to_spherical( xyz ) 
	da <- 2*pi 
    db <- 2*pi

    outside <- not_inside(abr) 

    rlim <- row.min( abr[,3] - R[1] , R[2] - abr[,3] )
    if( sum(outside) < N )  {
        if( ( ( A[2] == A[1] + 2*pi ) & 
              ( B[2] == B[1] +   pi ) &
              ( R[1] == 0 ) ) ) { 
                rlim <- R[2] - abr[,3]
        } else {
            if( A[2] != A[1] + 2*pi ) 
                da <- row.min( abr[,1] - A[1], A[2] - abr[,1] )
            db <- row.min( abr[,2] - B[1], B[2] - abr[,2] )  
            ind <- da < pi/2
            rlim[ ind ] <- row.min( rlim[ind], 
                                    abr[ind,3] * cos(abr[ind,2])*sin(da[ind])  )
            ind <- db < pi/2
            rlim[ ind  ] <- row.min( rlim[ind], 
                                     abr[ind,3] * sin(db[ind]) )
        }
    }
    rlim[ outside ] <- -1
    return( rlim )
}

wt <- function( x, file='/tmp/tmp.dat' ) {
    write.table( x, file, col.names=F, row.names=F )
}

not_inside <- function(abr) {
    abr <- matrix(abr, ncol=3)
    index <-    abr[,1] > A[2] | abr[,1] < A[1] |
                abr[,2] > B[2] | abr[,2] < B[1] |
                abr[,3] > R[2] | abr[,3] < R[1]
    return( (index) ) 
}

random.sphere <- function( R = 1, N = 1e3  ) {
    return( cbind( runif(N, 0, 2*pi), 
                   asin(runif(N,-1,1)), 
                   (runif(N,0,R**3))**(1/3)) )
}

