# Make VL samples from r and m 
foo <- read.table('/tmp/mr.dat')
index <- foo[,1] > 0 & foo[,2] > 0 
r <- foo[index, 1]
M <- foo[index, 2]


mr_info <- function(r, M, br=30) {
#{{{ mr_info plot mr and find VL sample with max N 
    histr <- hist(r,br=br, plot=F)
    binM <- histr$mids
    binCount <- histr$mids
    for( i in 1:length(histr$mids) ) {
        index <- r > histr$breaks[i] & r < histr$breaks[i+1] 
        if( sum(index) == 0 ) {
            binM[i] <- NA
            binCount[i] <- NA
        }
        else {
            binM[i] <- max(M[index])
            binCount[i] <- sum( M < binM[i] & r < histr$mids[i] )
        }
    }
    # plot(histr$mids, binM)
    # plot(histr$mids, binCount)
    # plot(histr$mids, binCount / histr$mids**3, t='l', log='y')
    
    i=which.max(binCount); index <- (M < binM[i] & r < histr$mids[i])
    plot(r, M, cex=.6, 
     ylim = c((tt <- range(M))[2], tt[1]), 
     main = paste( 'r_max = ', histr$breaks[i+1]), 
     xlab ='Radial distance', 
     ylab ='Absolute magnitude')
    lines(histr$mids, binM, col='blue', t='b')
    points(r[index], M[index], col='red')
#}}}
}

#plot(r[index], M[index], col='red')
#plot(histr)
#histr2 <- histr; histr2$counts[-i] = NA; plot(histr2, add=T, col='green')

#plot(r, M, cex=.6)

plot(r, M, cex=.6, 
    ylim = c((tt <- range(M))[2], tt[1]), 
    xlab ='Radial distance', 
    ylab ='Absolute magnitude')
#r.lim <- c( 50, histr$breaks[i+1], 200, 250 ) 
r.lim <- 100 
M.lim <- sapply( r.lim, function(x) { max( M[r > x] ) } )
VL.index <- lapply( 1:length(r.lim), 
            function(x) { 
                index <- r < r.lim[x] & M < M.lim[x]; 
                points(r[index], M[index], col=rainbow(length(r.lim))[x] ) 
                X <- r[index] * cos(a[index]*pi/180) * cos(b[index]*pi/180)
                Y <- r[index] * sin(a[index]*pi/180) * cos(b[index]*pi/180)
                Z <- r[index] * sin(b[index]*pi/180)
                write.table(cbind(X,Y,Z), file=paste0('VL', x-1, '.dat'), col.names=F, row.names=F)
                return( index )
            } ) 

plot(a, b, col='gray')
lapply( 1:length(r.lim), function(x) { points(a[ VL.index[[x]] ], b[ VL.index[[x]] ], col=rainbow(length(r.lim))[x], pch='.') }) 


