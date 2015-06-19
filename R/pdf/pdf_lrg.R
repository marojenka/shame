pgum <- function(y, alpha, betta ) {
    return( 1/betta * exp( - (y - alpha)/betta - exp( - (y - alpha)/betta ) ) )
}

pgumsc <- function( x ) {
    return( exp( -x - exp( -x ) ) ) 
}

gauss <- function(x, mu=0, sd=1) { 
    return( 1/( sd * sqrt(2*pi) ) * exp( - (x-mu)**2  / (2 * sd**2)  )  )
}

file.gamma <- '10_220-230.dat'
foo.gamma <- read.table( file.gamma ) 
ggplot(foo.gamma[ foo.gamma[,2] > 1, ], aes(V1, V2)) + geom_point() + scale_x_log10() + scale_y_log10()

files <- dir( patt = 'foo_lrg_1000-1050*' ) 
ang_scales <- read.table('scales.txt')
scales <- function ( scale = 0 ) {
    if( scale != 0 ) { 
        k = which.min( abs(ang_scales$V1 - scale) ) 
        cat( k ,' : ', ang_scales$V1[k], "\n" ) 
        return( k ) 
    } 
}


for( file in files ) { 
    print( file ) 
    foo <- read.table(file)
    
    dat <- foo[,3]
    hs <- hist( dat, br=20, plot=F ) 
    y <- hs$mids
    P <- hs$density
    mu <- mean( dat )
    del <- sd( dat )
    
    gamma <- - digamma(1)
    betta <- del * sqrt(6) / pi
    alpha <- mu - gamma * betta
    
    x <- (y - alpha) / betta 
    quu <- hist( ydat <- (dat - alpha) / betta , br=20, plot=F ) 
    
    chi <- chisq.test(quu$density, p = pgumsc( quu$mids ), res=T )
    print(chi$p.value)

    png(paste( sep='', './plots_pdf/', file, '.png'), width=1024, height=680 )
    plot( quu$mids, quu$density, log='y', main=paste(file, '\n p = ', chi$p.value ) )
    lines( x, pgumsc(x) , t='b', col='red')
    
    data <- data.frame( y = hs$mids, P = hs$density ) 
    fit <- nls( P ~ 1/betta * exp( - (y - alpha)/betta - exp( - (y - alpha)/betta ) ), 
               data = data, 
               start = list( alpha = alpha, betta = betta ) )
    alpha2 <- coefficients( fit )[1]
    betta2 <- coefficients( fit )[2]
    plot( hs$mids, hs$density, t='b' ) 
    lines( hs$mids, pgum( hs$mids, alpha2, betta2), t='b', col='red' ) 
    chi <- chisq.test(hs$density, p =  pgum( hs$mids, alpha2, betta2), res=T )
    fisher.test(hs$density, pgum( hs$mids, alpha2, betta2))
    x2 <- (y - alpha2) / betta2
    lines( x2, pgumsc(x2) , t='b', col='green')

    dev.off()
}

#lines( x, pgumsc(x)/betta , col='red')
#plot( hs$mids, P)
#lines( hs$mids, pgum( hs$mids, alpha, betta ), col='red', t='b' )

ab.ave <- function( data , br = 20 ) {
    hs <- hist( data, br, plot=F ) 
    y <- hs$mids
    P <- hs$density
    mu <- mean( data)
    del <- sd( data)
    
    gamma <- - digamma(1)
    betta <- del * sqrt(6) / pi
    alpha <- mu - gamma * betta
    
    x <- (y - alpha) / betta 
    quu <- hist( ydata <- (data - alpha) / betta , br=br, plot=F ) 
    
    return(list(alpha = alpha, 
               betta = betta, 
               x = quu$mids, 
               d = quu$density)) 
}

PDF <- function( i, br=20, GAUSS=FALSE ) { 
    foo1 <- read.table( files1[i] )
    
    s1 <- ab.ave( foo1[,3] , br = br)
    
    data <- data.frame()
    # s1.data 
    data <- rbind( data, data.frame(x=s1$x, P=s1$d, G=pgumsc(s1$x), sample="PDF" ) )

    pl <- ggplot(data=data, aes(x=x, y=P, col=sample)) + scale_y_log10() + geom_point() + geom_line()
    pl <- pl + geom_line( aes(x=x, y=G, colour="Gumbel"))
    pl <- pl + ggtitle(paste(files1[i],'\n', ang_scales$V1[i] * R0, 'Mpc'))
    if( GAUSS ) { 
        pl <- pl + geom_line( data=data.frame(x=s1$x, P=gauss(s1$x), sample="Gauss") ) 
    }
    print(pl)
}

PDFs <- function( index, br = 20 ) { 
    #print( cbind( index , ang_scales$V1[index] * R0) )
    scales <- cbind( index , ang_scales$V1[index] * R0) 
    sample_names <- cbind( index, paste("PDF", index, ang_scales$V1[index] * R0) )
    data <- data.frame()
    for( i in index ) { 
        s <- ab.ave( read.table( files[i] )[,3], br = br)
        data <- rbind( data, data.frame(x=s$x, P=s$d, sample = paste("PDF", sample_names[index == i,2]) ) )
        
        pl <- ggplot(data=data, aes(x=x, y=P, col=sample)) + scale_y_log10() + geom_point() + geom_line()
        #print( pl ) 
    }
    
    pl <- pl + scale_fill_discrete(breaks=sample_names[ sort( as.double(sample_names[,1]), index=T)$ix  ,2])

    x <- sort( data$x ) 
    gam <- pgumsc( x ) 
    gau <- gauss( x ) 
    gam[gam<1e-4] = NaN
    gau[gau<1e-4] = NaN
    ndata <- data.frame( x=x, gam=gam, gau=gau ) 
    pl <- pl + geom_line( data=ndata, aes(x=x, y=gam, col='Gumbel')) 
    pl <- pl + geom_line( data=ndata, aes(x=x, y=gau, col='Gauss' )) 
    print( pl )
    #return( pl ) 
}


R0 <- 220
files <- dir( patt = 'foo_S0_220-230_*' )
R0 <- 400
files <- dir( patt = 'foo_S0_400-410_*' )
files1 <- dir( patt = 'foo_s1_400-410_*' ) 
files2 <- dir( patt = 'foo_s2_400-410_*' ) 

for( i in length( files1 ) 
two.fields <- function( i, br=20 ) { 
    foo1 <- read.table( files1[i] )
    foo2 <- read.table( files2[i] )
    
    s1 <- ab.ave( foo1[,3] , br = br)
    s2 <- ab.ave( foo2[,3] , br = br)
    
    data <- data.frame()
    # s1.data 
    data <- rbind( data, data.frame(x=s1$x, P=s1$d, G=pgumsc(s1$x), sample="s1" ) )
    # s2.data <- data.frame(x=s2$x, P=s2$d, G=pgumsc(s2$x) )
    data <- rbind( data, data.frame(x=s2$x, P=s2$d, G=pgumsc(s2$x), sample="s2" ) )

    pl <- ggplot(data=data, aes(x=x, y=P, col=sample)) + scale_y_log10() + geom_point() + geom_line()
    pl <- pl + geom_line( aes(x=x, y=G, colour="Gumbel"))
    pl <- pl + ggtitle(paste(files1[i], files2[i]))
    print(pl)
}

