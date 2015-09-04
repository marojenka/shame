read_mad <- function(file = ('./VL1_north_mad.dat')) {
    foo <- read.table(file)
    boo <- as.matrix(foo)   
    res <- list(
                gr  = boo[1,-6:-1],
                r   = boo[-3:-1,3],
                rmax =  boo[-3:-1,2],
                rmaxind = boo[-3:-1,1],
                quu = as.matrix(boo[-3:-1,-6:-1])
    )
    return(res)
}
ScLe <- function(dd, scale) {
    rR_ind <- which.min(abs(dd$gr-scale))
    grid_lim <- dd$gr[rR_ind] 
    index <- dd$rmax > grid_lim 
    SL <- sapply(1:length(index)*index, function(ii) { sum(dd$quu[ii,1:rR_ind]) } )
    #SL <- sapply(1:length(index)*index, function(ii) { dd$quu[ii,rR_ind] } )
    #SL[SL == 0] = NA
    index <- SL != 0
    return( data.frame(r=dd$r[index], SL=SL[index]))
}
PDF <- function(dd, scale) {
    rR_ind <- which.min(abs(dd$gr-scale))
    grid_lim <- dd$gr[rR_ind] 
    index <- dd$rmax > grid_lim 
    SL <- sapply(1:length(index)*index, function(ii) { sum(dd$quu[ii,1:rR_ind]) } )
    #SL[SL == 0] = NA
    index <- SL != 0
    return( data.frame(r=dd$r[index], SL=SL[index]))

}
SL_mean <- function(SL, R, len=10) { 
    if( max(SL[,1]) - min(SL[,1]) < R ) {
        rbin <- c( mean(SL[,1])) 
    } else {
        rbin_tmp <- range(SL[,1]) + c(R, -R)/2
        rbin <- seq(rbin_tmp[1], rbin_tmp[2], len=len)
    }
    result <- data.frame()
    for( ii in 1:length( rbin ) ) {
        index <- (SL[,1] > rbin[ii] - R/2) & (SL[,1] < rbin[ii] + R/2); 
        m <- mean(SL[index,2])
        s <- sd(  SL[index,2])
        result <- rbind( result,   cbind(r=rbin[ii], mean=m, tsd=sqrt(m), sd=s  ) )
    }
    return( result )
}

# usage: 
f <- './VL2_south_mad.dat'
d <- read_mad(f)
SL <- ScLe(d, R <- 30)
m  <- SL_mean(SL, R)
plot(SLn, col='blue', pch='.', cex=3)
lines( mn[,1], mn[,2], col='red', t='b', pch=16 ) 
arrows(mn[,1], mn[,2] + mn[,3], mn[,1], mn[,2] - mn[,3], length=0.05, angle=90, code=3, col='red')
write.table(file=paste0('./VL2_north_SL_',R,'.dat'),      SLn, col.names=F, row.names=F)
write.table(file=paste0('./VL2_north_SL_',R,'_mean.dat'), mn, col.names=F, row.names=F)

###
#  PDF
library(ggplot2)
d <- read_mad(f)
xx <- data.frame()
pp <- data.frame()
for( R in seq(20,50,by=2) ) {
    SL <- ScLe(d, R)
    x <- (SL[,2] - mean(SL[,2])) / sd(SL[,2])
    xx <- rbind(xx, cbind(R=R, x=x, r=SL[,1], SL=SL[,2]) )
    foo <- hist(x, br=50)
    pp <- rbind(pp, cbind(x=foo$mids, P=foo$counts / sum(foo$counts), R=R) )
}
p <- ggplot(data=pp, aes(x=x, y=P, color=as.factor(R))) + geom_line()
print(p)
# PDF
########

################################
files <- c(
    '~/work/projects/2mrs/calc/gamma/VL2_north_mad.dat',
    '~/work/projects/2mrs/calc/gamma/VL2_south_mad.dat',
    '~/work/projects/millenium/data/g3_VL2/zspace_h.7_250_VL2_mad.dat',
    #'~/work/projects/tests/cantor/f2mrs/g3_VL2/Q.1_1_VL2_q1_mad.dat',
    '~/work/projects/tests/cantor/f2mrs/g3_VL2/Q.1_2_VL2_q1_mad.dat'
    #'~/work/projects/tests/cantor/f2mrs/g3_VL2/Q.1_3_VL2_q1_mad.dat'
)
rn <- array(dim=4)
for( f in files ) {
    d <- read_mad(f)
    SL <- ScLe(d, R <- 30)
    m  <- SL_mean(SL, R)
    plot(SL, col='blue', pch='.', cex=3, main=f, xlim=c(0,175), ylim=c(0,250))
    lines( m[,1], m[,2], col='red', t='b', pch=16 ) 
    arrows(m[,1], m[,2] + m[,3], m[,1], m[,2] - m[,3], length=0.05, angle=90, code=3, col='red')
    rn <- rbind(rn,c(range(SL$r), range(SL$SL)))
}

####
####
# Сделай рисунок канторовской статистики SL в виде нормированном на фактор f,
# равный f = N_millennium / N_cantor, где  N_millennium - среднее число точек в
# шарах выборки millennium (по всем бинам),   N_cantor -среднее по первым двум
# бинам SL кантора, тогда и средние значения, дисперсия и все точки SL для
# кантора и для millenium  можно сравнивать!  Сейчас там число точек в шарах
# этих выборок сильно отличаются.  Таким образом  ось y надо просто умножить на
# этот фактор y'  =  f y
foo <- lapply( files[3:4], function(x){read_mad(file=x)} )
R <- 30
SL <- lapply( foo, function(x){ScLe(dd=x, R)})
m  <- lapply( SL , function(x){SL_mean(SL=x, R)})

N_millennium <- mean(SL[[1]][,2])
N_cantor <- mean( SL[[2]] [ SL[[2]][,1] < 100,2] )
f <- N_millennium /  N_cantor
SLf <-cbind(SL[[2]][,1], SL[[2]][, 2] * f) 
mf <- SL_mean(SLf, R)
write.table(file=paste0('./VL2_cantor_f_',R,'.dat'),      SLf, col.names=F, row.names=F)
write.table(file=paste0('./VL2_cantor_f_',R,'_mean.dat'),  mf, col.names=F, row.names=F)

####
####
# Сравнивать SL статистики для различных распределений было бы правильнее
# в нормированном виде, когда за 1 принимается среднее число точек в тестовой сфере.
# Поэтому и надо разные графики SL  для случаев 2mrs-cantor-mock сначала
# отнормировать на среднее число точек в тестовом шаре данного распределения.
# Это тоже самое, что перенормировать ось у.
foo <- lapply( files[-2], function(x){read_mad(file=x)} )
R <- 30
SL <-  lapply( foo, function(x){ScLe(dd=x, R)})
f <-   lapply( SL , function(x){mean(x[,2])} )
SLf <- lapply( SL , function(x){return(cbind(x[,1], (x[,2] / mean(x[,2]))))} )
mf  <- lapply( SLf , function(x){SL_mean(SL=x, R)})

write.table(file=paste0('./VL2_2mrs_f_',R,'.dat'),      SLf[[1]], col.names=F, row.names=F)
write.table(file=paste0('./VL2_2mrs_f_',R,'_mean.dat'),  mf[[1]], col.names=F, row.names=F)
write.table(file=paste0('./VL2_crot_f_',R,'.dat'),      SLf[[2]], col.names=F, row.names=F)
write.table(file=paste0('./VL2_crot_f_',R,'_mean.dat'),  mf[[2]], col.names=F, row.names=F)
write.table(file=paste0('./VL2_cant_f_',R,'.dat'),      SLf[[3]], col.names=F, row.names=F)
write.table(file=paste0('./VL2_cant_f_',R,'_mean.dat'),  mf[[3]], col.names=F, row.names=F)





