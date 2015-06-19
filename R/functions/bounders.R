bounders <- function (a, b, r, limit) {
	Rsp <- 1:length(r)
	for( i in 1:length(r) ) 
	{
		if( r[i] > limit[1] & r[i] < limit[2] )  
		{	
			alpha1 <- a[i] - limit[3]
			alpha2 <- limit[4] - a[i]
			if( alpha1 > 90 ) { r1 <- Inf } else { 
				r1 <- r[i] * cos( b[i]*pi/180 ) * sin( alpha1*pi/180 )
		    }	
			if( alpha2 > 90 ) { r2 <- Inf } else {
				r2 <- r[i] * cos( b[i]*pi/180 ) * sin( alpha2*pi/180 )
			}
			
			delta1 <- b[i] - limit[5]
			delta2 <- limit[6] - b[i]
			# paranoia
			if( delta1 < 0 || delta2 < 0) { cat("delta < 0", delta1, delta2, i, "\n"); return(42) }
			# next ifs a useless tbh. 
			if( delta1 > 90 ) { r3 <- Inf } else {
				r3 <- r[i] * sin( delta1 * pi/180 )
			}
			if( delta2 > 90 ) { r4 <- Inf } else {
				r4 <- r[i] * sin( delta2 * pi/180 )
			}
			
			r5 <-  r[i] - limit[1]
			r6 <- -r[i] + limit[2]
			
			Rsp[i] <- min(r1,r2,r3,r4,r5,r6)
		} else {
			Rsp[i] <- 0
		}
	}
	return(Rsp)
}

bounders.old <- function( a, b, r, limit) {
	max.R <- 1:length(r)
	for( i in 1:length(r) ) {
		if( (limit[1]<r[i]) * (r[i]<limit[2]) ) { 
			max.R[i] <- min(abs(cos(b[i]*pi/180)*sin((limit[3:4] - a[i])*pi/180)) * r[i], 
							abs(sin((limit[5:6] - b[i])* pi / 180)) * r[i], 
							abs(limit[1:2] - r[i]))
		}  else
			max.R[i] = 0
	}
	return(max.R)
}


bounders.test <- function (a, b, r, limit) {
	Rsp <- matrix(ncol=6, nrow=length(r))
	p   <- matrix(ncol=3, nrow=7)
	for( i in 1:length(r) ) 
	{
		if( r[i] > limit[1] & r[i] < limit[2] )  
		{	
			r1 <-  r[i] - limit[1]
			p[1,] <- CoordinatesSphericalToDecart(a[i], b[i], limit[1])[4:6]
			r2 <- -r[i] + limit[2]
			p[2,] <- CoordinatesSphericalToDecart(a[i], b[i], limit[2])[4:6]
			
			alpha1 <- a[i] - limit[3]
			alpha2 <- limit[4] - a[i]
			if( alpha1 > 90 ) { r3 <- Inf } else { 
				r3 <- r[i] * cos( b[i]*pi/180 ) * sin( alpha1*pi/180 )
				p[3,] <- CoordinatesSphericalToDecart(limit[3], b[i], r[i])[4:6]
		    }	
			if( alpha2 > 90 ) { r4 <- Inf } else {
				r4 <- r[i] * cos( b[i]*pi/180 ) * sin( alpha2*pi/180 )
				p[4,] <- CoordinatesSphericalToDecart(limit[4], b[i], r[i])[4:6]
			}
			
			delta1 <- b[i] - limit[5]
			delta2 <- limit[6] - b[i]
			# paranoia
			if( delta1 < 0 || delta2 < 0) { cat("delta < 0", delta1, delta2, "\n"); return(42) }
			# next ifs a useless tbh. 
			if( delta1 > 90 ) { r5 <- Inf } else {
				r5 <- r[i] * sin( delta1 * pi/180 )
				p[5,] <- CoordinatesSphericalToDecart(a[i], limit[5], r[i])[4:6]
			}
			if( delta2 > 90 ) { r6 <- Inf } else {
				r6 <- r[i] * sin( delta2 * pi/180 )
				p[6,] <- CoordinatesSphericalToDecart(a[i], limit[6], r[i])[4:6]
			}

			Rsp[i,1] <- min(r1,r2,r3,r4,r5,r6)
			j <- which.min(c( r1,r2,r3,r4,r5,r6))
			Rsp[i,2:4] <- p[j,]
			Rsp[i,5] <- j
		} else {
			Rsp[i,1] <- 0
			Rsp[i,2] <- c(00,00,00)
		}
		p[7,] <- CoordinatesSphericalToDecart(a[i], b[i], r[i])[4:6]
		l <- sqrt( sum((p[j,] - p[7,])^2) )
		Rsp[i,6] <- l
	}
	return(Rsp)
}


