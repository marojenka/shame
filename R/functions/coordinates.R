CoordinatesSphericalToDecart <- function( a, b, r ) {
	a <- a * pi / 180
	b <- b * pi / 180
	xyz <- cbind(
			cos(a) * cos(b) * r,
			sin(a) * cos(b) * r,
					 sin(b) * r 
	)
	a <- a / pi * 180
	b <- b / pi * 180
	return(cbind(a, b, r, xyz))
}


