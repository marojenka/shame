# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
#from numpy import 
import numpypy
import numpy as np
# make spherical coordinates
def cartesian2spherical( foo ) :
#   void
#spherical_coordinats(int n, double *x, double *y, double *z, double *a, double *b, double *r) {
#	int i; 
#	for( i=0; i<n; i++ ) {
#		r[i] = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
#       a[i] = atan2(y[i], x[i]) ;  
#		if( a[i] < 0 ) 
#			a[i] += 2*M_PI; 
#		if( x[i] == 0 & y[i] == 0 ) {
#			if( z[i] > 0 ) b[i] = +M_PI/2.;
#			if( z[i] < 0 ) b[i] = -M_PI/2.;
#		}
#        b[i] = asin(z[i]/r[i]) ; 
#	}
#}
    r = np.sqrt(foo[:,0]**2+foo[:,1]**2+foo[:,2]**2)
    a = np.arctan2( foo[:,1], foo[:,0])
    b = np.arcsin(foo[:,2]/r)
    return() 








# read data file 
xyz = np.loadtxt('file.dat')
abr = cartesian2spherical(xyz)
