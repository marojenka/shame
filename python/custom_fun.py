# from __future__ import print_function
# import numpy as np
# import sys
# 
# x = np.random.gamma(2,2,100000)
# y = np.random.gamma(1,1,100000)
# np.savetxt( '/tmp/tmp.txt', hist2d_gnuplot(x,y,(1000,1000)) )
# 
# xedges = np.arange( 1, 2000, 5 )
# yedges = np.arange( -22, -5, 0.01 )
# M = np.loadtxt('abr5Mmz.txt')[:,[2]+[5]]
# bins_x = xedges.size
# bins_y = yedges.size
# H, xedges, yedges = np.histogram2d(M[:,0], M[:,1], [xedges, yedges])
# #H, xedges, yedges = np.histogram2d(M[:,0], M[:,1], [1000, 1000])
# 
# # remap xedges and yedges to contain the bin center coordinates
# xedges = xedges[:-1] + 0.5*(xedges[1] - xedges[0])
# yedges = yedges[:-1] + 0.5*(yedges[1] - yedges[0])
# 
# tmp = np.vstack( (xedges, H.T) )
# result = np.insert( tmp, 0, np.array( [bins_y] + list(yedges) ), axis = 1 )
# np.savetxt( '/tmp/tmp.dat', result )

def hist2d_gnuplot(x, y, bins ) :
    #from numpy import histogram2d, vstack, insert
    H, xedges, yedges = np.histogram2d(x, y, [bins[0], bins[1]])
    # remap xedges and yedges to contain the bin center coordinates
    xedges = xedges[:-1] + 0.5*(xedges[1] - xedges[0])
    yedges = yedges[:-1] + 0.5*(yedges[1] - yedges[0])
    tmp = np.vstack( (xedges, H.T) )
    result = np.insert( tmp, 0, np.array( [yedges.size] + list(yedges) ), axis = 1 )
    return result

