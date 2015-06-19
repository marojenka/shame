import numpy as np
import pyfits as pf

# DAME_SDSS_DR9_testset.fits
# DAME_SDSS_DR9_trainset.fits
FILENAME = "DAME_SDSS_DR9_testset.fits"

f = pf.pyfits.open( FILENAME )

cols = f[1].columns
print cols.names
data = f[1].data
np.savetxt( FILENAME.replace("fits", "ascii"), data,  header = str(cols.names[:-1]).strip('[]') )
data.field(-1)[:] = 0
# ['objid', 'specObjID', 'ra', 'dec', 'psfMag_u', 'psfMag_g', 'psfMag_r', 'psfMag_i', 'psfMag_z', 'U-G', 'G-R', 'R-I', 'I-Z', 'zspec', 'zphot', 'type']
np.savetxt('/tmp/tmp.dat', np.array([data['ra'], data['dec']]).T)


