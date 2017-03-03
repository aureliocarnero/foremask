import numpy as np
import pylab as pl
from astropy.wcs import WCS
from astropy.io import fits
import sys

filein = sys.argv[1]

sp = fits.open(filein)
sp.info()
exit()
header = sp[0].header

wcs = WCS(header)
print wcs
index = np.arange(header['NAXIS1'])

#print index
ra,dec = wcs.wcs_pix2world(index[:,np.newaxis],index[:,np.newaxis],0)

ramin = float(min(ra))
ramax = float(max(ra))
decmin = float(min(dec))
decmax = float(max(dec))

middlerac = (ramin+ramax)/2.
middledec = (decmin+decmax)/2.



fileout = open(filein + '.coords','w')

fileout.write('%.4f %.4f \n' %(ramin,decmin))
fileout.write('%.4f %.4f \n' %(ramax,decmin))
fileout.write('%.4f %.4f \n' %(ramin,decmax))
fileout.write('%.4f %.4f \n' %(ramax,decmax))

fileout.write('%.4f %.4f \n' %(middlerac,decmax))
fileout.write('%.4f %.4f \n' %(ramax,middledec))
fileout.write('%.4f %.4f \n' %(middlerac,decmin))
fileout.write('%.4f %.4f \n' %(ramin,middledec))
fileout.close()




'''
bitmask = sp[0].data

ra = np.array([item for sublist in ra for item in sublist])
dec = np.array([item for sublist in dec for item in sublist])
print min(ra),max(ra)
print min(dec),max(dec)

bitmask = np.array([item for sublist in bitmask for item in sublist])
print len(bitmask),len(ra),len(dec)

ratata,detata = [],[]
for i in range(len(sp[0].data)):
	de = [dec[i]]*len(sp[0].data)
	ratata.append(ra)
	detata.append(de)
ratata = np.array([item for sublist in ratata for item in sublist])
detata = np.array([item for sublist in detata for item in sublist])


print ratata[0],ratata[-1]
print detata[0],detata[-1]
mask_ = (bitmask!=0)
ratata = ratata[mask_]
detata = detata[mask_]
bitmask = bitmask[mask_]
pl.plot(ratata,detata,'.')
pl.savefig('3584p121.png')
#print ratata
#print detata	
print len(ratata),len(detata),len(bitmask)
'''
