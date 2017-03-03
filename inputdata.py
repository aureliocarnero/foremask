from __future__ import print_function

import numpy as np
import esutil
import astropy.io.fits as aft
#from . import standardValues


class InputMask(object):
    """
    A class to set the input mask

    parameters
    ----------
    maskin: Input mask file or list of files
        
    """
    def __init__(self, maskin):
        self.files = maskin

	self.file_dim = len(self.files)
	#Check file type, if one or various, what type: image? healpix? mangle?

        self.type = TYPE


	'''
        # get the band indices
        st=np.argsort(self.bands)
        temp=np.searchsorted(self.bands[st],self.zpts.zpts['BAND'][:])
        self._bind = st[temp]
	'''

    def checktpye(self, maskin):
	if isinstance(maskin, list):
	    print 'is a list of files'
	    self.checkfitstype(maskin, _list=True)
	else:
	    if maskin.endswidth('ply') or maskin.endswidth('pol'):
		print 'is a mangle file'
	    elif maskin.endswidth('fits'):
		print self.checkfitstype(maskin)

    def checkfitstype(self, maskin, _list=False):
	if not _list:
	    hdulist = aft.open(maskin)

	    hdulist.close()

	else:
	    for masc in maskin:
		hdulist = aft.open(masc)

		
