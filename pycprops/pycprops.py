from astrodendro import Dendrogram, ppv_catalog
from spectral_cube import SpectralCube
from skimage.segmentation import watershed
from skimage.morphology import disk
from astropy.io import fits
import astropy.units as u
import numpy as np
import astrodendro.pruning as pruning
import scipy.ndimage as nd
from .cloudalyze import cloudalyze
from .decomposition import cube_decomp
np.seterr(all='ignore')
import warnings
import os
#Shut up, I know what I'm doing
warnings.simplefilter("ignore")


def fits2props(cube_file,
               datadir=os.getcwd(),
               output_directory=os.getcwd(),
               mask_file=None,
               noise_file=None,
               distance=None,
               asgnname=None,
               propsname=None,
               delta=None,
               verbose=True,
               alphaCO=6.7,
               channelcorr=0.0,
               asgn=None, **kwargs):

    if asgnname is None:
        asgnname = cube_file.replace('.fits', '_asgn.fits')
    if propsname is None:
        propsname = cube_file.replace('.fits', '_props.fits')

    s = SpectralCube.read(datadir + '/' + cube_file)
    mask = fits.getdata(datadir + '/' + mask_file)
    noise = SpectralCube.read(datadir + '/' +noise_file)

    distance_Mpc = distance.to(u.Mpc).value

    s = s.with_mask(mask.astype(np.bool))    
    
    if delta is None:
        delta = 2 * noise.median().value
    
    if asgn is None:
        asgn = cube_decomp(s, delta=delta, verbose=verbose, **kwargs)
        asgn.write(output_directory + '/' + asgnname, overwrite=True)

    props = cloudalyze(s, asgn.filled_data[:].value,
                       distance=distance,
                       verbose=verbose,
                       noise=noise,
                       alphaCO=alphaCO,
                       channelcorr=channelcorr,
                       **kwargs)
    
    props.write(output_directory + '/' + propsname, overwrite=True)
