from astrodendro import Dendrogram, ppv_catalog
from spectral_cube import SpectralCube
from skimage.segmentation import watershed
from skimage.morphology import disk
from astropy.io import fits
from astropy.stats import mad_std
import astropy.units as u
import numpy as np
import astrodendro.pruning as pruning
import scipy.ndimage as nd
from .cloudalyze import cloudalyze
from .decomposition import cube_decomp
import warnings
import os

#Shut up, I know what I'm doing
warnings.simplefilter("ignore")
np.seterr(all='ignore')

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
               allow_huge=False,
               asgn=None, **kwargs):

    if asgnname is None:
        asgnname = cube_file.replace('.fits', '_asgn.fits')
    if propsname is None:
        propsname = cube_file.replace('.fits', '_props.fits')

    s = SpectralCube.read(datadir + '/' + cube_file)
    mask = fits.getdata(datadir + '/' + mask_file)

    if allow_huge:
      s.allow_huge_operations = True
    
    if noise_file is None:
        print("No noise file found.  Calculating noise from Med. Abs. Dev of Data")
        noise = s.mad_std().value
        if delta is None:
            delta = 2 * noise.value
    else:
        noise = SpectralCube.read(datadir + '/' +noise_file)
        if delta is None:
            delta = 2 * noise.median().value

    distance_Mpc = distance.to(u.Mpc).value

    # Cast to boolean
    nanmask = np.isnan(mask)
    mask = mask.astype(np.bool)
    mask[nanmask] = False
    
    s = s.with_mask(mask)
    
    if asgn is None:
        asgn = cube_decomp(s, delta=delta, verbose=verbose,
                           **kwargs)
        asgn.writeto(output_directory + '/' + asgnname, overwrite=True)

    props = cloudalyze(s, asgn.data,
                       distance=distance,
                       verbose=verbose,
                       noise=noise,
                       alphaCO=alphaCO,
                       channelcorr=channelcorr,
                       **kwargs)
    
    props.write(output_directory + '/' + propsname, overwrite=True)
