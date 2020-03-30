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
               asgnsuffix='_asgn',
               propsuffix='_props',
               delta=None,
               asgn=None, **kwargs):

    s = SpectralCube.read(datadir + '/' + cube_file)
    mask = fits.getdata(datadir + '/' + mask_file)
    noise = SpectralCube.read(datadir + '/' +noise_file)

    distance_Mpc = distance.to(u.Mpc).value

    s = s.with_mask(mask.astype(np.bool))    
    
    if delta is None:
        delta = 2 * noise.median().value
    
    if asgn is None:
        asgn = cube_decomp(s, delta=delta, **kwargs)
        output_asgn_name = cube_file.replace('.fits', asgnsuffix + '.fits')
        asgn.write(output_directory + '/' + output_asgn_name, overwrite=True)

    props = cloudalyze(s, asgn.filled_data[:].value, distance=distance, **kwargs)
    output_props_name = cube_file.replace(
        '.fits', propsuffix + '.fits')
    props.write(output_directory + '/' + output_props_name, overwrite=True)
