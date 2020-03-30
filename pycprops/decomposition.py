import warnings
from astrodendro import Dendrogram, ppv_catalog
from spectral_cube import SpectralCube
from skimage.segmentation import watershed
from skimage.morphology import disk
from astropy.io import fits
import astropy.units as u
import numpy as np
import astrodendro.pruning as pruning
import scipy.ndimage as nd
np.seterr(all='ignore')
#Shut up, I know what I'm doing
warnings.simplefilter("ignore")



def tuples_to_arrays(peaks):
    x0 = []
    x1 = []
    x2 = []
    for peak in peaks:
        x0 += [peak[0]]
        x1 += [peak[1]]
        x2 += [peak[2]]
    x0 = np.asarray(x0)
    x1 = np.asarray(x1)
    x2 = np.asarray(x2)
    return((x0, x1, x2))


def alllocmax(cube, friends=1, specfriends=1):
    data = np.nan_to_num(cube.filled_data[:].value)
    struct = disk(friends)
    maxfilt = nd.maximum_filter(data, footprint=struct[np.newaxis, :])
    struct = np.ones((2 * friends + 1, 1, 1), dtype=np.bool)
    maxfilt = nd.maximum_filter(maxfilt, footprint=struct)
    lmaxes = np.where((data == maxfilt) * (data != 0))
    return(lmaxes)


def deriv_decimate_leaves(d, cube, meta,
                          fscale=2, sigdiscont=0.5, nredun=2):
    goodleaf = np.ones(len(d.leaves), dtype=np.bool)

    for i, leaf in enumerate(d.leaves):
        if not goodleaf[i] or (leaf.parent is None):
            continue
        parentobj = ppv_catalog([leaf.parent], meta, verbose=False)
        children = ppv_catalog(leaf.parent.children, meta, verbose=False)
        thisleaf = np.array(children['_idx']) == leaf.idx
        dmajor = np.array((parentobj['major_sigma'] - children['major_sigma'])
                          / children['major_sigma']) > sigdiscont
        dminor = np.array((parentobj['minor_sigma'] - children['minor_sigma'])
                          / children['minor_sigma']) > sigdiscont
        dsigv = np.array((parentobj['v_rms'] - children['v_rms'])
                         / children['v_rms']) > sigdiscont
        dflux = np.array((parentobj['v_rms'] - children['v_rms'])
                         / children['v_rms']) > (sigdiscont * fscale)
        disc = np.any(np.sum(np.c_[dmajor, dminor, dsigv, dflux], axis=1)
                      >= nredun)
        goodleaf[i] = disc
    leaflist = []
    for i, leaf in enumerate(d.leaves):
        if goodleaf[i]:
            leaflist += [leaf]
    peaks = [leaf.get_peak()[0] for leaf in leaflist]
    indices = tuples_to_arrays(peaks)

    return(indices)


def cube_decomp(s, 
                minpix=5, 
                delta=2, 
                specfriends=1, 
                friends=3, 
                compactness=1, 
                **kwargs):
    maxes = alllocmax(s, friends=int(friends), specfriends=1)

    indep = pruning.all_true([pruning.min_npix(minpix),
                              pruning.min_delta(delta),
                              pruning.contains_seeds(maxes)])

    d = Dendrogram.compute(s.filled_data[:].value,
                           verbose=True,
                           is_independent=indep)

    meta = {}
    meta['spatial_scale'] = s.header['CDELT2'] * u.deg
    meta['beam_major'] = s.beam.major
    meta['beam_minor'] = s.beam.minor
    meta['data_unit'] = s.unit
    meta['vaxis'] = 0
    meta['wcs'] = s.wcs

    lam = np.median(s.with_spectral_unit(u.mm,
                                         velocity_convention='radio').spectral_axis)
    meta['wavelength'] = lam

    leaves = deriv_decimate_leaves(d, s, meta)
    label = np.zeros(s.shape, dtype=np.int)

    for i in np.arange(len(leaves[0])):
        label[leaves[0][i], leaves[1][i], leaves[2][i]] = i + 1
    lbimage = -s.filled_data[:].value
    maxval = np.nanmax(lbimage)
    baddata = np.isnan(lbimage)
    lbimage[baddata] = maxval+1
    wslabel = watershed(lbimage, markers=label, compactness=compactness)
    wslabel[baddata] = 0
    asgn = SpectralCube(wslabel, s.wcs, header=s.header)
    return(asgn)