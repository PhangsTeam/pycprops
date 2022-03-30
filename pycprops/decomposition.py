import warnings
from astrodendro import Dendrogram, Structure, ppv_catalog
from spectral_cube import SpectralCube
from skimage.segmentation import watershed, random_walker, relabel_sequential
from skimage.morphology import disk
from astropy.io import fits
import astropy.units as u
import numpy as np
import astrodendro.pruning as pruning
import scipy.ndimage as nd
np.seterr(all='ignore')
#Shut up, I know what I'm doing
warnings.simplefilter("ignore")


def get_leaves(struct):
    if type(struct) is Dendrogram:
        return(struct.leaves)
    if type(struct) is Structure:
        return([desc for desc in struct.descendants if desc.is_leaf])

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
    struct = np.ones((2 * specfriends + 1, 1, 1), dtype=np.bool)
    maxfilt = nd.maximum_filter(maxfilt, footprint=struct)
    lmaxes = np.where((data == maxfilt) * (data != 0))
    return(lmaxes)


def deriv_decimate_leaves(d, cube, meta,
                          fscale=2, sigdiscont=0.5,
                          nredun=2, **kwargs):
    leaves = get_leaves(d)
    goodleaf = np.ones(len(leaves), dtype=np.bool)

    for i, leaf in enumerate(leaves):
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
    for i, leaf in enumerate(leaves):
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
                method='watershed',
                sigdiscont=0.5,
                verbose=True,
                **kwargs):

    if method=='nodecomp':
        lbimage = s.mask.include()
        wslabel, _ = nd.label(lbimage)
        wslabel, _, _  = relabel_sequential(wslabel)
        asgn = SpectralCube(wslabel, s.wcs, header=s.header)
        return(asgn)

    maxes = alllocmax(s, friends=int(friends), specfriends=specfriends)

    indep = pruning.all_true([pruning.min_npix(minpix),
                              pruning.min_delta(delta),
                              pruning.contains_seeds(maxes)])

    d = Dendrogram.compute(s.filled_data[:].value,
                           verbose=verbose,
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

    wslabel = np.zeros(s.shape, dtype=np.int)    
    runmax = 0
    for trunk in d.trunk:
        v, y, x = trunk.indices()
        if trunk.is_leaf:
            wslabel[v,y,x] = runmax + 1
            runmax += 1
        else:
            if sigdiscont > 0:
                leaves = deriv_decimate_leaves(trunk, s, meta, **kwargs)
            else:
                leaves = get_leaves(trunk)
                peaks = [leaf.get_peak()[0] for leaf in leaves]
                leaves = tuples_to_arrays(peaks)

            vals = s.filled_data[v, y, x].value
            maxval = np.nanmax(vals)
            lbimage = np.zeros(s.shape)
            lbimage[v,y,x] = vals 
            slcs = (nd.find_objects((lbimage != 0).astype(np.int)))[0]
            

            label = np.zeros(s.shape, dtype=np.int)
            for i in np.arange(len(leaves[0])):
                label[leaves[0][i], leaves[1][i], leaves[2][i]] = i + 1

            sublbimage = lbimage[slcs]
            labelmask = sublbimage != 0
            sublabel = label[slcs]
            if method == 'random_walker':
                sublabel[~labelmask] = -1
                sublabel = random_walker(sublbimage, sublabel)
                sublabel[~labelmask] = 0
            if method == 'watershed':
                sublabel = watershed(sublbimage, markers=sublabel, 
                                     compactness=compactness,
                                     mask=labelmask)
            sublabel, _,_ = relabel_sequential(sublabel)
            mxlabel = np.max(sublabel)
            sublabel[sublabel != 0] += runmax
            wslabel[slcs] += sublabel
            runmax += mxlabel

    # This shouldn't be needed, but I'm paranoid.
    wslabel, _, _  = relabel_sequential(wslabel)
    ncld = np.max(wslabel)

    # Be parsimonious about bit depth
    if ncld > 32767:
        wslabel = wslabel.astype(np.int32)
    else:
        wslabel = wslabel.astype(np.int16)

    asgn = fits.PrimaryHDU(wslabel, header=s.header)
    return(asgn)
