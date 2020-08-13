import numpy as np
from spectral_cube import SpectralCube
from astropy.table import Table, Column
from astropy.io import fits
import astropy.units as u
import scipy.optimize as opt
from astropy.stats import mad_std
import matplotlib.pyplot as plt
import uuid
from astropy.utils.console import ProgressBar
from radio_beam import Beam
from collections import OrderedDict
import sys
# I too like to live dangerously
np.seterr(all='ignore')

sig2fwhm = np.sqrt(8 * np.log(2))

def cloudmom(x, y, v, t, target=0):
    moments = {}
    order = np.argsort(t)[::-1]
    t = t[order]
    x = x[order]
    y = y[order]
    v = v[order]

    meanx = np.cumsum(x*t)/np.cumsum(t)
    meany = np.cumsum(y*t)/np.cumsum(t)
    meanv = np.cumsum(v*t)/np.cumsum(t)
    mom0t = np.cumsum(t)

    term1x = np.cumsum(t * x**2)
    term2x = np.cumsum(t * x)**2/mom0t
    mom2x = np.sqrt((term1x - term2x)/mom0t)

    term1y = np.cumsum(t * y**2)
    term2y = np.cumsum(t * y)**2/mom0t
    mom2y = np.sqrt((term1y - term2y)/mom0t)

    term1v = np.cumsum(t * v**2)
    term2v = np.cumsum(t * v)**2/mom0t
    mom2v = np.sqrt((term1v - term2v)/mom0t)

    xcen = np.sum(x * t) / np.sum(t)
    ycen = np.sum(y * t) / np.sum(t)
    vcen = np.sum(v * t) / np.sum(t)


    moments['xcen'] = xcen
    moments['ycen'] = ycen
    moments['vcen'] = vcen
    
    moments['rmsx_noex'] = mom2x[-1]
    moments['rmsy_noex'] = mom2y[-1]
    moments['rmsv_noex'] = mom2v[-1]
    moments['flux_noex'] = mom0t[-1]
    slc = slice(int(len(t) / 2), None, None)
    moments['rmsx_ex'] = extrap(t, mom2x, order=1, 
                                targett=target, slc=slc)
    moments['rmsy_ex'] = extrap(t, mom2y, order=1, 
                                targett=target, slc=slc)
    moments['rmsv_ex'] = extrap(t, mom2v, order=1, 
                                targett=target, slc=slc)
    moments['flux_ex'] = extrap(t, mom0t, order=2, 
                                targett=target, slc=slc)

    # fig, axlist = plt.subplots(2, 2)
    # fig.set_size_inches(6, 6)
    # axlist = axlist.flatten()
    # axlist[0].plot(t, mom2x, 'ro')
    # axlist[0].plot(0, moments['rmsx_ex'], 'bx')
    # axlist[1].plot(t, mom2y, 'ro')
    # axlist[1].plot(0, moments['rmsy_ex'], 'bx')
    # axlist[2].plot(t, mom2v, 'ro')
    # axlist[2].plot(0, moments['rmsv_ex'], 'bx')
    # axlist[3].plot(t, mom0t, 'ro')
    # axlist[3].plot(0, moments['flux_ex'], 'bx')
    
    # unique_filename = str(uuid.uuid4())
    # plt.tight_layout()
    # plt.savefig('figs/' + unique_filename + '.png')
    # plt.close()
    # plt.clf()
    return(moments)

def polyloss(p, x, y, order=1):
    order = np.argsort(x)[::-1] + 1
    wt = np.sqrt(1 / order)
    return((y - polynomial(p, x, order=order)) * wt)


def polynomial(p, x, order=1):
    y = 0
    for i in range(len(p)):
        y = y * x + p[i] 
    return(y)


def extrap(x, y, targett=0, order=1, slc=slice(5, None, None)):
    good = np.isfinite(x) * np.isfinite(y)
    if len(x[good][slc]) == 0:
        return(np.nan)
    coeffs = np.polyfit(x[good][slc], y[good][slc], order)
    result = opt.least_squares(polyloss, coeffs, loss='arctan',
                               args=(x[good][slc], y[good][slc]),
                               kwargs={'order':order})
    if result.success:
        extrapval = polynomial(result.x, targett, order=order)
        if extrapval < y[good][slc].max():
            return(y[good][slc].max())
        return(polynomial(result.x, targett, order=order))
    else:
        return(np.nan)
    

def pa_moment(x, y, t):
    major = 0
    minor = 0
    x0 = np.nansum(t * x)/np.nansum(t)
    y0 = np.nansum(t * y)/np.nansum(t)
    wt = t
    mat = 1./(np.nansum(wt)) * np.array(
        [[np.nansum(wt * (x - x0)**2), np.nansum(wt*(x-x0) * (y - y0))],
         [np.nansum(wt * (x - x0) * (y - y0)), np.nansum(wt * (y - y0)**2)]])
    det = np.linalg.det(mat)
    if (det == 0) or np.isnan(det):
        return([np.nan, np.nan, np.nan])
    evals, evec = np.linalg.eig(mat)
    # if (np.nansum(transpose(mat) eq mat) eq 4 and mat[1, 0] eq 0) then $
    #     evec = [[1., 0], [0., 1.]] else $
    #     evec = eigenvec(mat, evals, residual = res, /double)
    big = np.argmax(evals)
    bigvec = evec[:, big]
    posang = np.arctan2(bigvec[1], bigvec[0])
    if posang < 0:
        posang += np.pi
    major = np.sqrt(evals[big])
    minor=np.sqrt(evals[1 - big])
    return(posang, major, minor)
    
def deconvolve_moments(moments, beam, pa, pixscale):
    if moments['rmsx_ex'] >= moments['rmsy_ex']:
        beam_ex = Beam(major=moments['rmsx_ex'] * pixscale * sig2fwhm,
                       minor=moments['rmsy_ex'] * pixscale * sig2fwhm,
                       pa=pa*u.rad)
    else:
        beam_ex = Beam(major=moments['rmsy_ex'] * pixscale * sig2fwhm,
                       minor=moments['rmsx_ex'] * pixscale * sig2fwhm,
                       pa=(pa + np.pi/2)*u.rad)
    beam_ex_dc = beam_ex.deconvolve(beam, failure_returns_pointlike=True)
    
    if moments['rmsx_noex'] >= moments['rmsy_noex']:
        beam_noex = Beam(major=moments['rmsx_noex'] * pixscale * sig2fwhm,
                         minor=moments['rmsy_noex'] * pixscale * sig2fwhm,
                         pa=pa*u.rad)
    else:
        beam_noex = Beam(major=moments['rmsy_noex'] * pixscale * sig2fwhm,
                          minor=moments['rmsx_noex'] * pixscale * sig2fwhm,
                          pa=(pa + np.pi/2)*u.rad)
    beam_noex_dc = beam_noex.deconvolve(beam, failure_returns_pointlike=True)
    
    outdict = {}
    outdict['rmsx_ex'] = (
        beam_ex_dc.major 
        / pixscale
        / sig2fwhm).to(u.dimensionless_unscaled).value
    outdict['rmsy_ex'] = (
        beam_ex_dc.minor 
        / pixscale
        / sig2fwhm).to(u.dimensionless_unscaled).value
    outdict['pa_ex'] = (
        beam_ex_dc.pa).to(u.rad).value
    outdict['rmsx_noex'] = (
        beam_noex.major 
        / pixscale
        / sig2fwhm).to(u.dimensionless_unscaled).value
    outdict['rmsy_noex'] = (
        beam_noex.minor 
        / pixscale
        / sig2fwhm).to(u.dimensionless_unscaled).value
    outdict['pa_noex'] = (
        beam_noex.pa).to(u.rad).value
    return(outdict)

def cloudalyze(cube, label, 
               distance=None,
               alphaCO=6.7, extrapolate=True,
               bootstrap=0, rmstorad=1.91,
               noise=None, verbose=True,
               channelcorr=0.0,
               **kwargs):
    cloudlist = []
    spectral_unit = cube.spectral_axis
    dv = np.abs((cube.spectral_axis[1] 
                 - cube.spectral_axis[0]).to(u.km / u.s).value)
    apix_sr =(np.abs(np.linalg.det(cube.wcs.celestial.pixel_scale_matrix)) 
              * u.deg**2).to(u.sr)
    pixscale = np.sqrt(np.abs(np.linalg.det(
        cube.wcs.celestial.pixel_scale_matrix))) * u.deg
    apix = (distance**2 * apix_sr / u.sr).to(u.pc**2).value
    dx = np.sqrt(apix)
    bmaj_pc = (cube.beam.major.to(u.rad) / u.rad * distance).to(u.pc).value
    bmin_pc = (cube.beam.minor.to(u.rad) / u.rad * distance).to(u.pc).value
    beamfwhm_pc = np.sqrt(bmaj_pc * bmin_pc)
    k = (0.47 * channelcorr
         -0.23 * channelcorr**2
         -0.16 * channelcorr**3
         +0.43 * channelcorr**4)
    
    sigchan = dv / np.sqrt(2 * np.pi) * (1 + 1.18 * k + 10.4 * k**2)
    uniqlabels = np.unique(label)
    # if verbose:
    #     bar = ProgressBar(len(uniqlabels))

    cloudnum = 0

    if verbose:
        barfile = sys.stdout
        print("Calculating cloud properties for {0} clouds".format(len(uniqlabels)))
    else:
        barfile = None
    
    for thislabel in ProgressBar(uniqlabels, file=barfile):
        
        if thislabel == 0:
            continue
        thiscloud = OrderedDict()
        (v, y, x,) = np.where(label == thislabel)

        t = cube.filled_data[v, y, x].value
        if noise is not None:
            if isinstance(noise, float):
                s2n = np.nanmax(t / noise)
                thiscloud['S2N'] = s2n
            elif (noise.shape == cube.shape):
                rms = noise.filled_data[v, y, x].value
                s2n = np.nanmax(t / rms)
                thiscloud['S2N'] = s2n
            else:
                thiscloud['S2N'] = np.nan
        cloudnum += 1
        thiscloud['NPIX'] = len(x)
        thiscloud['CLOUDNUM'] = cloudnum
        thiscloud['DISTANCE_PC'] = distance.to(u.pc).value
        thiscloud['BEAMFWHM_PC'] = beamfwhm_pc
        thiscloud['BEAMMAJ_PC'] = bmaj_pc
        thiscloud['BEAMMIN_PC'] = bmin_pc
        thiscloud['RMSTORAD'] = rmstorad
        thiscloud['PPBEAM'] = cube.pixels_per_beam
        thiscloud['SIGCHAN_KMS'] = sigchan
        thiscloud['TMAX_K'] = np.nanmax(t)
    
        moments = cloudmom(x, y, v, t)
        thiscloud['XCTR_PIX'] = moments['xcen']
        thiscloud['YCTR_PIX'] = moments['ycen']
        thiscloud['VCTR_PIX'] = moments['vcen']
        a, d, vnative = cube.wcs.wcs_pix2world(moments['xcen'],
                                               moments['ycen'],
                                               moments['vcen'], 0)
        vkms = u.Quantity(vnative, cube.spectral_axis.unit).to(u.km / u.s).value
        thiscloud['XCTR_DEG'] = a
        thiscloud['YCTR_DEG'] = d
        thiscloud['VCTR_KMS'] = vkms
        
        thiscloud['MOMYPIX_NOEX'] = moments['rmsy_noex']
        thiscloud['MOMYPIX'] = moments['rmsy_ex']
        thiscloud['MOMVPIX_NOEX'] = moments['rmsv_noex']
        thiscloud['MOMVPIX'] = moments['rmsv_ex']
        thiscloud['MOMXPIX_NOEX'] = moments['rmsx_noex']
        thiscloud['MOMXPIX'] = moments['rmsx_ex']

        thiscloud['FLUX_NOEX'] = moments['flux_noex'] * dv * apix
        thiscloud['FLUX_KKMS_PC2'] = moments['flux_ex'] * dv * apix
        
        thiscloud['SIGV_KMS'] = moments['rmsv_ex'] * dv
        
        pa, _, _ = pa_moment(x, y, t)

        xrot = x * np.cos(pa) + y * np.sin(pa)
        yrot =-x * np.sin(pa) + y * np.cos(pa)
        moments_rot = cloudmom(xrot, yrot, v, t)
        moments_rot_dc = deconvolve_moments(moments_rot,
                                            cube.beam,
                                            pa,
                                            pixscale)
        resolved = (moments_rot_dc['rmsx_ex'] > 0)
        noex_resolved = (moments_rot_dc['rmsx_noex'] > 0)

        thiscloud['MOMMAJPIX'] = moments_rot_dc['rmsx_ex']
        thiscloud['MOMMAJPIX_NOEX'] = moments_rot_dc['rmsx_noex']
        thiscloud['MOMMAJ_PC'] = moments_rot_dc['rmsx_ex'] * dx

        thiscloud['MOMMAJ_NODC_PC'] = moments_rot['rmsx_ex'] * dx
        thiscloud['MOMMAJPIX_NODC'] = moments_rot['rmsx_ex']
        thiscloud['MOMMAJPIX_NODC_NOEX'] = moments_rot['rmsx_noex']


        thiscloud['MOMMINPIX'] = moments_rot_dc['rmsy_ex']
        thiscloud['MOMMINPIX_NOEX']=moments_rot_dc['rmsy_noex']
        thiscloud['MOMMIN_PC'] = moments_rot_dc['rmsy_ex'] * dx

        thiscloud['MOMMAJ_NODC_PC'] = moments_rot['rmsx_ex'] * dx
        thiscloud['MOMMINPIX_NODC'] = moments_rot['rmsy_ex']
        thiscloud['MOMMINPIX_NODC_NOEX'] = moments_rot['rmsy_noex']


        thiscloud['POSANG'] = pa
        if not resolved:
            thiscloud['PA_DC'] = np.nan
            thiscloud['FWHM_MIN_DC'] = np.nan
            thiscloud['FWHM_MAJ_DC'] = np.nan
            thiscloud['RAD_PC'] = np.nan
        else:        
            thiscloud['PA_DC'] = moments_rot_dc['pa_ex']
            thiscloud['FWHM_MIN_DC'] = moments_rot_dc['rmsy_ex'] * dx * sig2fwhm
            thiscloud['FWHM_MAJ_DC'] = moments_rot_dc['rmsx_ex'] * dx * sig2fwhm
            thiscloud['RAD_PC'] = dx * rmstorad * np.sqrt(
                moments_rot_dc['rmsy_ex'] * moments_rot_dc['rmsx_ex']
            )
        if not noex_resolved:
            thiscloud['RAD_NOEX'] = np.nan
        else:
            thiscloud['RAD_NOEX'] = dx * rmstorad * np.sqrt(
                moments_rot_dc['rmsy_noex'] * moments_rot_dc['rmsx_noex']
            )

        thiscloud['RAD_NODC'] = (np.sqrt(moments_rot['rmsx_ex'] 
                                        * moments_rot['rmsy_ex'])
                                 * rmstorad * dx)
        thiscloud['RAD_NODC_NOEX'] = (np.sqrt(moments_rot['rmsx_noex']
                                             * moments_rot['rmsy_noex']) 
                                      * rmstorad * dx)

        thiscloud['SIGV_KMS'] = np.sqrt(moments_rot['rmsv_ex']**2 
                                        - sigchan**2) * dv
        thiscloud['SIGV_NOEX'] = np.sqrt(moments_rot['rmsv_noex']**2
                                        - sigchan**2) * dv
        thiscloud['SIGV_NODC_NOEX'] = moments_rot['rmsv_noex'] * dv
        thiscloud['SIGV_NODC'] = moments_rot['rmsv_ex'] * dv

        if callable(alphaCO):
            thisalphaCO = alphaCO(thiscloud['XCTR_DEG'],
                                  thiscloud['YCTR_DEG'])
            thisalphaCO = thisalphaCO.to(u.M_sun
                                         / (u.K * u.km
                                            / u.s * u.pc**2)).value
            thiscloud['ALPHA_CO'] = thisalphaCO
        else:
            thisalphaCO = alphaCO
            thiscloud['ALPHA_CO'] = thisalphaCO
            
        thiscloud['MLUM_MSUN'] = thisalphaCO * thiscloud['FLUX_KKMS_PC2']
        thiscloud['MVIR_MSUN'] = (1040 
                                  * thiscloud['RAD_PC']
                                  * thiscloud['SIGV_KMS']**2)
        
        if bootstrap > 0:
            bootlist = []
            bootlist_rot = []
            N = len(x)
            idx = np.arange(N)
            for i in range(bootstrap):
                subset = np.random.choice(idx, size=N, replace=True)
                bootlist += [cloudmom(x[subset],
                                      y[subset],
                                      v[subset],
                                      t[subset])]
                pa, _, _ = pa_moment(x[subset], y[subset], t[subset])
                xrot = x * np.cos(pa) + y * np.sin(pa)
                yrot = -x * np.sin(pa) + y * np.cos(pa)
                bootlist_rot += [cloudmom(xrot[subset],
                                          yrot[subset],
                                          v[subset],
                                          t[subset])]
            boottable = Table(bootlist)
            boottable_rot = Table(bootlist_rot)
            uc_dict = {}
            uc_dict_rot = {}
            for k in boottable.keys():
                ucval = (mad_std(boottable[k]) 
                         / np.median(boottable[k]) 
                         * cube.pixels_per_beam**0.5)
                ucval_rot = (mad_std(boottable_rot[k]) 
                            / np.median(boottable_rot[k]) 
                            * cube.pixels_per_beam**0.5)
                uc_dict[k] = ucval
                uc_dict_rot[k] = ucval_rot
            thiscloud['MLUM_UC'] = uc_dict['flux_ex']
            thiscloud['FLUX_UC'] = uc_dict['flux_ex']
            thiscloud['FLUX_NOEX_UC'] = uc_dict['flux_noex']
            thiscloud['SIGV_NOEX_UC'] = uc_dict['rmsv_noex']
            thiscloud['SIGV_NODC_NOEX_UC'] = uc_dict['rmsv_noex']
            thiscloud['SIGV_UC'] = uc_dict['rmsv_ex']
            thiscloud['SIGV_NODC_UC'] = uc_dict['rmsv_ex']
            thiscloud['MOMXPIX_UC'] = uc_dict['rmsx_ex']
            thiscloud['MOMXPIX_NOEX_UC'] = uc_dict['rmsx_noex']
            thiscloud['MOMYPIX_UC'] = uc_dict['rmsy_ex']
            thiscloud['MOMYPIX_NOEX_UC'] = uc_dict['rmsy_noex']
            thiscloud['MOMVPIX_UC'] = uc_dict['rmsv_ex']
            thiscloud['MOMVPIX_NOEX_UC'] = uc_dict['rmsv_noex']
            thiscloud['MOMMAJPIX_UC'] = uc_dict_rot['rmsx_ex']
            thiscloud['MOMMAJPIX_NOEX_UC'] = uc_dict_rot['rmsx_noex']
            thiscloud['MOMMINPIX_UC']=uc_dict_rot['rmsy_ex']
            thiscloud['MOMMINPIX_NOEX_UC']=uc_dict_rot['rmsy_noex']
            thiscloud['RAD_UC']=np.sqrt(uc_dict_rot['rmsx_ex']**2 +
                                        uc_dict_rot['rmsx_ex']**2) * 0.5
            thiscloud['RAD_NODC_UC'] = np.sqrt(uc_dict_rot['rmsx_ex']**2 +
                                               uc_dict_rot['rmsx_ex']**2) * 0.5
            thiscloud['RAD_NOEX_UC']=np.sqrt(uc_dict_rot['rmsx_noex']**2 +
                                             uc_dict_rot['rmsx_noex']**2) * 0.5
            thiscloud['RAD_NODC_NOEX_UC'] = np.sqrt(uc_dict_rot['rmsx_noex']**2 +
                                                    uc_dict_rot['rmsx_noex']**2) * 0.5
            thiscloud['MVIR_UC']=np.sqrt(thiscloud['RAD_UC']**2 
                                         + 2*thiscloud['SIGV_UC']**2)

        cloudlist += [thiscloud]
    outtable = Table(cloudlist)
    return(outtable)

