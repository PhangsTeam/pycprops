# pycprops
Python implementation of the CPROPS algorithm from [Rosolowsky & Leroy (2006)](https://ui.adsabs.harvard.edu/abs/2006PASP..118..590R/abstract).

This implements the same approaches developed in that paper, but there are a few changes in the shift to python:

- Local maximum identification and decimation uses the `astrodendro` dendrogram generation package, which is much master than the approach used in the original IDL version. 
- The assignment of pixels to clouds uses the `scikit-image` segmentation toolkit, which leads to different (better) pixel assignment around a local maxium.  The boundaries between clouds are more "natural" than they were in the original `cprops` code.
- Property extrapolation uses actual robust regression leading to <10% differences in some properties.

# Installation

`pycprops` is an installable python package.  Clone or download the github repository.

```
git clone https://github.com/PhangsTeam/pycprops.git
cd pycprops
python setup.py install
```

After installation, you can run in python:

```
import pycprops
import astropy.units as u
cubefile = 'mycube.fits'  # Your cube
mask = 'mycube.mask.fits' # Mask defining where to find emission
d = 8.0 * u.kpc           # Distance (with units)

pycprops.fits2props(cubefile,
                    mask_file=mask,
                    distance=d,               
                    asgnname='mycube.asgn.fits',
                    propsname='mycube.props.fits')
```

The results of the `pycprops` code is an assignment file (`mycube.asgn.fits`) and a FITS BINTABLE (`mycube.props.fits`) catalog of GMC properties.  The tags in that catalog are described in the table below.

 | Tag Name | Description |
 |----------|------------|
 |NPIX | The total number of pixels in the cloud. |
 |CLOUDNUM | The cloud number from this data file. |
 | XCTR_PIX | The x position of the center of the cloud in pixels.|
 | XCTR_DEG | The x position of the center of the cloud in decimal degrees.|
 |YCTR_PIX            | The y position of the center of the cloud in pixels.                                         | 
 |YCTR_DEG            | The y position of the center of the cloud in decimal degrees.                                | 
 |VCTR_PIX            | The v position of the center of the cloud in pixels.                                         | 
 |VCTR_KMS            | The v position of the center of the cloud in km/s.                                           | 
 |POSANG              | The position angle of the cloud (measured from the X-axis)                                   | 
 |MOMXPIX             | The 2nd moment of the emission the x-direction in pixels.                                    | 
 |MOMXPIX_UC          | The fractional uncertainty of the 2nd moment in the x-direction.                             | 
 |MOMXPIX_NOEX        | The 2nd moment of the emission the x-direction in pixels with no extrapolation               | 
 |MOMXPIX_NOEX_UC     | The fractional uncertainty of the 2nd moment in the x-direction with no extrapolation.       | 
 |MOMYPIX             | The 2nd moment of the emission the y-direction in pixels.                                    | 
 |MOMYPIX_UC          | The fractional uncertainty of the 2nd moment in the y-direction.                             | 
 |MOMYPIX_NOEX        | The 2nd moment of the emission the y-direction in pixels with no extrapolation.              | 
 |MOMYPIX_NOEX_UC     | The fractional uncertainty of the 2nd moment in the y-direction with no extrapolation.       | 
 |MOMMAJPIX           | The 2nd moment of the emission along the major axis in pixels.                               | 
 |MOMMAJPIX_UC        | The fractional uncertainty of the 2nd moment in along the major axis.                        | 
 |MOMMAJPIX_NOEX      | The 2nd moment of the emission along the major axis in pixels with no extrapolation.         | 
 |MOMMAJPIX_NOEX_UC   | The fractional uncertainty of the 2nd moment along the major axis with no extrapolation.     | 
 |MOMMINPIX           | The 2nd moment of the emission along the minor axis in pixels.                               | 
 |MOMMINPIX_UC        | The fractional uncertainty of the 2nd moment in along the minor axis.                        | 
 |MOMMINPIX_NOEX      | The 2nd moment of the emission along the minor axis in pixels with no extrapolation.         | 
 |MOMMINPIX_NOEX_UC   | The fractional uncertainty of the 2nd moment along the minor axis with no extrapolation.     | 
 |MOMVPIX             | The 2nd moment of the emission the v-direction in pixels.                                    | 
 |MOMVPIX_UC          | The fractional uncertainty of the 2nd moment in the v-direction.                             | 
 |MOMVPIX_NOEX        | The 2nd moment of the emission the v-direction in pixels with no extrapolation.              | 
 |MOMVPIX_NOEX_UC     | The fractional uncertainty of the 2nd moment in the v-direction with no extrapolation.       | 
 |RAD_PC              | The cloud radius in parsecs.                                                                 | 
 |RAD_UC              | The fractional uncertainty in the cloud radius.                                              | 
 |RAD_NODC            | The cloud radius without deconvolution in parsecs.                                           | 
 |RAD_NODC_UC         | The fraction uncertainty in the cloud radius without deconvolution.                          | 
 |RAD_NOEX            | The radius without extrapolation in parsecs.                                                 | 
 |RAD_NOEX_UC         | The fractional uncertainty in the radius without extrapolation.                              | 
 |RAD_NODC_NOEX       | The radius without deconvolution or extrapolation in parsecs.                                | 
 |RAD_NODC_NOEX_UC    | The fraction uncertainty in the radius without deconvolution or extrapolation.               | 
 |SIGV_KMS            | The velocity dispersion in km/s.                                                             | 
 |SIGV_UC             | The fractional uncertainty in the velocity dispersion.                                       | 
 |SIGV_NODC           | The velocity dispersion without deconvolution in km/s.                                       | 
 |SIGV_NODC_UC        | The fraction uncertainty in the velocity dispersion without deconvolution.                   | 
 |SIGV_NOEX           | The velocity dispersion without extrapolation in km/s.                                       | 
 |SIGV_NOEX_UC        | The fraction uncertainty in the velocity dispersion without extrapolation.                   | 
 |SIGV_NODC_NOEX      | The velocity dispersion without deconvolution or extrapolation in km/s.                      | 
 |SIGV_NODC_NOEX_UC   | The fraction uncertainty in the velocity dispersion without extrapolation or deconvolution.  | 
 |FLUX_KKMS_PC2       | The luminosity of the cloud in K km/s pc2.                                                   | 
 |FLUX_UC             | The fractional uncertainty of the luminosity.                                                | 
 |FLUX_NOEX           | The luminosity of the cloud in K km/s pc2 without extrapolation.                             | 
 |FLUX_NOEX_UC        | The fractional uncertainty of the luminosity without extrapolation.                          | 
 |MLUM_MSUN           | The luminous (CO) mass of the cloud in M .                                                   | 
 |MLUM_UC             | The fractional uncertainty in the luminous mass.                                             | 
 |MVIR_MSUN           | The virial mass of the cloud in M .                                                          | 
 |MVIR_UC             | The fractional uncertainty in the virial mass.                                               | 
 |NAME                | The file name that was used to generate the properties.                                      | 
 |BEAMFWHM_PC         | The FWHM of the beam projected into parsecs.                                                 | 
 |SIGCHAN_KMS         | The equivalent velocity dispersion of a channel in km/s                                      | 
 |GALNAME             | The name of the galaxy that hosts the molecular cloud (can be defined later).                | 
 | RMSTORAD | The scaling factor used between 2nd moment in the position and the radius.|
 | PPBEAM |  The number of pixels per beam (area). | 
 | S2N | The peak signal-to-noise ratio in the cloud. |
 | TMAX_K | The peak antenna temperature in Kelvins. | 
 | DISTANCE_PC | The distance to the cloud in parsecs. |
