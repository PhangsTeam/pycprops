# pycprops
Python implementation of the CPROPS algorithm from [Rosolowsky & Leroy (2006)](https://ui.adsabs.harvard.edu/abs/2006PASP..118..590R/abstract).

This implements the same approaches developed in that paper, but there are a few changes in the shift to python:

- Local maximum identification and decimation uses the `astrodendro` dendrogram generation package, which is much master than the approach used in the original IDL version. 
- The assignment of pixels to clouds uses the `scikit-image` segmentation toolkit, which leads to different (better) pixel assignment around a local maxium.  The boundaries between clouds are more "natural" than they were in the original `cprops` code.
- Property extrapolation uses actual robust regression leading to <10% differences in some properties.
