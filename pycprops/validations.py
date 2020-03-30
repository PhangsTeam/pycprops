
# cube = SpectralCube.read('ngc3137_12m+7m+tp_co21_pbcorr_round_k.fits')
# rms = SpectralCube.read('ngc3137_12m+7m+tp_co21_noise_pbcorr_round_k.fits')
# label = fits.getdata(
#     'ngc3137_12m+7m+tp_co21_pbcorr_round_k_cprops_v3p4_nosigdiscont_asgn.fits')

# t = cloudalyze(cube, label, distance=10.569999694824219 * u.Mpc,
#                noise=rms, bootstrap=30)
# t2 = Table.read(
#     'ngc3137_12m+7m+tp_co21_pbcorr_round_k_cprops_v3p4_nosigdiscont_props.fits')
# t.write('test.fits', overwrite=True)
# notimplist = []
# for k in t2.keys():
#     try:
#         vec = t[k]
#         print(k)
#         print('    Median:  ', np.nanmedian(np.abs(vec - t2[k])))
#         print('    MedFrac: ', np.nanmedian(np.abs(vec - t2[k])/t2[k]))
#         print('    Max:     ', np.nanmax(np.abs(vec - t2[k])))
#         print('    MaxFrac: ', np.nanmax(np.abs(vec - t2[k])/t2[k]))
#     except KeyError:
#         # if np.any(np.isfinite(t2[k])):
#         notimplist += [k]

# print('Not implemented', notimplist)
