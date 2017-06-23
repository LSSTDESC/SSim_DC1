import sys
import time
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sklearn.neighbors
import desc.pserv
import desc.pserv.utils as pserv_utils
plt.ion()

def app_mag(flux, flux_err):
    # From butler.get('deepCoaad').getCalib().getFluxMag0():
    flux_mag0 = 63095734448.0194
    return -2.5*np.log10(flux/flux_mag0), 2.5*np.log10(flux_mag0)*flux_err/flux

table_name = 'GalaxyTruth'
#projectId = 0
#projectName = 'DC1 imsim undithered'
projectId = 1
projectName = 'DC1 imsim dithered'

conn = desc.pserv.DbConnection(host='nerscdb04.nersc.gov',
                               port=3306,
                               database='DESC_DC1_Level_2')

query = '''select id, coord_ra, coord_dec,
           modelfit_CModel_flux, modelfit_CModel_fluxSigma
           from Coadd_Object where deblend_nChild=0
           and base_ClassificationExtendedness_value=1
           and projectId=%i and patch="'10,10'" limit 300000''' % projectId
print(query)
tstart = time.time()
imsim_galaxies = conn.get_pandas_data_frame(query)
print('Coad_Object query time:', time.time() - tstart)
sys.stdout.flush()
imsim_galaxies.to_pickle('imsim_galaxies.pkl')

#imsim_galaxies = pd.read_pickle('imsim_galaxies.pkl')

ra_min = min(imsim_galaxies['coord_ra'])*180./np.pi
ra_max = max(imsim_galaxies['coord_ra'])*180./np.pi
dec_min = min(imsim_galaxies['coord_dec'])*180./np.pi
dec_max = max(imsim_galaxies['coord_dec'])*180./np.pi
query = '''select * from GalaxyTruth where
           raICRS > %f and raICRS < %f and decICRS > %f and decICRS < %f
        ''' % (ra_min, ra_max, dec_min, dec_max)
print(query)
tstart = time.time()
galaxy_truth = conn.get_pandas_data_frame(query)
print("galaxy_truth query time:", time.time() - tstart)
sys.stdout.flush()
galaxy_truth.to_pickle('galaxy_truth.pkl')

#galaxy_truth = pd.read_pickle('galaxy_truth.pkl')

tree = sklearn.neighbors.KDTree(
    np.array(((galaxy_truth['raICRS'].values*np.pi/180.,
               galaxy_truth['decICRS'].values*np.pi/180.))).transpose())


candidates = np.array(((imsim_galaxies['coord_ra'].values,
                        imsim_galaxies['coord_dec'].values))).transpose()

offset, index = tree.query(candidates, k=1)

# Since the KDTree is Euclidean for 2 dimensions by default, account
# for the effect of a non-zero latitude and convert from radians to
# arcsec.
offset = np.array([x[0] for x in offset])
offset *= np.cos(imsim_galaxies['coord_dec'].values)*180./np.pi*3600.

# Build the summary data frame.
index = (tuple([x[0] for x in index]),)
mag, mag_err = app_mag(imsim_galaxies['modelfit_CModel_flux'].values,
                       imsim_galaxies['modelfit_CModel_fluxSigma'].values)
df = pd.DataFrame(dict(objectId=imsim_galaxies['id'].values,
                       coord_ra=imsim_galaxies['coord_ra']*180./np.pi,
                       coord_dec=imsim_galaxies['coord_dec']*180./np.pi,
                       r_mag=mag, r_mag_err=mag_err,
                       coord_ra_true=galaxy_truth['raICRS'].values[index],
                       coord_dec_true=galaxy_truth['decICRS'].values[index],
                       r_mag_true=galaxy_truth['r_mag'].values[index],
                       offset=offset))
plt.figure()
plt.hist(df['offset'].values, range=(0, 0.5), bins=100, histtype='step')
plt.xlabel('offset (arcsec)')
plt.title(projectName + ', galaxies, patch 10,10')
plt.savefig('%s_galaxy_offsets.png' % projectName.replace(' ', '_'))

df_match = df.query("offset < 0.2")

plt.figure()
plt.hist(df['r_mag'], bins=100, range=(14, 40), label='r mag (meas)',
         histtype='step', color='blue')
plt.hist(df['r_mag_true'], bins=100, range=(14, 40), label='r mag (true)',
         histtype='step', color='green')
plt.hist(df_match['r_mag'], bins=100, range=(14, 40),
         label='r mag (meas), offset<0.2"', linestyle='dashed', color='blue',
         histtype='step')
plt.hist(df_match['r_mag_true'], bins=100, range=(14, 40),
         linestyle='dashed', color='green',
         label='r mag (true), offset<0.2"', histtype='step')
plt.xlabel('r mag')
plt.ylabel('entries / bin')
plt.title(projectName + ', galaxies, patch 10,10')
plt.legend(loc=0, fontsize='small')
plt.savefig('%s_galaxy_r_mag_dists.png' % projectName.replace(' ', '_'))

plt.figure()
#plt.errorbar(df_match['r_mag_true'], df_match['r_mag'],
#             yerr=df_match['r_mag_err'], fmt='.')
plt.errorbar(df_match['r_mag_true'], df_match['r_mag'], fmt='.')
plt.xlabel('r mag (true)')
plt.ylabel('r mag (meas)')
plt.title(projectName + ', galaxies, patch 10,10, (offset<0.2")')
axis = list(plt.axis())
axis[2:] = axis[:2]
plt.axis(axis)
plt.savefig('%s_r_mag_galaxy_comparison.png' % projectName.replace(' ', '_'))
