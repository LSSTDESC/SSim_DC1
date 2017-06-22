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

table_name = 'StarTruth'
#projectId = 0
#projectName = 'DC1 imsim undithered'
projectId = 1
projectName = 'DC1 imsim dithered'

conn = desc.pserv.DbConnection(host='nerscdb04.nersc.gov',
                               port=3306,
                               database='DESC_DC1_Level_2')

star_truth = conn.get_pandas_data_frame('select * from StarTruth')
tree = sklearn.neighbors.KDTree(
    np.array(((star_truth['raICRS'].values*np.pi/180.,
               star_truth['decICRS'].values*np.pi/180.))).transpose())

query = '''select id, coord_ra, coord_dec,
           base_PsfFlux_flux, base_PsfFlux_fluxSigma
           from Coadd_Object where deblend_nChild=0
           and base_ClassificationExtendedness_value=0
           and projectId=%i limit 100000''' % projectId
tstart = time.time()
imsim_stars = conn.get_pandas_data_frame(query)
print('query time:', time.time() - tstart)
candidates = np.array(((imsim_stars['coord_ra'].values,
                        imsim_stars['coord_dec'].values))).transpose()

offset, index = tree.query(candidates, k=1)

# Since the KDTree is Euclidean for 2 dimensions by default, account
# for the effect of a non-zero latitude and convert from radians to
# arcsec.
offset = np.array([x[0] for x in offset])
offset *= np.cos(imsim_stars['coord_dec'].values)*180./np.pi*3600.

# Build the summary data frame.
index = (tuple([x[0] for x in index]),)
mag, mag_err = app_mag(imsim_stars['base_PsfFlux_flux'].values,
                       imsim_stars['base_PsfFlux_fluxSigma'].values)
df = pd.DataFrame(dict(objectId=imsim_stars['id'].values,
                       coord_ra=imsim_stars['coord_ra']*180./np.pi,
                       coord_dec=imsim_stars['coord_dec']*180./np.pi,
                       r_mag=mag, r_mag_err=mag_err,
                       coord_ra_true=star_truth['raICRS'].values[index],
                       coord_dec_true=star_truth['decICRS'].values[index],
                       r_mag_true=star_truth['r_mag'].values[index],
                       offset=offset))
plt.figure()
plt.hist(df['offset'].values, range=(0, 0.5), bins=100, histtype='step')
plt.xlabel('offset (arcsec)')
plt.title(projectName)
plt.savefig('%s_star_offsets.png' % projectName.replace(' ', '_'))

df_match = df.query("offset < 0.2")

plt.figure()
plt.errorbar(df_match['r_mag_true'], df_match['r_mag'],
             yerr=df_match['r_mag_err'], fmt='.')
plt.errorbar(df_match['r_mag_true'], df_match['r_mag'], fmt='.')
plt.xlabel('r mag (true)')
plt.ylabel('r mag (meas)')
plt.title(projectName)
axis = list(plt.axis())
axis[2:] = axis[:2]
plt.axis(axis)
plt.savefig('%s_r_mag_comparison.png' % projectName.replace(' ', '_'))
