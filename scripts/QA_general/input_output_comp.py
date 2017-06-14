import numpy as np
import matplotlib.pyplot as plt
import astropy.table
import astropy.io.fits as fits
from astropy.wcs import WCS
import pandas as pd
import pickle
from sklearn.neighbors import KDTree
from optparse import OptionParser

# Options

parser = OptionParser()

parser.add_option("--truth-table", dest="truth_tab", default=None,
                 help="Path to PhoSim instance catalog", type="string")
parser.add_option("--input-catalog", dest="input_catalog", default=None,
                 help="Path to L2 output catalog to analyze", type="string")
parser.add_option("--calexp", dest="calexp", default=None,
                 help="Path to the calibrated exposure corresponding to the \
                 catalog to be analyzed", type="string")
parser.add_option("--mag-truth", dest="mag_truth", default=None,
                 help="Path to magnitude truth table (this will be deprecated \
                 in DC1", type="string")
(o, args) = parser.parse_args()


# We read the PhoSim instance catalog as a pandas dataframe
dataframe = pd.read_table(o.truth_tab, skiprows=20, delim_whitespace=True, header=None, names=['object', 'id', 'ra','dec', 'mag_norm', 'sed_name', 'redshift', 'gamma1', 'gamma2', 'kappa', 'delta_ra', 'delta_dec', 'source_type', 'a', 'b', 'theta', 'n', 'dust_rest_name', 'A_v', 'R_v', 'dust_lab_name'])


# The next steps should be probably done with the butler but I wanted
# to allow the users work without the DM stack for now

# We read the L2 output
mytab = astropy.table.Table.read(o.input_catalog)
# We read the calibrated exposure to obtain the zeropoints
hdulist = fits.open(o.calexp)
# We get the hdu containing the image which we will use later
reference = hdulist[1]
# We get the WCS on the calibrated exposure to calculate the positions on the chip of the input sources
w = WCS(reference.header)
#xcent, ycent = w.all_world2pix(w.wcs.crval[0],w.wcs.crval[1],0.,ra_dec_order=True)
#print xcent, ycent, w.wcs.crval[0], w.wcs.crval[1], dataframe['ra'][0], dataframe['dec'][0], ' Central pixels'
x, y = w.all_world2pix(dataframe['ra']+dataframe['delta_ra'],dataframe['dec']+dataframe['delta_dec'],0.,ra_dec_order=True)
# We get the zeropoints from the calibrated exposure
zeropoint = 2.5 * np.log10(hdulist[0].header["FLUXMAG0"])
# Calculating the magnitudes
psfmag_tot = zeropoint - 2.5 * np.log10(mytab['base_PsfFlux_flux'])
# Comparing to the truth table 
catsim = astropy.table.Table.read(o.mag_truth)

# Selecting stars and galaxies using the truth table
star_sel = dataframe['sed_name'].str.contains('star')
gal_sel = dataframe['sed_name'].str.contains('galaxy')

# Some useful numbers

print 'Number of stars ', np.count_nonzero(star_sel.values)
print 'Number of galaxies ', np.count_nonzero(gal_sel.values)

# Routine to match sources using KDTree. It matches to the closest detected
# source (sometimes this is incorrect)

def match(x,y,mytab):
    """Routine that matches the truth catalog
    with the input table
    
    Args:
    ----
        x: `float` RA of the truth objects to match (in degrees)
        y: `float` dec of the truth objects to match (in degrees)
        mytab: `astropy.table.Table` table containing the L2
            input catalog.

    Returns:
    -------
        ind: `int` array of indices to select the truth objects
            that match the detected objects
    """
    X = np.zeros((len(x),2))
    X[:,0]=x
    X[:,1]=y
    tree = KDTree(X,leaf_size=40)
    Y = np.zeros((len(mytab),2))
    Y[:,0]=mytab['coord_ra']*180/np.pi
    Y[:,1]=mytab['coord_dec']*180/np.pi
    dist, ind = tree.query(Y,k=1)
    print 'Matches with distance > 1 px, ', np.count_nonzero(dist>1)
    return ind

# Routine to match sources using KDTree. It matches to the closest bright source
# it looks more accurate than match(x,y,mytab)

def match_bright(x,y,x2,y2,mags,dist=1./3600.):
    """Routine that matches the truth catalog
    with the input table
    
    Args:
    ----
        x: `float` RA of the truth objects to match (in degrees)
        y: `float` dec of the truth objects to match (in degrees)
        x2: `float` RA of detected objects to match (in degrees)
        y2: `float` dec of detected objects to match (in degrees)
        mags: `float` array containing the true input magnitudes
        dist: `float` maximum distance in degrees considered to match
            the objects, the default is 1 arcsecond.
    Returns:
    -------
        brightest_ind: `int` array of indices to select the truth objects
            that match the detected objects, returns -1 if no match has
            been found for a particular object
    """
    X = np.zeros((len(x),2))
    X[:,0]=x
    X[:,1]=y
    Y = np.zeros((len(x2),2))
    Y[:,0]=x2
    Y[:,1]=y2
    tree = KDTree(X,leaf_size=40)
    ind = tree.query_radius(Y, r=dist)
    brightest_indices = np.zeros(len(ind),dtype=np.int64)
    for i,ii in enumerate(ind):
        sorted_indices = np.argsort(mags[ii])
        if(len(sorted_indices)>0):
            brightest_indices[i] = ii[sorted_indices[0]]
        else:
            brightest_indices[i]=-1 
    return brightest_indices


# Routine to compute the ellipse parameters for each object from the input catalog
# a and b are in units of pixels!

def compute_a_b_theta_from_table(cat_table):
    Mxx = cat_table['base_SdssShape_xx']
    Myy = cat_table['base_SdssShape_yy']
    Mxy = cat_table['base_SdssShape_xy']
    Muu_p_Mvv = Mxx + Myy
    Muu_m_Mvv = np.sqrt((Mxx - Myy)**2 + 4*Mxy**2)
    Muu = 0.5*(Muu_p_Mvv + Muu_m_Mvv)
    Mvv = 0.5*(Muu_p_Mvv - Muu_m_Mvv)
    theta = 0.5*np.arctan2(2*Mxy, Mxx - Myy)*180/np.pi
    a = np.sqrt(Muu)
    b = np.sqrt(Mvv)
    return a,b,theta

# We match by distance
matched_indices = match(dataframe['ra'],dataframe['dec'],mytab)

# We match by distance and magnitude
matched_mag = match_bright(catsim['ra'],catsim['dec'],mytab['coord_ra']*180/np.pi,mytab['coord_dec']*180/np.pi,catsim['r'])

matched_gals = gal_sel.values[matched_indices]
matched_gal_indices = matched_indices[gal_sel.values[matched_indices]]
matched_stars = star_sel.values[matched_indices]
matched_star_indices = matched_indices[star_sel.values[matched_indices]]

number_gal = np.count_nonzero(matched_gals)
number_star =  np.count_nonzero(matched_stars)

# Compute the ellipse parameters
a,b,theta = compute_a_b_theta_from_table(mytab[matched_gals.flatten()])
mask_1 = a<b
aux = a
a[mask_1] = b[mask_1]
b[mask_1] = aux[mask_1]

# Plots showing the input and output ellipse parameters parameters
# The factor 0.7 comes from the PSF FWHM
plt.figure()
plt.scatter(dataframe['a'][matched_gal_indices].values*0.7,a*0.2,c=b*0.2,s=4, vmin=0., vmax=2.)
plt.plot(np.linspace(-1,6,7),np.linspace(-1,6,7),'-')
plt.colorbar(label='$b_{meas}$ [arcsec]')
plt.xlim(-.1,2.1)
plt.ylim(-.1,2.1)
plt.xlabel('$a_{true}$ [arcsec]')
plt.ylabel('$a_{meas}$ [arcsec]')
plt.savefig('bias_ellipse_a.png')
print 'Generated bias_ellipse_a.png'

plt.figure()
plt.scatter(dataframe['b'][matched_gal_indices].values*0.7,b*0.2,c=a*0.2,s=4, vmin=0, vmax=2.)
plt.plot(np.linspace(-1,6,7),np.linspace(-1,6,7),'-')
plt.ylim(-.1,2.1)
plt.xlim(-.1,2.1)
plt.xlabel('$b_{true}$ [arcsec]')
plt.ylabel('$b_{meas}$ [arcsec]')
plt.colorbar(label='$a_{meas}$ [arcsec]')
plt.savefig('bias_ellipse_b.png')
print 'Generated bias_ellipse_b.png'

# -----------------------------

# Plot checking the accuracy of extendedness as S/G classifier

plt.figure()
plt.hist(mytab['base_ClassificationExtendedness_value'][matched_gals.flatten()],weights=np.ones(number_gal)/number_gal,label='galaxies',range=(-0.05,1.05),bins=10,alpha=0.5)
plt.hist(mytab['base_ClassificationExtendedness_value'][matched_stars.flatten()],weights=np.ones(number_star)/number_star,label='stars', range=(-0.05,1.05),bins=10,alpha=0.5)
plt.xlabel('Classification Extendedness')
plt.ylabel(r'$P$ (class)')
plt.xlim(-0.2,1.4)
plt.legend()
plt.savefig('SG_classification.png')
print 'Generated SG_classification.png'
# -----------------------------

# Plot to check the redshift distribution

plt.figure()
plt.hist(dataframe['redshift'].values[matched_indices], alpha=0.5, bins=100, range=(0,5.))
plt.xlabel(r'$z$')
plt.ylabel(r'$N(z)$')
plt.savefig('redshift_hist.png')
print 'Generated redshift_hist.png'
# -----------------------------

# Plot to check the magnitude distribution

plt.figure()
plt.hist(psfmag_tot[matched_gals.flatten()], label='galaxies', range=(20,28), bins=40)
plt.savefig('mag_hist_galaxies.png')
print 'Generatd mag_hist_galaxies.png'
# ------------------------------

# Comparison of input and measured magnitudes

plt.figure()
plt.plot(catsim['r'][matched_mag[matched_mag>-1]],psfmag_tot[matched_mag>-1],',')
plt.plot(np.linspace(0,35,3),np.linspace(0,35,3),'-')
plt.xlim(14,30)
plt.ylim(14,30)
plt.xlabel(r'$r_{input,AB}$')
plt.ylabel(r'$r_{obs,AB}$')
plt.savefig('mag_scatter_one_chip.png')
print 'Generated mag_scatter_one_chip.png'
# Routine to count the number of objects in a given distance
# We are going to use this to check stellar obscuration
 
def count_close(x,y,x2,y2,distances):
    """Routine that counts the number of 
    objects that are within certain radius
    
    Args:
    ----
        x: `float` position X of objects to count
        y: `float` position Y of objects to count
        x2: `float` position X of the objects that serve as the center
            of the circle where we look for neighbors 
        y2: `float` position Y of the objects that serve as the center
            of the circle where we look for neighbors  
        distances: `float` array of radii where to count the objects
    Returns:
    -------
        neighbors: `float` the mean number of neighbors in a circle of radii
        corresponding to each entry of distances
        err: `float` standard deviation of the number of neighbors in a circle
        of radii corresponding to each entry of distances
    """
    X = np.zeros((len(x),2))
    X[:,0]=x
    X[:,1]=y
    Y = np.zeros((len(x2),2))
    Y[:,0]=x2
    Y[:,1]=y2
    tree = KDTree(X,leaf_size=40)
    neighbors = np.zeros(len(distances))
    err = np.zeros(len(distances))
    for i,distance in enumerate(distances):
        neighbors[i], err[i] = np.nanmean(tree.query_radius(Y, r=distance, count_only=True)), np.nanstd(tree.query_radius(Y, r=distance, count_only=True))
    return neighbors, err

# We are going to count how many galaxies surround each star given their flux percentile 
distances = np.linspace(1,2201,50)
stars_flux = mytab['base_SdssShape_flux'][matched_stars.flatten()]
plt.figure()
for i in range(0,5):
    flux_cuts = stars_flux > np.percentile(stars_flux[~np.isnan(stars_flux)],20*i+15)
    n, err = count_close(mytab['base_SdssCentroid_x'][matched_gals.flatten()], mytab['base_SdssCentroid_y'][matched_gals.flatten()], mytab['base_SdssCentroid_x'][matched_stars.flatten()][flux_cuts],mytab['base_SdssCentroid_y'][matched_stars.flatten()][flux_cuts],distances) 
    name='Stars in the %d-th flux percentile' % (20*i+15)
    plt.plot(distances[:-1]*0.2,(n[1:]-n[:-1])/(distances[1:]-distances[:-1]),'-',label=name)

plt.legend(loc=4)
plt.grid()
plt.xlabel(r'$\theta$ [arcsec]')
plt.ylabel(r'$\Delta\bar{N}_{gal}/\Delta\theta$')
plt.xlim(0,100)
plt.ylim(0,4.)
plt.savefig('dn_dtheta_gal_star.png')
print 'Generated dn_dtheta_gal_star.png'
flux_cuts = stars_flux > np.percentile(stars_flux[~np.isnan(stars_flux)],95)
mag_cuts = psfmag_tot > 25.5

# This probably is implemented in the DM stack as well using DS9
def plot_ref_image(xmin, xmax, ymin, ymax,reference,x,y,x2,y2,x3,y3,savename, vmin=-5, vmax=5):
    """Routine to produce plots of the image in a region from xmin
    to xmax, and ymin to ymax of the reference image
    and annotating the position of three more different catalogs
    (for example input objects, detected stars, and detected galaxies)
    
    Args:
    ----
        xmin: `float` minimum X position in the chip to be shown
        xmax: `float` maximum X position in the chip to be shown
        ymin: `float` minimum Y position in the chip to be shown
        ymax: `float` maximum Y position in the chip to be shown
        reference: `HDU` HDU containing the image to be analyzed
        x, x2, x3: `float` arrays of X positions to be marked on the image
        y, y2, y3: `float` arrays of Y positions to be marked on the image
        vmin: `float` minimum of the color scale
        vmax: `float` maximum of the color scale
    """
    
    fig, ax = plt.subplots(ncols=1,figsize=(12,12)) 
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    ax.plot(x+1,y+1,'+',c='r')
    ax.plot(x2+1,y2+1,'+',c='c')
    ax.plot(x3+1,y3+1,'x',markersize=10,c='b')
    ax.grid()
    im = ax.imshow(reference[ymin:ymax,xmin:xmax],extent=[xmin,xmax,ymin,ymax],cmap='gray',vmin=vmin,vmax=vmax, origin="lower")
    fig.colorbar(im, ax=ax, shrink=1)
    fig.savefig(savename)
# We will refer to the 
x2, y2 = w.all_world2pix(mytab['coord_ra']*180/np.pi,mytab['coord_dec']*180/np.pi,0.,ra_dec_order=True)
# Plotting around one very bright star (in the 95-th flux percentile)
plot_ref_image(x2[matched_stars.flatten()][flux_cuts][2]-100,x2[matched_stars.flatten()][flux_cuts][2]+100,y2[matched_stars.flatten()][flux_cuts][2]-100,y2[matched_stars.flatten()][flux_cuts][2]+100,reference.data,x2[matched_stars.flatten()][flux_cuts],y2[matched_stars.flatten()][flux_cuts],x2[matched_gals.flatten()],y2[matched_gals.flatten()],x,y,savename='bright_star_example.png',vmin=0,vmax=15)
print 'Generated bright_star_example.png'
# Plotting around one faint object (r > 25.5)

plt.figure()
plot_ref_image(x2[mag_cuts][1]-100,x2[mag_cuts][1]+100,y2[mag_cuts][1]-100,y2[mag_cuts][1]+100,reference.data,x2[mag_cuts],y2[mag_cuts],x2[matched_gals.flatten()],y2[matched_gals.flatten()],x,y,savename='faint_example.png')
print 'Generated faint_example.png'
print 'Done'
