import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import treecorr as tc
import astropy.table
from optparse import OptionParser
from astropy.cosmology import FlatLambdaCDM

#Options
parser = OptionParser()
parser.add_option("--truth-file", dest="truth_tab", default=None,
                 help="Path to PhoSim instance catalog", type="string")
parser.add_option("--plot-map", dest="plot_map", default=None, action="store_true",
                 help="Adding this option will generate a 2D histogram of the \
                 input file")
parser.add_option("--random-multiplier", dest="f_random", default=10,
                 help="Ratio of random points to input galaxies", type="float")

(o, args) = parser.parse_args()
 # Input cosmology of the simulation
omega_matter = 0.25
Omega_baryon = 0.045
Omega_curvature = 0
H0 = 73
sigma_8 = 0.9
n_s = 1
cosmo=FlatLambdaCDM(H0=H0,Om0=omega_matter)
# PhoSim instance catalog ("truth table")
dataframe = pd.read_table(o.truth_tab, skiprows=20, delim_whitespace=True, header=None, names=['object', 'id', 'ra','dec', 'mag_norm', 'sed_name', 'redshift', 'gamma1', 'gamma2', 'kappa', 'delta_ra', 'delta_dec', 'source_type', 'dust_rest_name', 'A_v', 'R_v', 'dust_lab_name'])
# We select the stars and galaxies reading the name of the SED
star_sel = dataframe['sed_name'].str.contains('star')
gal_sel = dataframe['sed_name'].str.contains('galaxy')
# We generate an array with ra, dec, z and the comoving distance for galaxies
redshift = dataframe['redshift'][gal_sel].values
ra = dataframe['ra'][gal_sel].values
dec = dataframe['dec'][gal_sel].values
rad=cosmo.comoving_distance(redshift)
# Generating the QA 2d plot of the map
if o.plot_map:
    plt.hist2d(ra,np.sin(dec*np.pi/180.),bins=25)
    plt.xlabel('RA [deg]')
    plt.ylabel('$sin$(dec)')
    plt.colorbar()
    plt.savefig('input_galaxies.png')
# Converting to a cartesian grid
catx=rad*np.cos(dec/180*np.pi)*np.cos(ra/180*np.pi)
caty=rad*np.cos(dec/180*np.pi)*np.sin(ra/180*np.pi)
catz=rad*np.sin(dec/180*np.pi)

# Coordinates of the center of the FOV (also present in PhoSim instance catalog)
cent_ra = 53.0091385
cent_dec = -27.4389488
radius=2.5

# Routine to generate randoms to compute the 2pcf

def generate_rnd(factor=10):
    #Creating random that follows N(z)
    r_rnd = np.random.choice(rad,size=factor*len(ra))
    ra_rnd = cent_ra-radius+2*radius*np.random.random(size=factor*len(ra))
    cth_rnd = np.sin(cent_dec*np.pi/180-radius*np.pi/180)+(np.sin(cent_dec*np.pi/180+radius*np.pi/180)-np.sin(cent_dec*np.pi/180-radius*np.pi/180))*np.random.random(size=factor*len(ra))
    dec_rnd = np.arcsin(cth_rnd)*180/np.pi
    myrand = (ra_rnd-cent_ra)**2+(dec_rnd-cent_dec)**2<=radius**2
    ra_rnd = ra_rnd[myrand]
    dec_rnd = dec_rnd[myrand]
    r_rnd = r_rnd[myrand]
    xrand=r_rnd*np.cos(ra_rnd*np.pi/180.)*np.cos(dec_rnd*np.pi/180.)
    yrand=r_rnd*np.sin(ra_rnd*np.pi/180.)*np.cos(dec_rnd*np.pi/180.)
    zrand=r_rnd*np.sin(dec_rnd*np.pi/180.)
    return xrand, yrand, zrand

xrand, yrand, zrand = generate_rnd(factor=o.f_random)
print 'Generated ', len(xrand), ' randoms'
xires=[]
icount=0

# Computing the correlation function in these bins, for now hardcoded

for zmin,zmax in [(0,0.2),(0.2,0.4),(0.4,0.6),(0.6,0.8), (0.8,1.0),(1.0,1.2),(1.2,1.4),(1.4,1.6),(1.6,1.8),(1.8,2.0),(2.0,5.0)]:
    print "Correlating",zmin,"-",zmax
    dimin,dimax=cosmo.comoving_distance([zmin,zmax])
    print dimin, dimax
    wh=np.logical_and(redshift>zmin,redshift<zmax)
    print 'Bin contains : ', np.count_nonzero(wh), ' galaxies'
    sigcat=tc.Catalog(x=catx[wh],y=caty[wh],z=catz[wh])
    whr=np.logical_and(xrand**2+yrand**2+zrand**2>dimin.value**2, xrand**2+yrand**2+zrand**2<dimax.value**2)
    print 'Bin contains : ', np.count_nonzero(whr), ' random objects'
    rancat=tc.Catalog(x=xrand[whr],y=yrand[whr],z=zrand[whr])
    dd=tc.NNCorrelation(min_sep=0.1,bin_size=0.1,max_sep=200.)
    dd.process(sigcat)
    dr=tc.NNCorrelation(min_sep=0.1,bin_size=0.1,max_sep=200.)
    dr.process(sigcat,rancat)
    rr=tc.NNCorrelation(min_sep=0.1,bin_size=0.1,max_sep=200.)
    rr.process(rancat,rancat)
    xi,xivar=dd.calculateXi(rr,dr)
    rdist=np.exp(dd.logr)
    tab = astropy.table.Table([rdist,xi,xivar],names=('r','xi','xivar'))
    filename = 'xi_r_%d.fits.gz' % (icount)
    tab.write(filename)
    icount=icount+1

