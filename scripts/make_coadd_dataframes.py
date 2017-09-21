"""
Convert Level 2 coadd catalogs to dask dataframes and write out to hdf.
"""
import sys
import re
import logging
import numpy as np
import pandas as pd
import lsst.daf.persistence as dp
import lsst.pex.exceptions as pexExcept

logging.basicConfig(format="%(message)s", stream=sys.stdout)
logger = logging.getLogger()
logger.setLevel(logging.INFO)
#logger.setLevel(logging.DEBUG)

def make_dataframe(catalog):
    """
    Create a data frame from a coadd catalog using the coadd catalog
    schema.
    """
    df = pd.DataFrame()
    for col in catalog.getSchema():
        name = col.field.getName()
        df[name] = catalog.get(name)
    return df

def add_calibrated_mags(df, catalog, coadd):
    """
    Add columns for calibrated magnitudes for psfFlux and cmodelFlux.
    """
    modelFlux = catalog.getModelFlux()
    psfFlux = catalog.getPsfFlux()
    coadd_calib = coadd.getCalib()
    modelMags = []
    psfMags = []
    for mflux, pflux in zip(modelFlux, psfFlux):
        try:
            modelMags.append(coadd_calib.getMagnitude(mflux))
        except pexExcept.DomainError:
            modelMags.append(np.nan)
        try:
            psfMags.append(coadd_calib.getMagnitude(pflux))
        except pexExcept.DomainError:
            psfMags.append(np.nan)
    df['cmodelMag'] = np.array(modelMags)
    df['psfMag'] = np.array(psfMags)
    return df

def get_patches(butler):
    """
    Get all of the patches in a discrete skymap.
    """
    skymap = butler.get('deepCoadd_skyMap')
    tract = [x for x in skymap][0]
    patches = dict()
    for tract in skymap:
        patches[tract.getId()] = ['%i,%i' % x.getIndex() for x in tract]
    return patches

dataset = 'DC1-imsim-dithered'

def key_prefix(patch, dataset=dataset):
    r = re.compile('[,-]')
    return r.sub('_', "%s_%s" % (dataset, patch))

hdf_file = 'coadd-%s.hdf' % dataset
repo = '/global/cscratch1/sd/descdm/DC1/rerun/%s' % dataset

butler = dp.Butler(repo)
patches = get_patches(butler)

tract = 0
band = 'r'

for patch in patches[tract]:
    dataId = dict(patch=patch, tract=tract, filter=band)
    try:
        coadd = butler.get('deepCoadd', dataId=dataId)
        catalog = butler.get('deepCoadd_meas', dataId=dataId)
        logging.info("processing dataId: %s", dataId)
        df = make_dataframe(catalog)
        df = add_calibrated_mags(df, catalog, coadd)
        df.to_hdf(hdf_file, key_prefix(patch), format='table')
    except RuntimeError as eobj:
        logger.debug("RuntimeError for patch %s:", patch)
        logger.debug(eobj)
        logger.debug("Skipping")
