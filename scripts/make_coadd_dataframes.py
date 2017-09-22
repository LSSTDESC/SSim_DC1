#!/usr/bin/env python
"""
Convert Level 2 coadd catalogs to dask dataframes and write out to hdf.
"""
import os
import sys
import re
import logging
import numpy as np
import pandas as pd
import lsst.daf.persistence as dp
import lsst.pex.exceptions as pexExcept

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

def key_prefix(dataset, patch, band):
    r = re.compile('[,-]')
    return r.sub('_', "%s_%s_%s" % (dataset, patch, band))

def get_datasets(infile):
    return np.loadtxt(infile, dtype=str)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('repo_list', help='File with list of repos and bands.')
    parser.add_argument('--log_level', default='INFO',
                        choices='DEBUG INFO WARN ERROR CRITICAL'.split(),
                        help='Logging level. Default: "INFO"')
    parser.add_argument('--outdir', default='.', type=str,
                        help='Output directory.')
    parser.add_argument('--clobber', default=False, action='store_true',
                        help='Flag to replace existing output files. Default: False.')
    args = parser.parse_args()

    logging.basicConfig(format="%(message)s", stream=sys.stdout)
    logger = logging.getLogger()
    logger.setLevel(args.log_level)

    tract = 0
    for repo, bands in get_datasets(args.repo_list):
        if not os.path.isdir(repo):
            raise RuntimeError("%s is not a valid repository", repo)

        dataset = os.path.basename(repo)
        logging.info('processing %s', dataset)

        hdf_file = os.path.join(args.outdir, 'coadd-%s.hdf' % dataset)
        if args.clobber and os.path.isfile(hdf_file):
            os.remove(hdf_file)

        butler = dp.Butler(repo)
        patches = get_patches(butler)

        for band in bands:
            for patch in patches[tract]:
                dataId = dict(patch=patch, tract=tract, filter=band)
                try:
                    coadd = butler.get('deepCoadd', dataId=dataId)
                    catalog = butler.get('deepCoadd_meas', dataId=dataId)
                    df = make_dataframe(catalog)
                    df = add_calibrated_mags(df, catalog, coadd)
                    df.to_hdf(hdf_file, key_prefix(dataset, patch, band),
                              format='table')
                    logger.info("processed dataId: %s", dataId)
                except RuntimeError as eobj:
                    logger.debug("RuntimeError for dataId %s:", dataId)
                    logger.debug(eobj)
                    logger.debug("Skipping")
        logger.info('')
