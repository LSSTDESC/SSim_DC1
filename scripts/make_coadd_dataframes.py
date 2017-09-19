"""
Convert Level 2 coadd catalogs to dask dataframes and write out to hdf.
"""
from __future__ import print_function
import pandas as pd
import lsst.daf.persistence as dp

def make_dataframe(catalog):
    df = pd.DataFrame()
    for col in catalog.getSchema():
        name = col.field.getName()
        df[name] = catalog.get(name)
    return df

def get_patches(butler):
    skymap = butler.get('deepCoadd_skyMap')
    tract = [x for x in skymap][0]
    patches = dict()
    for tract in skymap:
        patches[tract.getId()] = ['%i,%i' % x.getIndex() for x in tract]
    return patches

dataset = 'DC1-imsim-dithered'

key_prefix = lambda i: "%s-%05i" % (dataset, i)
hdf_file = 'coadd-%s.hdf' % dataset
repo = '/global/cscratch1/sd/descdm/DC1/rerun/%s' % dataset

butler = dp.Butler(repo)
patches = get_patches(butler)

tract = 0
band = 'r'

i = -1
for patch in patches[tract]:
    dataId = dict(patch=patch, tract=tract, filter=band)
    try:
        catalog = butler.get('deepCoadd_meas', dataId=dataId)
        i += 1
    except RuntimeError as eobj:
        print("RuntimeError for patch %s:" % patch)
        print(eobj)
        print("skipping")
        continue

    df = make_dataframe(catalog)
    df.to_hdf(hdf_file, key_prefix(i), format='table')
