from __future__ import with_statement
import argparse
import os
import numpy as np


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Generate an InstanceCatalog')
    parser.add_argument('--db', type=str,
                        default='minion_1016_sqlite_new_dithers.db',
                        help='path to the OpSim database to query')
    parser.add_argument('--out', type=str,
                        default='.',
                        help='directory where output will be written')
    parser.add_argument('--id', type=int, nargs='+',
                        default=None,
                        help='obsHistID to generate InstanceCatalog for (a list)')
    parser.add_argument('--dither', type=str,
                        default='True',
                        help='whether or not to apply dithering (true/false; default true)')
    args = parser.parse_args()

    obshistid_list = args.id
    opsimdb = args.db
    out_dir = args.out
    dither_switch = True
    if args.dither.lower()[0] == 'f':
        dither_switch = False

    from lsst.sims.catUtils.utils import ObservationMetaDataGenerator

    if not os.path.exists(opsimdb):
        raise RuntimeError('%s does not exist' % opsimdb)

    obs_generator = ObservationMetaDataGenerator(database=opsimdb, driver='sqlite')

    from lsst.sims.catUtils.exampleCatalogDefinitions import PhoSimCatalogPoint
    from lsst.sims.catUtils.exampleCatalogDefinitions import PhoSimCatalogSersic2D
    from lsst.sims.catUtils.exampleCatalogDefinitions import PhoSimCatalogZPoint
    from lsst.sims.catUtils.exampleCatalogDefinitions import DefaultPhoSimHeaderMap
    from lsst.sims.catUtils.baseCatalogModels import StarObj
    from lsst.sims.catUtils.baseCatalogModels import GalaxyBulgeObj, GalaxyDiskObj
    from lsst.sims.catUtils.baseCatalogModels import GalaxyAgnObj
    from lsst.sims.utils import _getRotSkyPos
    import copy

    star_db = StarObj(database='LSSTCATSIM', host='fatboy.phys.washington.edu',
                      port=1433, driver='mssql+pymssql')

    bulge_db = GalaxyBulgeObj(connection=star_db.connection)
    disk_db = GalaxyDiskObj(connection=star_db.connection)
    agn_db = GalaxyAgnObj(connection=star_db.connection)

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    phosim_header_map = copy.deepcopy(DefaultPhoSimHeaderMap)
    phosim_header_map['nsnap'] = 1
    phosim_header_map['vistime'] = 30.0
    phosim_header_map['camconfig'] = 1
    for obshistid in obshistid_list:

        obs_list = obs_generator.getObservationMetaData(obsHistID=obshistid,
                                                        boundType='circle',
                                                        boundLength=2.0)

        obs = obs_list[0]
        if dither_switch:
            print 'dithering'
            obs.pointingRA = np.degrees(obs.OpsimMetaData['randomDitherFieldPerVisitRA'])
            obs.pointingDec = np.degrees(obs.OpsimMetaData['randomDitherFieldPerVisitDec'])
            rotSky = _getRotSkyPos(obs._pointingRA, obs._pointingDec, obs,
                                   obs.OpsimMetaData['ditheredRotTelPos'])

            obs.rotSkyPos = np.degrees(rotSky)
            obs.OpsimMetaData['rotTelPos'] = obs.OpsimMetaData['ditheredRotTelPos']

        cat_name = os.path.join(out_dir,'phosim_cat_%d.txt' % obshistid)
        star_name = os.path.join('star_cat_%d.txt' % obshistid)
        gal_name = os.path.join('gal_cat_%d.txt' % obshistid)
        agn_name = os.path.join('agn_cat_%d.txt' % obshistid)

        cat = PhoSimCatalogPoint(star_db, obs_metadata=obs)

        cat.phoSimHeaderMap = phosim_header_map
        with open(cat_name, 'w') as output:
            cat.write_header(output)
            output.write('includeobj %s\n' % star_name)
            output.write('includeobj %s\n' % gal_name)
            output.write('includeobj %s\n' % agn_name)

        cat.write_catalog(os.path.join(out_dir, star_name), write_header=False,
                          chunk_size=100000)

        cat = PhoSimCatalogSersic2D(bulge_db, obs_metadata=obs)
        cat.write_catalog(os.path.join(out_dir, gal_name), write_header=False,
                          chunk_size=100000)
        cat = PhoSimCatalogSersic2D(disk_db, obs_metadata=obs)
        cat.write_catalog(os.path.join(out_dir, gal_name), write_header=False,
                          write_mode='a', chunk_size=100000)

        cat = PhoSimCatalogZPoint(agn_db, obs_metadata=obs)
        cat.write_catalog(os.path.join(out_dir, agn_name), write_header=False,
                          chunk_size=100000)
