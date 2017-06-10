"""
Script to ingest DC1 Level 2 coadd catalogs into a MySQL database.
"""
from __future__ import print_function
import os
import sys
import glob
import time
import argparse
import desc.pserv
import desc.pserv.utils as pserv_utils

def add_project(connection, projectId, projectName, table_name='Project',
                dry_run=True):
    "Insert a project into the Project table."
    sql = "insert into %s (projectId, projectName) values (%i, '%s')" \
        % (table_name, projectId, projectName)
    if dry_run:
        print(sql)
    else:
        connection.apply(sql)

#db_info = dict(host='nerscdb04.nersc.gov',
#               database='DESC_DC1_Level_2')
db_info = dict(host='scidb1.nersc.gov',
               database='DESC_Twinkles_Level_2')

projectId = 0
projectName = 'DC1 imsim undithered'

dry_run = False
#dry_run = True

connection = desc.pserv.DbConnection(**db_info)

pserv_utils.create_table(connection, 'Project', dry_run=dry_run)
add_project(connection, projectId, projectName, dry_run=dry_run)

output_repo = '/global/cscratch1/sd/descdm/DC1/full_focalplane_undithered'
coadd_files = sorted(glob.glob(os.path.join(output_repo, 'deepCoadd-results',
                                            'merged', '0', '*', 'ref-*.fits')))

hdunum = 1
table_name = 'Coadd_Object'
sql_file = 'create_%s.sql' % table_name
desc.pserv.create_schema_from_fits(coadd_files[0], hdunum, sql_file, table_name,
                                   primary_key='id, projectId, patch',
                                   add_columns=('projectId INT',
                                                'patch char(20)'))
connection.run_script(sql_file, dry_run=dry_run)

added_columns = dict(projectId=projectId)
for i, fits_file in enumerate(coadd_files):
    added_columns['patch'] \
        = fits_file[:-len('.fits')].split('-')[-1].replace(',', '\,')
    csv_file = os.path.basename(fits_file).replace('.fits', '.csv')
    print("% 3i: %s - processing %s" % (i, time.asctime(), csv_file))
    sys.stdout.flush()
    if not dry_run:
        desc.pserv.create_csv_file_from_fits(fits_file, hdunum, csv_file,
                                             added_columns=added_columns)
        connection.load_csv(table_name, csv_file)
        os.remove(csv_file)
