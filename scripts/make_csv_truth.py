import os
import sys
import csv
from desc.imsim import parsePhoSimInstanceFile as catalog_parser
from desc.imsim import imsim_truth

def imsim_star_truth_csv_file(object_cat, nrows, skiprows, csv_file):
    commands, objects = catalog_parser(object_cat, numRows=nrows,
                                       skiprows=skiprows, check_commands=False)
    objects = objects.query("galSimType == 'pointSource'")
    seds = [x[0] for x in objects.groupby('sedFilepath').sedFilepath.unique()]
    columns = 'id raICRS decICRS u_mag g_mag r_mag i_mag z_mag y_mag'.split()
    with open(csv_file, 'w') as output:
        writer = csv.writer(output, delimiter=',', lineterminator='\n')
        writer.writerow(columns)
        for sed in seds:
            app_mags = imsim_truth.ApparentMagnitudes(sed)
            my_objects = objects.query("sedFilepath == '%s'" % sed)
            for i in range(len(my_objects)):
                obj = my_objects.iloc[i]
                mags = app_mags(obj)
                row = [long(obj.uniqueId), obj.raJ2000, obj.decJ2000]
                row.extend(mags.values())
                writer.writerow(row)

def imsim_galaxy_truth_csv_file(object_cat, nrows, skiprows, csv_file):
    commands, objects = catalog_parser(object_cat, numRows=nrows,
                                       skiprows=skiprows, check_commands=False)
    objects = objects.query("galSimType == 'sersic'")
    seds = [x[0] for x in objects.groupby('sedFilepath').sedFilepath.unique()]
    columns = '''id raICRS decICRS u_mag g_mag r_mag i_mag z_mag y_mag
                 redshift majorAxis minorAxis positionAngle sindex'''.split()
    with open(csv_file, 'w') as output:
        writer = csv.writer(output, delimiter=',', lineterminator='\n')
        writer.writerow(columns)
        for sed in seds:
            app_mags = imsim_truth.ApparentMagnitudes(sed)
            my_objects = objects.query("sedFilepath == '%s'" % sed)
            for i in range(len(my_objects)):
                obj = my_objects.iloc[i]
                mags = app_mags(obj)
                row = [long(obj.uniqueId), obj.raJ2000, obj.decJ2000]
                row.extend(mags.values())
                row.extend([obj.redshift, obj.majorAxis, obj.minorAxis,
                            obj.positionAngle, obj.sindex])
                writer.writerow(row)

funcs = dict(star=imsim_star_truth_csv_file,
             galaxy=imsim_galaxy_truth_csv_file)

object_catalog = sys.argv[1]
object_type = sys.argv[2]
row_min = int(sys.argv[3])
row_max = int(sys.argv[4])

visit = os.path.basename(object_catalog).split('_')[1]

csv_file = '%s_%s_%08i_%08i.csv' % (object_type, visit, row_min, row_max)

funcs[object_type](object_catalog, row_max - row_min, row_min, csv_file)
