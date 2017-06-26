import os
import pickle
import numpy as np
from lsst.afw.coord.coordLib import IcrsCoord
import lsst.afw.fits
import lsst.afw.geom as afw_geom
import matplotlib.pyplot as plt

plt.ion()

# Nominal field of view radius in radians (= 1.7475 deg)
_fov_radius = 0.0305

class DC1Regions(object):
    def __init__(self, wcs, fov_radius=_fov_radius):
        self.wcs = wcs
        self.fov_radius = fov_radius
        self.ra = np.array([1.625039, 1.679307, 1.58809, 1.641324])
        self.dec = np.array([-0.523011, -0.513259, -0.4789, -0.46990])
        self.field_centers \
            = [IcrsCoord(afw_geom.Angle(ra), afw_geom.Angle(dec)) for
               ra, dec in zip(self.ra, self.dec)]

    def plot_regions(self, radius=None, ls='-', color='black'):
        if radius is None:
            radius = self.fov_radius*180./np.pi
        phi = np.linspace(0, 2.*np.pi, 100)
        for ra, dec in zip(self.ra, self.dec):
            x, y = (radius*np.sin(phi) + ra*180./np.pi,
                    radius*np.cos(phi) + dec*180./np.pi)
            plt.plot(x, y, color='black', ls=ls)

    def eff_radius(self, patch):
        bbox = patch.getOuterBBox()
        half_side = self.wcs.pixelScale().asRadians()*max(bbox.getHeight(),
                                                          bbox.getWidth())/2.
        return np.sqrt(self.fov_radius**2 + half_side**2)

    def patch_vertices(self, patch):
        bbox = patch.getOuterBBox()
        vertices = [(bbox.getMinX(), bbox.getMinY()),
                    (bbox.getMinX(), bbox.getMaxY()),
                    (bbox.getMaxX(), bbox.getMinY()),
                    (bbox.getMaxX(), bbox.getMaxY())]
        return [self.wcs.pixelToSky(*vertex) for vertex in vertices]

    def include_patch(self, patch):
        radius = self.eff_radius(patch)
        vertices = self.patch_vertices(patch)
        for vertex in vertices:
            for field_dir in self.field_centers:
                if vertex.angularSeparation(field_dir).asRadians() < radius:
                    return True
        return False


def plot_tract_patches(sim, repo, filter='r', apply_dc1_filter=False):
    skymap = pickle.load(open(os.path.join(repo, 'deepCoadd', 'skyMap.pickle')))
    tract = skymap[0]
    tract_label = '%i' % tract.getId()
    wcs = tract.getWcs()

    dc1_regions = DC1Regions(wcs)

    # plot patch centers with objects.
    ra, dec = [], []
    patches = []
    for patch in tract:
        if not os.path.isfile(os.path.join(repo, 'deepCoadd', filter,
                                           tract_label,
                                           '%s,%s.fits' % patch.getIndex())):
            continue
        if apply_dc1_filter and not dc1_regions.include_patch(patch):
            continue
        patches.append(patch)
        bbox = patch.getInnerBBox()
        x = (bbox.getMinX() + bbox.getMaxX())/2.
        y = (bbox.getMinY() + bbox.getMaxY())/2.
        sky_coord = wcs.pixelToSky(x, y)
        ra.append(sky_coord.getLongitude().asDegrees())
        dec.append(sky_coord.getLatitude().asDegrees())
    plt.errorbar(ra, dec, fmt='.', label='%s' % sim)

    # plot tract boundary
    vertex_list = [(vtx.getLongitude().asDegrees(),
                    vtx.getLatitude().asDegrees())
                   for vtx in tract.getVertexList()]
    vertex_list.append(vertex_list[0])
    x, y = zip(*vertex_list)
    plt.plot(x, y, label='%s tract' % sim)
    plt.xlabel('RA')
    plt.ylabel('Dec')
    return tract, patches

if __name__ == '__main__':
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    plot_tract_patches('phosim',
                       '/global/cscratch1/sd/descdm/DC1/DC1-phoSim-3a')
    plot_tract_patches('imsim-dith',
                       '/global/cscratch1/sd/descdm/DC1/DC1-imsim-dithered')
    plot_tract_patches('imsim-undith',
                       '/global/cscratch1/sd/descdm/DC1/full_focalplane_undithered')
    tract, patches \
        = plot_tract_patches('phosim_trimmed',
                             '/global/cscratch1/sd/descdm/DC1/DC1-phoSim-3a',
                             apply_dc1_filter=True)

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1.01), fontsize='x-small')

    dc1_regions = DC1Regions(tract.getWcs())
    dc1_regions.plot_regions()

    with open('phosim_DC1_region_patches.txt', 'w') as output:
        for patch in patches:
            output.write('%s,%s\n' % patch.getIndex())
