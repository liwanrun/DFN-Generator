import os
import shapely
import numpy as np
import shapely.geometry
from fracture import SmoothFracture, RoughFracture
from bridson_sampling import Bridson_sampling
from space_partition import CellPartition

class RandomDFN:
    '''Randomly distributed discrete fracture network'''
    def __init__(self, doi=(10.0, 10.0), num=100, msh=None) -> None:
        self.doi = doi
        self.num = num
        self.msh = msh
        self.fractures = []
        self.frac_bbox = []
        self._minimum_gap = 0.0
        self._guess_times = 0

    @property
    def minimum_gap(self):
        return self._minimum_gap
    
    @minimum_gap.setter
    def minimum_gap(self, value):
        if not isinstance(value, float) or value <= 0:
            raise ValueError("positive number expected")
        self._minimum_gap = value

    @property
    def guess_times(self):
        return self._guess_times
    
    @guess_times.setter
    def guess_times(self, value):
        if not isinstance(value, int) or value <= 0:
            raise ValueError("positive number expected")
        self._guess_times = value

    def poisson_sample(self, dist=1.0):
        loc = Bridson_sampling(self.doi[0], self.doi[1], dist)
        self.num = len(loc)
        return loc
    
    def random_sample(self, num):
        loc = np.random.rand(num, 2) * self.doi
        self.num = len(loc)
        return loc

    # Rough fractures
    def generate_RDFN(self, loc, dip, scf, amp=(0.0218, 0.0111, 0.0019)):
        self.fractures = []
        for i in range(self.num):
            frac = RoughFracture(d238=amp)
            frac.scale(scf[i])
            frac.rotate(dip[i])
            frac.translate(loc[i])
            self.fractures.append(frac)

    # Smooth fractures
    def generate_SDFN(self, loc, dip, scf):
        self.fractures = []
        for i in range(self.num):
            frac = SmoothFracture(npts=2)
            frac.Points[0][0] = loc[i][0] - 0.5*np.cos(dip[i]*np.pi/180)*scf[i]
            frac.Points[0][1] = loc[i][1] - 0.5*np.sin(dip[i]*np.pi/180)*scf[i]
            frac.Points[1][0] = loc[i][0] + 0.5*np.cos(dip[i]*np.pi/180)*scf[i]
            frac.Points[1][1] = loc[i][1] + 0.5*np.sin(dip[i]*np.pi/180)*scf[i]
            self.fractures.append(frac)

    # Fracture-by-fracture
    def has_collision(self, frac):
        bbox = shapely.LineString(frac.Points).buffer(0.5*self._minimum_gap)
        pt1 = shapely.Point(frac.Points[0])
        pt2 = shapely.Point(frac.Points[-1])
        cells = self.msh.get_incident_cells(bbox)

        checked = []
        for bid in cells:
            i, j = bid
            points = self.msh.Points[i][j]
            for p in points:
                if bbox.contains(p):
                    return True
                
            curves = self.msh.Curves[i][j]
            for c in curves:
                if c not in checked:
                    if c.contains(pt1) or c.contains(pt2):
                        return True
                checked.append(c)
        return False

    def generate_FBYF(self, dip, scf):
        domain_box = shapely.geometry.box(0.0, 0.0, self.doi[0], self.doi[1])
        shrink_box = shapely.geometry.box(self._minimum_gap, self._minimum_gap, 
                                          self.doi[0]-self._minimum_gap * 0.5, 
                                          self.doi[1]-self._minimum_gap * 0.5)
        for _ in range(self._guess_times):
            fracture = RoughFracture(d238=(0.0118, 0.0000, 0.0))
            fracture.scale(np.random.uniform(scf[0], scf[1]))
            fracture.rotate(np.random.uniform(dip[0], dip[1]))

            max_iters = 100
            for _ in range(max_iters):
                fracture.translate(np.random.uniform(0.0, self.doi[0], 2))
                if not shrink_box.contains(shapely.LineString(fracture.Points)):
                    continue

                if not self.has_collision(fracture):
                    line = shapely.LineString(fracture.Points)
                    bbox = line.buffer(0.5 * self._minimum_gap)
                    self.frac_bbox.append(bbox)     # optional
                    self.fractures.append(fracture)
                    self.msh.assign_curve_to_cells(bbox)

                    point = shapely.Point(fracture.Points[0])
                    if domain_box.contains(point):
                        self.msh.assign_point_to_cell(point)

                    point = shapely.Point(fracture.Points[-1])
                    if domain_box.contains(point):
                        self.msh.assign_point_to_cell(point)
                    break

    def visualize(self, axs):
        for frac in self.fractures:
            x = frac.Points[:, 0]
            y = frac.Points[:, 1]
            axs.plot(x, y, lw=1.5)

    def write_to_gmsh_model(self, fname):
        suffix = os.path.splitext(os.path.basename(fname))[1]
        if suffix != '.geo':
            print('Invalid file format!')
            return
        
        with open(fname, 'wt') as file:
            file.write('SetFactory("OpenCASCADE");\n')
            file.write(f'lc = {self._minimum_gap};\n')
            point_off = 1
            curve_off = 1
            surface_off = 1
            # write DFNs
            file.write('/* Add discrete fracture network (DFN) */\n')
            for i, frac in enumerate(self.fractures):
                for j, vert in enumerate(frac.Points):
                    file.write(f'Point({j + point_off}) = {{{vert[0]}, {vert[1]}, 0.0, 1.0}};\n')
                file.write(f'BSpline({i + curve_off}) = {{{point_off}:{point_off + len(frac.Points)-1}}};\n')
                point_off = point_off + len(frac.Points)
            # write domain
            file.write('/* Add domain of interest (DOI) */\n')
            file.write(f'Rectangle({surface_off}) = {{0.0, 0.0, 0.0, {self.doi[0]}, {self.doi[1]}, 0.0}};\n')
            file.write(f'BooleanFragments {{ Surface{{:}}; Delete; }}{{ Curve{{:}}; Delete; }}\n')
            # physical group
            file.write('/* Add physical groups */\n')
            file.write(f'Physical Surface("block-1", newreg) = {{ Surface{{:}} }};\n')
            file.write(f'Physical Curve("constraint-1", newreg) = {{ Curve{{:}} }};\n')
            # mesh size
            tol = 1.0e-06
            file.write('/* Set mesh size */\n')
            file.write(f'dfnPoints[] = Point In BoundingBox {{{tol}, {tol}, 0.0, {self.doi[0] - tol}, {self.doi[1] - tol}, 0.0}};\n') 
            file.write(f'bndPoints[] = Point{{:}};\n')
            file.write(f'bndPoints[] -= dfnPoints[];\n')
            file.write(f'MeshSize {{ dfnPoints[] }} = lc;\n')
            file.write(f'MeshSize {{ bndPoints[] }} = 2*lc;\n')

# Unit test
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon

    # boundary
    w = 10.0; h = 10.0
    bbox = np.array([[0.0, 0.0], [w, 0], [w, h], [0.0, h]])

    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(8, 4), layout='tight', sharey=True)
    axs[0].set_title('Rough DFN')
    axs[0].set_xlim([-1, 11])
    axs[0].set_ylim([-1, 11])
    axs[0].set_aspect(1)
    axs[0].grid()

    # fractures
    rdfn = RandomDFN(doi=(10.0, 10.0), num=100)
    coord = rdfn.poisson_sample(dist=1.0)
    # pos = dfn.random_sample(doi=(10.0, 10.0), num=100)
    angle = np.random.uniform(30.0, 60.0, rdfn.num)
    scale = np.random.uniform(2.0, 2.0, rdfn.num)
    rdfn.generate_RDFN(coord, angle, scale, amp=(0.0150, 0.005, 0.001))
    rdfn.visualize(axs[0])
    axs[0].add_patch(Polygon(bbox, ec='k', ls='--', lw=1.0, fill=False))

    axs[1].set_title('Smooth DFN')
    axs[1].set_xlim([-1, 11])
    axs[1].set_ylim([-1, 11])
    axs[1].set_aspect(1)
    axs[1].grid()

    sdfn = RandomDFN(num=len(angle), doi=(10.0, 10.0))
    sdfn.generate_SDFN(coord, angle, scale)
    sdfn.visualize(axs[1])
    axs[1].add_patch(Polygon(bbox, ec='k', ls='--', lw=1.0, fill=False))

    plt.savefig('./img/dfn_couple.png', dpi=330)
    plt.show()