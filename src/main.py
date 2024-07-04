import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.ticker import MaxNLocator
import shapely
import shapely.plotting
from fracture import RoughFracture
from space_partition import CellPartition
from random_dfn import RandomDFN

doi = (10.0, 10.0)      # Problem domain
spacing = 0.5
bkgmesh = CellPartition(doi, spacing/np.sqrt(2))
dfn = RandomDFN(num=100, doi=doi, msh=bkgmesh)
# Set_1
dfn.guess_times = 100
dfn.minimum_gap = 0.5
dfn.generate_FBYF(dip=(30., 30.), scf=(1.0, 2.0))
print(f'[ INFO ]: generate DFN of {len(dfn.fractures)}.\n')
# Set_2
dfn.guess_times = 500 
dfn.minimum_gap = 0.5
dfn.generate_FBYF(dip=(120., 120.), scf=(1.0, 2.0))
print(f'[ INFO ]: generate DFN of {len(dfn.fractures)}.\n')

fig, axs = plt.subplots(figsize=(6, 6))
axs.set_xlim([-1.0, 11.0])
axs.set_ylim([-1.0, 11.0])
axs.set_aspect(1)
width = spacing / np.sqrt(2)
nrows = int(np.ceil(doi[0] / width))
ncols = int(np.ceil(doi[1] / width))
axs.xaxis.set_major_locator(MaxNLocator(ncols))
axs.yaxis.set_major_locator(MaxNLocator(nrows))
# axs.grid()

w = doi[0]; h = doi[1]
bbox = np.array([[0.0, 0.0], [w, 0], [w, h], [0.0, h]])
axs.add_patch(Polygon(bbox, ec='k', ls='--', lw=1.0, fill=False))
dfn.visualize(axs)

# for poly in dfn.frac_bbox:
#     shapely.plotting.plot_polygon(poly, axs, add_points=False, fill=False)
plt.show()

# Output
dfn.write_to_gmsh_model('dfn.geo')
print('[ INFO ]: Output DFN for meshing.\n')


