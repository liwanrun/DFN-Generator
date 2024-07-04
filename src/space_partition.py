import numpy as np

class CellPartition:
    def __init__(self, doi, spa) -> None:
        nrows = int(np.ceil(doi[0] / spa))
        ncols = int(np.ceil(doi[1] / spa))
        self.doi = doi
        self.spa = spa
        self.shape = (nrows, ncols)
        self.Points = [[[] for _ in range(ncols)] for _ in range(nrows)]
        self.Curves = [[[] for _ in range(ncols)] for _ in range(nrows)]

    def assign_point_to_cell(self, point):
        cells = self.get_incident_cells(point)
        for bin in cells:
            i, j = bin
            self.Points[i][j].append(point)

    def assign_curve_to_cells(self, curve):
        cells = self.get_incident_cells(curve)
        for bin in cells:
            i, j = bin
            self.Curves[i][j].append(curve)

    def get_incident_cells(self, geometry):
        nrows, ncols = self.shape
        minx, miny, maxx, maxy = geometry.bounds
        r_beg = max(int(np.floor(miny / self.spa)), 0)
        r_end = min(int(np.ceil(maxy / self.spa)), nrows)
        c_beg = max(int(np.floor(minx / self.spa)), 0)
        c_end = min(int(np.ceil(maxx / self.spa)), ncols)
        I = np.dstack(np.mgrid[r_beg:r_end, c_beg:c_end])
        I = I.reshape(I.size//2, 2).tolist()
        return I

# Unit test
if __name__ == '__main__':
    
    from shapely import Point, LineString, Polygon
    from shapely.plotting import plot_line, plot_polygon
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator

    searcher = CellPartition((10.0, 10.0), 1.0)

    point = Point(np.array([5.0, 5.0]))
    P = searcher.get_incident_cells(point)
    print(P)

    line = LineString(np.array([[0.5, 0.5], [4.5, 5.1]]))
    L = searcher.get_incident_cells(line)
    print(L)
    
    fig, axs = plt.subplots()
    axs.xaxis.set_major_locator(MaxNLocator(10))
    axs.yaxis.set_major_locator(MaxNLocator(10))
    axs.set_xlim([0.0, 10.0])
    axs.set_ylim([0.0, 10.0])
    axs.grid()
    axs.set_aspect(1)
    plot_line(line, axs)
    plot_polygon(line.envelope)
    plt.show()


    