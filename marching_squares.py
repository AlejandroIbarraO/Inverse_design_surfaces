import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry.polygon import Polygon
from shapely import is_ccw

def extractPolygons(edges):
    """
        Extract connected polygonal chains from set of edges.
        The polygonal chains are grown in both directions from one starting edge,
        appending och prepending edges as they are encountered while traversing
        the list of unconnected edges.
        The head and tail vertex indices are updated as the chain is grown until
        no more edges to connect are left.
        A completed polygonal chain is stored in a list of such polygons.
        copy from : https://nils-olovsson.se/articles/marching_squares/
    """
    polygons = []
    while edges:
        chain = []
        e = edges.pop()
        chain.append(e)
        tail = e[0]
        head = e[1]
        cont = True
        while cont and edges:
            e = None
            index = -1
            for i, ne in enumerate(edges):
                t = ne[0]
                h = ne[1]
                if head == t:
                    e = ne
                    chain.append(e)
                    head = h
                    index = i
                    break
                elif tail == h:
                    e = ne
                    chain.insert(0, e)
                    tail = t
                    index = i
                    break
                elif head == h:
                    e = ne
                    chain.append(e)
                    head = t
                    index = i
                    break
                elif tail == t:
                    e = ne
                    chain.insert(0, e)
                    tail = h
                    index = i
                    break
            if e is not None:
                edges.pop(index)
            else:
                cont = False
        polygons.append(chain)
    return polygons


class marching:
    def __init__(self,data,iso):
        self.data = data
        self.mask = data>iso
        self.iso = iso
        self.edges = []
        self.vertices = []
        self.polygons = []
        j,i = np.meshgrid(np.arange(0,data.shape[1]),np.arange(0,data.shape[0]))
        self.index = np.array([j.flatten(),i.flatten()]).T
        self.flat = np.arange(0,data.shape[0]*data.shape[1])
        self.flat2ij = np.array([[0,0],[1,0],[1,1],[0,1]])
        self.lookup = {0:[0,1],1:[1,2],2:[2,3],3:[0,3]}
        self.edge_k = {1:[0,3],
                       2:[0,1],
                       3:[3,1],
                       4:[1,2],
                       6:[0,2],
                       7:[3,2],
                       8:[2,3],
                       9:[2,0],
                       11:[2,1],
                       12:[1,3],
                       13:[1,0],
                       14:[0,3]}
        self.points = []


        self.march()

        edges = [[tail,head] for tail,head in self.vertices]
        self.polygons = extractPolygons(edges)
        self.polycolection = []

        for poly in self.polygons:
            self.polycolection.append([list(el[0]) for el in poly])
            self.polycolection[-1].append(list(poly[-1][1]))
        self.edgetopoint()

    def edgetopoint(self):
        for polyline in self.polycolection:
            vert = np.array(polyline)
            refp0 = vert[:,0]
            refp1 = vert[:,1]
            v0 = self.data.flat[refp0]
            v1 = self.data.flat[refp1]
            t = (self.iso-v0)/(v1-v0)
            p0 = np.unravel_index(refp0,self.data.shape)
            p0 = np.array([p0[0],p0[1]])

            p1 = np.unravel_index(refp1,self.data.shape)
            p1 = np.array([p1[0],p1[1]])
            pint = p0*(1-t)+t*p1
            pint = pint.T
            pint = np.flipud(pint)
            poly = Polygon(pint.tolist())
            print(is_ccw(poly))
            self.points.append(pint)
        
    def edge2vertice(self,edge,ij):
        self.lookup = {0:[0,1],1:[1,2],2:[2,3],3:[0,3]}
        v1 = self.flat2ij[self.lookup[edge[0]]]+ij
        v2 = self.flat2ij[self.lookup[edge[1]]]+ij
        v1f = np.sort(self.ij2flat(v1[:,1],v1[:,0]))
        v2f = np.sort(self.ij2flat(v2[:,1],v2[:,0]))
        self.vertices.append([(v1f[0],v1f[1]),(v2f[0],v2f[1])])

    def ij2flat(self,i,j):
        flat = self.data.shape[1]*j+i
        return flat
    def lut(self,coord):
        j,i = coord
        cord = np.array(coord)
        k = 0
        for l,ij in zip([1,2,4,8],self.flat2ij):
            k += l*self.mask[j+ij[0],i+ij[1]]
        if self.edge_k.get(int(k)) is not None:
            e = self.edge_k[k]
            self.edge2vertice(e,cord)
        elif k == 5:
            for ki in [7,13]:
                e = self.edge_k[ki]
                self.edge2vertice(e,cord)
        elif k == 10:
            for ki in [11,14]:
                e = self.edge_k[ki]
                self.edge2vertice(e,cord)

    def march(self):
        # march over the matrix and save the crossed levels in the vertices array
        m,n = self.data.shape
        for j in range(m-1):
            for i in range(n-1):
                self.lut([j,i])



        
def main():
    d = 100
    x,y = np.meshgrid(np.linspace(0,d,d),np.linspace(0,d,d))
    data = np.exp(-((x-0)**2+(y-0)**2)/(d))
    for i in range(50):
        x0,y0 = np.random.rand(2)
        s = np.random.rand(1)
        data = data + np.exp(-((x-x0*d)**2+(y-y0*d)**2)/(d*s))
    mr = marching(data,0.9)
    #print(mr.vertices)
    plt.imshow(mr.mask)
    #print(mr.mask)
    for poly in mr.points:
        plt.plot(poly[:,1],poly[:,0],'-k')   
    #points = np.array(mr.points)
    #plt.plot(points[:,0],points[:,1],'.k')
    plt.show()
    return mr
   
if __name__ == '__main__':
    mr = main()