import matplotlib.pyplot as plt
import numpy as np
from scipy import spatial
from skimage import measure
from marching_squares import marching
image = np.load('data.npy')
spacing = 1.0

#xx ,yy = np.meshgrid(np.arange(0,(image.shape[1])*spacing,spacing),np.arange(0,(image.shape[0])*spacing,spacing))
#plt.imshow(image)
#print(ident)

contours = measure.find_contours(image, 0.5,fully_connected = 'high')

mr = marching(image,0.5)
contours = mr.points
for poly in mr.points:
    plt.plot(poly[:,1],poly[:,0],'-')   
#contours = [np.delete(con,-1,axis=0) for con in contours]
points = np.concatenate(contours,axis = 0)
#points = points*spacing
tree = spatial.KDTree(points)

def edge_energy(I,J):
    l = np.array([points[I[0]]-points[J[1]],
                  points[I[1]]-points[J[0]],
                  points[I[0]]-points[J[0]],
                  points[I[1]]-points[J[1]],
                  points[I[0]]-points[I[1]],
                  points[J[0]]-points[J[1]]])
    ld = np.linalg.norm(l,axis=1)
    #stich options
    options = np.array([ld[0]+ld[1],ld[2]+ld[3]])
    E = min(options)-ld[4]-ld[5]
    return E, np.argmin(options)

start = []
end = []
cycle = []
n = 0

for j,con in enumerate(contours):
    for i,_ in enumerate(con):
        cycle.append(j)
        start.append(i+n)
        
        if i < len(con)-1:
            end.append(i+1+n)
        else:
            end.append(n)  
    n = n + len(con)
C = np.array([cycle,start,end]).T



print(np.unique(C[:,0]))
cj = None
corners = None
while len(np.unique(C[:,0]))>1:
    count = [np.sum(C[:,0] == c) for c in np.unique(C[:,0])]
    k = np.unique(C[:,0])[np.argmin(count)]
    print(k)
    index = C[C[:,0] == k,:]
    pts = tree.query_ball_point(points[index[:,1],:],4.0)
    neig = []
    energy = 1e10
    for i,pt in enumerate(pts):
        I = index[i,1:]
        N = [el for el in pt if el not in index[:,1]]
        if len(N)> 0:
            for j in N:
                J = C[C[:,1] == j,1:][0]
                cost,option = edge_energy(I,J)
                if cost < energy:
                    energy = cost
                    corners = [I[0],I[1],J[0],J[1]]
                    cj = C[C[:,1] == j,0][0]
                neig.append(j)
    if cj is None:
        break
    else:
        # stich the cycles !!!
        if option == 0:
            C[C[:,1]== corners[0],2] = corners[3]
            C[C[:,1]== corners[2],2] = corners[1]
        else:
            C[C[:,1]== corners[0],2] = corners[2]
            C[C[:,1]== corners[1],1] = corners[1]            
        C[C[:,0]== cj,0] = k

count = [np.sum(C[:,0] == c) for c in np.unique(C[:,0])]
#print(count)
#print(corners,k,cj)


#plt.imshow(image)

kloop = C[C[:,0] ==k,1:]
for edge in kloop:
    plt.plot(points[edge,1],points[edge,0],':k')





#repulsion stage
if False:
    new_points = points.copy()
    for i,pt in enumerate(points):
        #print(pt)
        neig = tree.query_ball_point(pt,spacing*2.0)
        if len(neig)>0:
            N = len(neig)
            new_points[i,:] = 0.5*(pt+points[neig,:].sum(axis = 0)/N)

    #order the edges to cordinates... 
    kloop_hash = {edge[0]:edge[1] for edge in kloop}
    ids = [kloop[0,0]]
    while len(ids)<len(kloop):
        ids.append(kloop_hash[ids[-1]])

    plt.plot(new_points[ids,1],new_points[ids,0],'--r')
    out = new_points[ids,:]
    print(out)
    np.save('path.npy', out)

plt.show()




