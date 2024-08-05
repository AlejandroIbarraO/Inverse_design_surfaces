import igl
import numpy as np
import scipy as sp

class local_global:
    def __init__(self,filename,flavor =None):
        v, f = igl.read_triangle_mesh(filename)
        self.v = v
        self.f = f
        self.u = v[:,:2]
        self.hist_u = [self.u]
        self.memory_lenght = 10
        self.precomputation_state = False
        self.flavor = flavor
        # by default ARAP
        self.lam1 = 1.0
        self.lam2 = 1.0
    def precomputation(self):
        self.grad = igl.grad(self.v,self.f)
        cot = igl.cotmatrix_entries(self.v,self.f)
        cot = np.roll(cot,1,axis = 1) # roll the index, so we get the edges in the order [0,1] [1,2] [2,0]
        # given a id of the vertices of an edge retrun the indice of the triangle
        n = self.v.shape[0]
        N = igl.adjacency_list(self.f)
        ij = -1*np.ones((n,n),dtype = np.int32)
        cot_ij = -1*np.ones((n,n))
        K = np.zeros((n,n))
        for tri_ind,tri in enumerate(self.f):
            for i in range(3):
                ij[tri[i],tri[(i+1)%3]] = tri_ind
                cot_ij[tri[i],tri[(i+1)%3]] = 2.0*cot[tri_ind,i]
        for i in range(n):
            for j in N[i]:
                t = ij[i,j]
                t_ji = ij[j,i]
                if t >-1:
                    K[i,j] -=cot_ij[i,j]
                    K[i,i] +=cot_ij[i,j]
                if t_ji >-1:
                    K[i,j] -=cot_ij[j,i]
                    K[i,i] +=cot_ij[j,i]    
        Kext = np.zeros((2*n,2*n))
        Kext[:n,:n] =K
        Kext[n:,n:] =K
        self.precomputation_state = True
        self.n = n 
        self.N = N
        self.cot = cot
        self.ij = ij
        self.cot_ij = cot_ij
        self.Kext = Kext

    def constrain(self,constrains):
        self.constraints = constrains
        if self.precomputation_state:
            if constrains.size>0:
                bn = np.append(constrains,constrains+self.n)
                self.bn = bn
                self.Kext[bn,:] = 0
                self.Kext[bn,bn] = 1 
                self.ufix = self.u.flatten('F')[bn]
            return True
        else:
            print('Make the precomputation step first')
            return False

    def covariance_matrix(self,u):
        J = self.grad.dot(u)
        J = np.reshape(J,(len(self.f),3,2), order = 'F')
        J = np.swapaxes(J, 1,2)    
        return J

    def cal_b(self,L):
        B = np.zeros((self.n,2))
        for i in range(self.n):
            for j in self.N[i]:
                t = self.ij[i,j]
                t_ji = self.ij[j,i]
                if t >-1:
                    B[i] +=self.cot_ij[i,j]*np.dot(L[t],self.v[i]-self.v[j])
                if t_ji >-1:
                    B[i] +=self.cot_ij[j,i]*np.dot(L[t_ji],self.v[i]-self.v[j])
        B = B.flatten('F')
        if self.constraints.size>0:
            B[self.bn] = self.ufix
        return B
    def local_min(self):
        J = self.covariance_matrix(self.u)
        L = np.zeros((self.f.shape[0],2,3))
        if self.flavor == None:
            lam1,lam2 = [self.lam1,self.lam2]
        elif self.flavor == 'Fabric':
            lam1,lam2 = [np.pi/2,1.0]
        elif self.flavor == 'Baromophs':
            lam1 = 1.0
        for i,_ in enumerate(self.f):
            U,c_sig,V = np.linalg.svd(J[i])
            if self.flavor == 'Baromophs':
                lam2 = (c_sig[1]<0.9)*c_sig[1]+(c_sig[1]>=0.9)*0.9
            sigm = np.array([[lam1,0,0],[0,lam2,0]])
            L[i,:,:] = np.dot(U,np.dot(sigm,V))
        return L,J
    def global_min(self,L, save = False):
        b = self.cal_b(L)
        u = np.linalg.solve(self.Kext, b).reshape((-1,2),order = 'F')
        self.u = u
        if save == True:
            self.hist_u.append(u)
            if len(self.hist_u)>self.memory_lenght:
                self.hist_u.pop(0)
    def local_global_step(self):
        if self.precomputation_state == True:
            L,_ = self.local_min()
            self.global_min(L)
    def local_global_n_step(self,n_steps):
        if self.precomputation_state == True:
            for _ in range(n_steps):
                L,_ = self.local_min()
                self.global_min(L)
                
    def strech_and_angle(self):
        J = self.covariance_matrix(self.u)
        lms = np.zeros((self.f.shape[0],2))
        angle = np.zeros((self.f.shape[0],))
        for i,_ in enumerate(self.f):
            U,sig,_ = np.linalg.svd(J[i])
            lms[i]=sig
            angle[i] = np.arctan2(-U[0][1],U[0][0])
        centers = igl.barycenter(self.u,self.f)
        return centers,lms,angle
    def angles(self):
        J = self.covariance_matrix(self.u)
        angle = np.zeros((self.f.shape[0],))
        for i,_ in enumerate(self.f):
            U,_,_ = np.linalg.svd(J[i])
            angle[i] = np.arctan2(-U[0][1],U[0][0])
        centers = igl.barycenter(self.u,self.f)
        return centers,angle
    def scalar_to_vertex(self,scalar):
        scalar_v = np.zeros((len(self.v),))
        count = np.zeros((len(self.v),))
        for i,tri in enumerate(self.f):
            for el in tri:
                scalar_v[el] += scalar[i]
                count[el] += 1.0
        scalar_v = scalar_v/count
        return scalar_v    
    def angle_to_vertice(self):
        _,angles = self.angles() 
        fx = np.cos(angles)
        fy = np.sin(angles)
        vx = self.scalar_to_vertex(fx)
        vy = self.scalar_to_vertex(fy)
        centers = igl.barycenter(self.u,self.f)
        return centers,np.arctan2(vy,vx)