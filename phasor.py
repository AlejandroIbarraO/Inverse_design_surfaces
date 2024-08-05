import numpy as np
from scipy import spatial

class phasor:
    def __init__(self):

        self.n_sync = 20
        self.phase_sync = True
    def phasor_compute(self,sources,angles,frec = None,r_neig =5.0):
        # compute the phasors parameters given the sources positions
        # setup the frequencies 
        self.p_phasors = sources
        self.phasor_quad = spatial.KDTree(sources)
        
        if frec is None:
            self.frec = np.ones((len(sources),))
        elif isinstance(frec, float):
            self.frec = np.ones((len(sources),))*frec
        elif isinstance(frec, list):
            self.frec = np.array(frec)
        else:
            self.frec = frec
        self.phasor_phase = np.zeros_like(angles)
        # compute the directions 
        self.source_uv = np.array([np.cos(angles),np.sin(angles)]).T
        # syc phases.
        if self.phase_sync:
            #a = np.zeros_like(frec)
            N = self.phasor_quad.query_ball_point(sources,r_neig)
            
            new_phase = np.zeros_like(self.phasor_phase,dtype=np.complex_)
            for _ in range(self.n_sync):
                for j,source in enumerate(sources):
                    i = N[j]
                    #print(i)
                    p = self.p_phasors[j]-self.p_phasors[i]
                    a = np.fmax(0*self.frec[i],np.dot(self.source_uv[i],self.source_uv[j]))*np.exp(2*np.pi*1j*self.frec[i]*np.sum(p*self.source_uv[i],axis = 1)+1j*self.phasor_phase[i])
                    new_phase[j] = np.angle(a.sum())
                self.phasor_phase = new_phase
    def eval(self,positions,range):
        gabor = np.zeros((len(positions),),dtype=np.complex_)
        for j,pos in enumerate(positions):
            near_phasors = self.phasor_quad.query_ball_point(pos,range)
            p = pos-self.p_phasors[near_phasors]
            l = np.linalg.norm(p,1)
            A = np.exp(-(l*self.frec[near_phasors]/3)**2)
            phi = 2*np.pi*np.sum(self.source_uv[near_phasors]*p,axis=1)*self.frec[near_phasors]+self.phasor_phase[near_phasors]
            #print(phi.shape,self.phasor_phase[near_phasors].shape)
            gabor[j] = np.sum(A*np.exp(phi*1j))
        phasor_noise = np.angle(gabor)
        return phasor_noise