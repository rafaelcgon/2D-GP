
# GPy kernel
# implement a new kernel
#   1) implement the new covariance as a GPy.kern.src.kern.Kern object
#   2) update the GPy.kern.src file


from GPy.kern import Kern
from GPy.core import Param
import numpy as np

class myKernel(Kern):

    def __init__(self,input_dim,active_dim=[0,1],l_df=1.,l_cf=1,ratio=1.):
        super(myKernel, self).__init__(input_dim,active_dim, 'myKern')
        assert input_dim == 2, "For this kernel we assume input_dim=2"
        self.length_df = Param('length_df', l_df)
        self.length_cf = Param('length_cf', l_cf)
        
        self.ratio = Param('ratio', ratio)
        self.link_parameters(self.length_df, self.length_cf, self.ratio)
    def parameters_changed(self):
        # nothing todo here
        pass

    def K(self,X,X2):
        if X2 is None: X2 = X
        p = 2 # number of dimensions   
        dx1 = X[:,0][:,None] - X2[:,0]    
        dx2 = X[:,1][:,None] - X2[:,1]    
        B11 = dx1*dx1
        B12 = dx1*dx2
        B22 = dx2*dx2
        norm = np.sqrt(np.square(dx1)+np.square(dx2))
        # divergence free (df)
        rdf2 = np.square(self.length_df) 
        Cdf = np.square(norm/self.length_df)   
        aux = (p-1) - Cdf
        Adf = np.concatenate([np.concatenate([B11/rdf2+aux,B12/rdf2],axis=1),
                         np.concatenate([B12/rdf2,B22/rdf2+aux],axis=1)],axis=0)
        Cdf = np.concatenate([np.concatenate([Cdf,Cdf],axis=1),
                         np.concatenate([Cdf,Cdf],axis=1)],axis=0)
        Kdf = np.square(1./self.length_df)*np.exp(-Cdf/2.)*Adf 
        # curl free (cf)
        rcf2 = np.square(self.length_cf) 
        Ccf = np.square(norm/self.length_cf)  
        Acf = np.concatenate([np.concatenate([1-B11/rcf2,-B12/rcf2],axis=1),
                         np.concatenate([-B12/rcf2,1-B22/rcf2],axis=1)],axis=0)
        Ccf = np.concatenate([np.concatenate([Ccf,Ccf],axis=1),
                         np.concatenate([Ccf,Ccf],axis=1)],axis=0)
        Kcf = np.square(1./self.length_cf)*np.exp(-Ccf/2.)*Acf 
        return (self.ratio*Kdf)+(1-self.ratio)*Kcf 

    def Kdiag(self,X):
        var = self.ratio*(1/self.length_df**2)+(1-self.ratio)*(1/self.length_cf**2)
        return np.ones(X.shape[0]*X.shape[1])*var

    def update_gradients_full(self, dL_dK, X, X2): # edit this###########3
        if X2 is None: X2 = X
        p = 2 # number of dimensions   
        dx1 = X[:,0][:,None] - X2[:,0]    
        dx2 = X[:,1][:,None] - X2[:,1]    
        B11 = dx1*dx1
        B12 = dx1*dx2
        B22 = dx2*dx2
        norm2 = np.square(dx1)+np.square(dx2)
# derivative of the length scale from divergent-free kernel
        ldf2 = np.square(self.length_df) 
        ldf3 = self.length_df**3
        ldf5 = self.length_df**5
        Cdf = norm2/ldf2   
        aux = (p-1) - Cdf
        Adf = np.concatenate([np.concatenate([B11/ldf2+aux,B12/ldf2],axis=1),
                         np.concatenate([B12/ldf2,B22/ldf2+aux],axis=1)],axis=0)

        dAdf = (2/ldf3)*np.concatenate([np.concatenate([norm2-B11,-B12],axis=1),
                         np.concatenate([-B12,norm2-B22],axis=1)],axis=0)

        dl_df = self.ratio*np.exp(-norm2/(2*ldf2))*(dAdf + Adf*(2*ldf2 - norm2)/(ldf5))  

# derivative of the length scale from curl-free kernel
        lcf2 = np.square(self.length_cf) 
        lcf3 = self.length_cf**3
        lcf5 = self.length_cf**5
        Ccf = norm2/lcf2
        Acf = np.concatenate([np.concatenate([1-B11/lcf2,-B12/lcf2],axis=1),
                         np.concatenate([-B12/lcf2,1-B22/lcf2],axis=1)],axis=0)

        dAcf = (2/lcf3)*np.concatenate([np.concatenate([B11,B12],axis=1),
                         np.concatenate([B12,B22],axis=1)],axis=0)

        dl_cf = (1-self.ratio)*np.exp(-norm2/(2*lcf2))*(dAcf + Acf*(2*lcf2 - norm2)/(lcf5))  

# derivative of the ratio 
        Cdf = np.concatenate([np.concatenate([Cdf,Cdf],axis=1),
                         np.concatenate([Cdf,Cdf],axis=1)],axis=0)
        Kdf = (1./ldf2)*np.exp(-Cdf/2.)*Adf 
        Ccf = np.concatenate([np.concatenate([Ccf,Ccf],axis=1),
                         np.concatenate([Ccf,Ccf],axis=1)],axis=0)
        Kcf = (1./lcf2)*np.exp(-Ccf/2.)*Acf 

        dr = Kdf - Kcf 
        
        self.length_df.gradient = np.sum(dl_df*dL_dK)
        self.length_cf.gradient = np.sum(dl_cf*dL_dK)
        self.ratio.gradient = np.sum(dr*dL_dK)

# 1/x^2exp(-r^2/(2x^2))

    def update_gradients_diag(self, dL_dKdiag, X):
        pass


    def gradients_X(self,dL_dK,X,X2):
        """derivative of the covariance matrix with respect to X."""
        if X2 is None: X2 = X

        p = 2 # number of dimensions   
        dx1 = X[:,0][:,None] - X2[:,0]    
        dx2 = X[:,1][:,None] - X2[:,1]    
        B11 = dx1*dx1
        B12 = dx1*dx2
        B22 = dx2*dx2
        norm2 = np.square(dx1)+np.square(dx2)
        norm = np.sqrt(norm)

    # derivative of divergent-free part 
        ldf2 = np.square(self.length_df) 
        Cdf = norm2/ldf2   
        aux = (p-1) - Cdf
        Adf = np.concatenate([np.concatenate([B11/ldf2+aux,B12/ldf2],axis=1),
                         np.concatenate([B12/ldf2,B22/ldf2+aux],axis=1)],axis=0)
        dAdf = (2/ldf2)*np.concatenate([np.concatenate([dx1+norm,(dx1+dx2)/2.],axis=1),
                         np.concatenate([(dx1+dx2)/2.,dx2+norm],axis=1)],axis=0)
        dX_df = self.ratio*np.exp(-norm2/(2*ldf2))*(dAdf - Adf*norm/ldf2)/ldf2  
    # derivative of curl-free part
        lcf2 = np.square(self.length_cf) 
        Ccf = norm2/lcf2
        Acf = np.concatenate([np.concatenate([1-B11/lcf2,-B12/lcf2],axis=1),
                         np.concatenate([-B12/lcf2,1-B22/lcf2],axis=1)],axis=0)
        dAcf = (1/lcf2)*np.concatenate([np.concatenate([-2*dx1,-dx1-dx2],axis=1),
                         np.concatenate([-dx1-dx2,-2*dx2],axis=1)],axis=0)
        dX_cf = (1-self.ratio)*np.exp(-norm2/(2*lcf2))*(dAcf - Acf*norm/lcf2)/lcf2  
        return np.sum(dL_dK*(dX_df+dX_cf),1)[:,None]

    def gradients_X_diag(self,dL_dKdiag,X):
        # no diagonal gradients
        pass

