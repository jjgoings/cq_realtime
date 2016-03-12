from __future__ import division
import numpy as np
import sys
from sklearn.linear_model import OrthogonalMatchingPursuit
from sklearn.linear_model import OrthogonalMatchingPursuitCV

def CSSK(h,const=5.0,noise=0.0000001):
    """Compressed Sensing replacement of Fourier Transform on 1D array h
       * REQUIRES CVXPY PACKAGE *
         h       = sampled time signal
         const   = scalar multiple dimension of h, larger values give greater
                     resolution albeit with increased cost.
         noise   = scalar constant to account for numerical noise

         returns:
         g       = fourier transform h to frequency domain using CS technique
    """

    h = np.asarray(h, dtype=float)
    Nt = len(h)
    Nw = int(const*Nt)
    t = np.arange(Nt)
    w = np.arange(Nw)
    #F = np.sin(2 * np.pi * np.outer(t,w) / Nw)
    F = (1/np.float(Nw))*np.sin(2.0*np.pi*np.outer(t,w)/np.float(Nw))

    #omp_cv = OrthogonalMatchingPursuit(n_nonzero_coefs=n_nonzero_coefs)
    #omp_cv = OrthogonalMatchingPursuitCV(verbose=True,normalize=True)
    omp_cv = OrthogonalMatchingPursuit(tol=noise)
    omp_cv.fit(F, h)
    coef = omp_cv.coef_
    #idx_r, = coef.nonzero()
    g = coef


    ### begin using cvxpy
    #g = cvx.Variable(Nw)
    ## min |g|_1 subject to |F.g - h|_2 < noise
    #objective = cvx.Minimize(cvx.norm(g,1))
    #constraints = [cvx.norm(F*g - h,2) <= noise]
    #prob = cvx.Problem(objective, constraints)
    #prob.solve(solver='SCS',verbose=True)
    #g = np.asarray(g.value)
    #g = g[:,0]
    ### end using cvxpy
    return g

if __name__ == '__main__':

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from scipy.fftpack import fftfreq, fft


    t = np.arange(0,200,1.0) 
    h = (np.sin(2.0*t) + 2.0*np.sin(t) + 0.5*np.sin(1.5*t))

    g = CSSK(h,10.0)
    #g = np.imag(fft(h))

    w = fftfreq(len(g),d=(t[1]-t[0]))*2.0*np.pi

    plt.plot(w,abs(g))
    plt.xlim(0,2.5)
    plt.savefig('cs_sklearn.pdf')


    
    



