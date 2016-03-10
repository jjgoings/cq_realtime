from __future__ import division
import numpy as np
import sys

def CS(h,const=5.0,noise=0.0000001):
    """Compressed Sensing replacement of Fourier Transform on 1D array h
       * REQUIRES CVXPY PACKAGE *
         h       = sampled time signal
         const   = scalar multiple dimension of h, larger values give greater
                     resolution albeit with increased cost.
         noise   = scalar constant to account for numerical noise

         returns:
         g       = fourier transform h to frequency domain using CS technique
    """

    try:
        import cvxpy as cvx
    except ImportError as e:
        print "You need CVXPY: http://www.cvxpy.org/"
        sys.exit()

    h = np.asarray(h, dtype=float)
    Nt = len(h)
    Nw = int(const*Nt)
    t = np.arange(Nt)
    w = np.arange(Nw)
    F = np.sin(2 * np.pi * np.outer(t,w) / Nw)
    g = cvx.Variable(Nw)
    objective = cvx.Minimize(cvx.norm(g,1))
    constraints = [cvx.norm(F*g - h,2) <= noise]
    prob = cvx.Problem(objective, constraints)
    prob.solve(solver='SCS',verbose=True)
    g = np.asarray(g.value)
    g = g[:,0]
    return g

