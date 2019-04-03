#!/usr/bin/env python

import numpy as np

def getEqns(N):
    """
    String representations of equations to solve for N external particles.
    """

    for i in reversed(range(2, N)):
        print("r_{} = {}u^{} - {}u^{}".format(i-1, N-i+1, N-i, N-i, N-i+1))


def mkData(N, npoint=100):
    U=np.linspace(0,1,npoint)
    ret = []
    for i in reversed(range(2, N)):
        ret.append( (N-i+1)*U**(N-i) - (N-i)*U**(N-i+1))
    return U, ret



if __name__ == "__main__":
    import sys

    NP = int(sys.argv[1])
    getEqns(NP)
    import matplotlib.pyplot as plt
    plt.style.use("ggplot")

    # U, YY = mkData(NP)

    # for num, y in enumerate(YY):
        # plt.plot(U, y, label="$r_{}$".format(num+1))
    # plt.legend()
    # plt.show()

    def u2(r1):
        """
        Integration variable in the n=3 case
        """
        return 1. -  np.sqrt(1.-r1)

    def u3(r):
        kappa = np.cbrt(  2*np.sqrt(r*(r-1) ) -2*r +1  )
        return 0.5*(kappa + 1./kappa + 1)

    def f(x, a, r):
        return a*x**(a-1) - (a-1)*x**a - r

    def fp(x, a):
        return a*(a-1)*(x**(a-2) - x**(a-1))

    def fpp(x, a):
        return a*(a-1)*((a-2)*x**(a-3) - (a-1)*x**(a-2))

    a=NP-1



    r = np.random.rand()
    print("Drawn random number r={}".format(r))

    from scipy import optimize
    if a>2:
        res = optimize.newton(lambda x : f(x,a,r), r, fprime=lambda x: fp(x,a), fprime2=lambda x: fpp(x,a))
    else:
        res = optimize.newton(lambda x : f(x,a,r), r, fprime=lambda x: fp(x,a))

    # from IPython import embed
    # embed()

    # X = np.linspace(0,1,100)
    # plt.plot(X, f(X,a,r),   label="f, r={}".format(r))
    # plt.plot(X, fp(X,a),  label="fp")
    # if a>2:
        # plt.plot(X, fpp(X,a), label="fpp")
    # plt.axvline(res, label="u={}".format(res))
    # plt.legend()
    # plt.show()

    # Numerical experimentation --- measure distance between r and u
    NSAMPLES=int(sys.argv[2])#10000

    import time
    U = np.empty(NSAMPLES)
    R = np.empty(NSAMPLES)
    t1=time.time()
    # U = [optimize.newton(lambda x : f(x,a,r), r, fprime=lambda x: fp(x,a), fprime2=lambda x: fpp(x,a)) for r in R]
    # for num, r in enumerate(R):
    for num in range(NSAMPLES):
        r = np.random.rand()
        u = optimize.newton(lambda x : f(x,a,r), r, fprime=lambda x: fp(x,a), fprime2=lambda x: fpp(x,a))
        R[num] = r
        U[num] = u

    t2=time.time()

    print("Evaluation of {} samples with a={} took {} seconds, residual = {}".format(NSAMPLES, a, t2-t1, sum([f(u,a,r) for u, r in zip(U,R)])))


    print("Distance between u and r = {} +/- {}".format(np.mean(R-U), np.std(R-U)))

    plt.scatter(R,U, s=0.01)
    plt.show()

    # plt.hist(U, bins=50, label="u", histtype="step")
    # plt.hist(R, bins=50, label="r", histtype="step")
    # plt.legend()
    # plt.show()

    t1=time.time()
    UNO = [u2(r) for r in R]
    t2=time.time()
    print("Evaluation of {} samples with a={} took {} seconds".format(NSAMPLES, a, t2-t1))
    # from IPython import embed
    # embed()

