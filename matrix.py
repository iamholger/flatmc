
import numpy as np

import math as m
import random as r

from vector import Vec4
from particle import Particle

class eetojj:

    def __init__(self,alphas,ecms=91.2):
        self.alphas = alphas
        self.ecms = ecms
        self.MZ2 = pow(91.1876,2.)
        self.GZ2 = pow(2.4952,2.)
        self.alpha = 1./128.802
        self.sin2tw = 0.22293

    def ME2(self,fl,s,t): #7.51 --- args sind 4 impulse
        qe = -1.
        ae = -0.5
        ve = ae - 2.*qe*self.sin2tw;
        qf = 2./3. if fl in [2,4] else -1./3.
        af = 0.5 if fl in [2,4] else -0.5
        vf = af - 2.*qf*self.sin2tw;
        kappa = 1./(4.*self.sin2tw*(1.-self.sin2tw))
        chi1 = kappa * s * (s-self.MZ2)/(pow(s-self.MZ2,2.) + self.GZ2*self.MZ2);
        chi2 = pow(kappa * s,2.)/(pow(s-self.MZ2,2.) + self.GZ2*self.MZ2);
        term1 = (1+pow(1.+2.*t/s,2.))*(pow(qf*qe,2.)+2.*(qf*qe*vf*ve)*chi1+(ae*ae+ve*ve)*(af*af+vf*vf)*chi2)
        term2 = (1.+2.*t/s)*(4.*qe*qf*ae*af*chi1+8.*ae*ve*af*vf*chi2)
        return pow(4.*m.pi*self.alpha,2.)*3.*(term1+term2)

    def GeneratePoint(self): # hier werden 4impulse gebaut
        ct = 2.*r.random()-1.
        st = m.sqrt(1.-ct*ct)
        phi = 2.*m.pi*r.random()
        p1 = Vec4(1,st*m.cos(phi),st*m.sin(phi),ct)*self.ecms/2
        p2 = Vec4(p1.E,-p1.px,-p1.py,-p1.pz)
        pa = Vec4(self.ecms/2,0,0,self.ecms/2)
        pb = Vec4(self.ecms/2,0,0,-self.ecms/2)
        fl = r.randint(1,5)
        lome = self.ME2(fl,(pa+pb).M2(),(pa-p1).M2())
        dxs = 5.*lome*3.89379656e8/(8.*m.pi)/(2.*pow(self.ecms,2))  # 1./8 ... ist das PS gewicht
        return ( [
            Particle(-11,-pa),
            Particle(11,-pb),
            Particle(fl,p1,[1,0]),
            Particle(-fl,p2,[0,1])
        ], dxs, lome )

    def GenerateLOPoint(self):
        lo = self.GeneratePoint()
        return ( lo[0], lo[1] )

if __name__ == "__main__":
    M = eetojj(0.128)
    import matplotlib.pyplot as plt


    # Spaghetti code for the three particle case
    Q=100 # total momentum
    n=3   # number of external particles

    # Following Simon's notation, in the 3particle case, "i" in Algorithm 1 is 2
    #
    # we solve r= (n+1-i)*u**(n-i) - (n-i)*u**(n-1 + 1)  = 2 * u - u**2 for u
    #
    # we can neglect the second solution and write u_i = u_2 = 1-sqrt(1-r_(i-1))
    #                                                        = 1 - sqrt(1-r_1)
    #
    # where r_1 is our first random number

    def u2(r1):
        """
        Integration variable in the n=3 case
        """
        return 1. -  np.sqrt(1.-r1)

    def rho(Mout, Min, mext=0.0):
        """
        Helper function for mass term eq (5)
        """
        M2 = Mout*Mout
        return 0.125 * np.sqrt( (M2 - (Min+mext)*(Min+mext)) * (M2 - (Min-mext)*(Min-mext))) / M2

    def rho_massless(Mout, Min):
        """
        Helper function for mass term eq (5)
        """
        M2  = Mout*Mout
        M22 = Min*Min
        return 0.125 * np.sqrt( M2*M2 - 2*M2*M22 + M22*M22) / M2

    Q_1 = Q
    M_1 = Q

    from sympy.solvers import solve
    from sympy import Symbol
    u = Symbol("u")

    NSAMPLES=1
    R_1 = np.random.rand(NSAMPLES)

    import time
    # t1 = time.time()
    # U2 = [min(solve(2*u - u**2 -r1, u, positive=True)) for r1 in R1]
    # t2 = time.time()
    # print("Method 1 took {} s".format(t2-t1))

    t1 = time.time()
    U_2 = u2(R_1)
    t2 = time.time()
    print("Method 2 took {} s".format(t2-t1))

    # Now the mass of the intermediary system
    # M_i = u2...u_i*sqrt(Q**2)
    # i.e. M_2 = u_2 * Q
    M_2 = U_2 * Q


    # Now build momenta
    # q_i  = 4*M_(i-1) * rho (M_(i-1), M_i, 0)
    # i.e. = 4 * M_1 * rho(M_1, M_2)

    t1 = time.time()
    q_2 = 4*M_1 * rho_massless(M_1, M_2)
    t2 = time.time()
    print("Method 1 took {} s".format(t2-t1))
    # t1 = time.time()
    # q_22 = 4*M_1 *rho(M_1, M_2)
    # t2 = time.time()
    # print("Method 2 took {} s".format(t2-t1))

    # Get costheta_i in [0,1] and phi_i in  [0,2pi]
    # using random numbers r_(n-5+2*i) and r_(n-4+2*i)
    # i.e. r_2 and r_3
    R_2 = np.random.rand(NSAMPLES)
    R_3 = np.random.rand(NSAMPLES)
    CT_2  = 2.*R_2-1.
    PHI_2 = 2.*np.pi * R_3

    # Three momenta come first
    # p_(i-1) = q_i * [
    #                  cos(phi_i)*sqrt(1-costheta_i**2)
    #                  sin(phi_i)*sqrt(1-costheta_i**2)
    #                  costheta_i
    #                  ]
    COSPHI_2 = np.cos(PHI_2)
    SINPHI_2 = np.sin(PHI_2)
    TMP = np.sqrt(1.-CT_2*CT_2)

    # 4-mom vectors as array
    MOM = np.empty((NSAMPLES,4))
    MOM[:,0] = q_2
    MOM[:,1] = q_2 * COSPHI_2 * TMP
    MOM[:,2] = q_2 * SINPHI_2 * TMP
    MOM[:,3] = q_2 * CT_2

    # Q_2 4-vectors
    Q_2  = np.empty((NSAMPLES, 4))
    Q_2[:,0] = np.sqrt(q_2*q_2 + M_2*M_2)
    Q_2[:,(1,2,3)] = - MOM[:,(1,2,3)]


    t1 = time.time()
    vP1 = [Vec4(*m) for m in MOM] # p1 --- but in CMS
    vQ2 = [Vec4(*q) for q in Q_2] # Q2 --- but in CMS
    t2 = time.time()
    print("Conversion took {} s".format(t2-t1))

    # TODO understand boost --- also if n>3 there is also a boost of Q_2 required to proceed correctly
    #
    # NOTE in the final step, i.e. no more intermediate steps, M_i becomes 0
    # and the final 4-momenta are (in the correponding CMS) simply +/- MOM
    q_3 = 4*M_2 * rho_massless(M_2, 0)

    # in the n=3 case, we are done with intermediate steps and only need to draw two more random numbers
    R_4 = np.random.rand(NSAMPLES)
    R_5 = np.random.rand(NSAMPLES)
    CT_3  = 2.*R_4-1.
    PHI_3 = 2.*np.pi * R_5

    COSPHI_3 = np.cos(PHI_3)
    SINPHI_3 = np.sin(PHI_3)
    TMP = np.sqrt(1.-CT_3*CT_3)

    # 4-mom vectors as array
    MOMLAST = np.empty((NSAMPLES,4))
    MOMLAST[:,0] = q_3
    MOMLAST[:,1] = q_3 * COSPHI_3 * TMP
    MOMLAST[:,2] = q_3 * SINPHI_3 * TMP
    MOMLAST[:,3] = q_3 * CT_3

    MOMLAST_minus = np.empty((NSAMPLES,4))
    MOMLAST_minus[:,0] = q_3
    MOMLAST_minus[:,(1,2,3)] = - MOMLAST[:,(1,2,3)]
    vP2 = [Vec4(*m) for m in MOMLAST] # p2 --- but in CMS
    vP3 = [Vec4(*m) for m in MOMLAST_minus] # p3 --- but in CMS
    from IPython import embed
    embed()


