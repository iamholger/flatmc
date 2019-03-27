import math as m
import random as r

from vector import Vec4
from particle import Particle, CheckEvent
from qcd import AlphaS, NC, TR, CA, CF

class Kernel:

    def __init__(self,flavs):
        self.flavs = flavs

class Pqq (Kernel):

    def Value(self,z,y):
        return CF*(2./(1.-z*(1.-y))-(1.+z))

    def Estimate(self,z):
        return CF*2./(1.-z)

    def Integral(self,zm,zp):
        return CF*2.*m.log((1.-zm)/(1.-zp))

    def GenerateZ(self,zm,zp):
        return 1.+(zp-1.)*m.pow((1.-zm)/(1.-zp),r.random())

class Pgg (Kernel):

    def Value(self,z,y):
        return CA/2.*(2./(1.-z*(1.-y))-2.+z*(1.-z))

    def Estimate(self,z):
        return CA/(1.-z)

    def Integral(self,zm,zp):
        return CA*m.log((1.-zm)/(1.-zp))

    def GenerateZ(self,zm,zp):
        return 1.+(zp-1.)*m.pow((1.-zm)/(1.-zp),r.random())

class Pgq (Kernel):

    def Value(self,z,y):
        return TR/2.*(1.-2.*z*(1.-z))

    def Estimate(self,z):
        return TR/2.

    def Integral(self,zm,zp):
        return TR/2.*(zp-zm)

    def GenerateZ(self,zm,zp):
        return zm+(zp-zm)*r.random()

class Shower:

    def __init__(self,alpha,t0):
        self.t0 = t0
        self.alpha = alpha
        self.alphamax = alpha(self.t0)
        self.kernels = [ Pqq([fl,fl,21]) for fl in [-5,-4,-3,-2,-1,1,2,3,4,5] ]
        self.kernels += [ Pgq([21,fl,-fl]) for fl in [1,2,3,4,5] ]
        self.kernels += [ Pgg([21,21,21]) ]

    def MakeKinematics(self,z,y,phi,pijt,pkt):
        Q = pijt+pkt
        rkt = m.sqrt(Q.M2()*y*z*(1.-z))
        kt1 = pijt.Cross(pkt)
        if kt1.P() < 1.e-6: kt1 = pijt.Cross(Vec4(0.,1.,0.,0.))
        kt1 *= rkt*m.cos(phi)/kt1.P()
        kt2cms = Q.Boost(pijt).Cross(kt1)
        kt2cms *= rkt*m.sin(phi)/kt2cms.P()
        kt2 = Q.BoostBack(kt2cms)
        pi = z*pijt + (1.-z)*y*pkt + kt1 + kt2
        pj = (1.-z)*pijt + z*y*pkt - kt1 - kt2
        pk = (1.-y)*pkt
        return [pi,pj,pk]

    def MakeColors(self,flavs,colij,colk):
        self.c += 1
        if flavs[0] != 21:
            if flavs[0] > 0:
                return [ [self.c,0], [colij[0],self.c] ]
            else:
                return [ [0,self.c], [self.c,colij[1]] ]
        else:
            if flavs[1] == 21:
                if colij[0] == colk[1]:
                    if colij[1] == colk[0] and r.random()>0.5:
                        return [ [colij[0],self.c], [self.c,colij[1]] ]
                    return [ [self.c,colij[1]], [colij[0],self.c] ]
                else:
                    return [ [colij[0],self.c], [self.c,colij[1]] ]
            else:
                if flavs[1] > 0:
                    return [ [colij[0],0], [0,colij[1]] ]
                else:
                    return [ [0,colij[1]], [colij[0],0] ]

    def GeneratePoint(self,event):
        while self.t > self.t0:
            t = self.t0
            for split in event[2:]:
                for spect in event[2:]:
                    if spect == split: continue
                    if not split.ColorConnected(spect): continue
                    for sf in self.kernels:
                        if sf.flavs[0] != split.pid: continue
                        m2 = (split.mom+spect.mom).M2()
                        if m2 < 4.*self.t0: continue
                        zp = .5*(1.+m.sqrt(1.-4.*self.t0/m2))
                        g = self.alphamax/(2.*m.pi)*sf.Integral(1.-zp,zp)
                        tt = self.t*m.pow(r.random(),1./g)
                        if tt > t:
                            t = tt
                            s = [ split, spect, sf, m2, zp ]
            self.t = t
            if t > self.t0:
                z = s[2].GenerateZ(1.-s[4],s[4])
                y = t/s[3]/z/(1.-z)
                if y < 1.:
                    f = (1.-y)*self.alpha(t)*s[2].Value(z,y)
                    g = self.alphamax*s[2].Estimate(z)
                    if f/g > r.random():
                        phi = 2.*m.pi*r.random()
                        moms = self.MakeKinematics(z,y,phi,s[0].mom,s[1].mom)
                        cols = self.MakeColors(s[2].flavs,s[0].col,s[1].col)
                        event.append(Particle(s[2].flavs[2],moms[1],cols[1]))
                        s[0].Set(s[2].flavs[1],moms[0],cols[0])
                        s[1].mom = moms[2]
                        return
    
    def Run(self,event,t):
        self.c = 1
        self.t = t
        while self.t > self.t0:
            self.GeneratePoint(event)
            
# build and run the generator

import sys
from matrix import eetojj
from durham import Analysis

alphas = AlphaS(91.1876,0.118)
hardxs = eetojj(alphas)
shower = Shower(alphas,t0=1.)
jetrat = Analysis()

r.seed(123456)
for i in range(100000):
    event, weight = hardxs.GenerateLOPoint()
    t = (event[0].mom+event[1].mom).M2()
    shower.Run(event,t)
    sys.stdout.write('\rEvent {0}'.format(i))
    sys.stdout.flush()
    jetrat.Analyze(event,weight)
jetrat.Finalize("myshower")
print ""
