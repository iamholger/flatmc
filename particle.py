
import math as m

from vector import Vec4

class Particle:

    def __init__(self,pdgid,momentum,col=[0,0]):
        self.Set(pdgid,momentum,col)

    def __repr__(self):
        return "{0} {1} {2}".format(self.pid,self.mom,self.col)

    def __str__(self):
        return "{0} {1} {2}".format(self.pid,self.mom,self.col)

    def Set(self,pdgid,momentum,col=[0,0]):
        self.pid = pdgid
        self.mom = momentum
        self.col = col

    def ColorConnected(self,p):
        return (self.col[0] > 0 and self.col[0] == p.col[1]) or \
               (self.col[1] > 0 and self.col[1] == p.col[0])

def CheckEvent(event):
    psum = Vec4()
    csum = {}
    for p in event:
        psum += p.mom
        if p.col[0] > 0: 
            csum[p.col[0]] = csum.get(p.col[0],0) + 1
            if csum[p.col[0]] == 0: del csum[p.col[0]]
        if p.col[1] > 0:
            csum[p.col[1]] = csum.get(p.col[1],0) - 1
            if csum[p.col[1]] == 0: del csum[p.col[1]]
    return (m.fabs(psum.E)<1.e-12 and \
            m.fabs(psum.px)<1.e-12 and \
            m.fabs(psum.py)<1.e-12 and \
            m.fabs(psum.pz)<1.e-12 and \
            len(csum) == 0)
