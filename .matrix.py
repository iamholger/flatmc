
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
        self.amin = 1.e-10
        self.ye = .5
        self.ze = .01
        self.ws = .25

    def ME2(self,fl,s,t):
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

    def GeneratePoint(self):
        ct = 2.*r.random()-1.
        st = m.sqrt(1.-ct*ct)
        phi = 2.*m.pi*r.random()
        p1 = Vec4(1,st*m.cos(phi),st*m.sin(phi),ct)*self.ecms/2.
        p2 = Vec4(p1.E,-p1.px,-p1.py,-p1.pz)
        pa = Vec4(self.ecms/2,0,0,self.ecms/2.)
        pb = Vec4(self.ecms/2,0,0,-self.ecms/2.)
        fl = r.randint(1,5)
        lome = self.ME2(fl,(pa+pb).M2(),(pa-p1).M2())
        dxs = 5.*lome*3.89379656e8/(8.*m.pi)/(2.*pow(self.ecms,2))
        return ( [
            Particle(-11,-pa),
            Particle(11,-pb),
            Particle(fl,p1,[1,0]),
            Particle(-fl,p2,[0,1])
        ], dxs, lome )

    def GenerateLOPoint(self):
        lo = self.GeneratePoint()
        return ( lo[0], lo[1] )

    def GeneratePOWHEGPoint(self):
        lo = self.GeneratePoint()
        mu = self.ecms
        l = m.log(mu*mu/(lo[0][2].mom+lo[0][3].mom).M2())
        V = [ -2., -3. - 2.*l, -8.+m.pi*m.pi - 3.*l - l*l ]
        I = [ 2., 3. + 2.*l, 19./2.-m.pi*m.pi + 3.*l + l*l ]
        if V[0]+I[0] != 0.:
            print "Pole check failed (double pole)"
        if V[1]+I[1] != 0.:
            print "Pole check failed (single pole)"
        K = self.alphas(mu*mu)/(2.*m.pi)*4./3.*(V[2]+I[2])
        return ( lo[0], lo[1]*(1.+K) )

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

    def Real(self,lome,p1,p2,p3,mu):
        s12 = 2.*p1*p2
        s13 = 2.*p1*p3
        s23 = 2.*p2*p3
        R = lome*(s23/s13+s13/s23+2.*s12*(s12+s13+s23)/(s13*s23))
        cpl = 8.*m.pi*self.alphas(mu*mu)*4./3.
        return cpl*R

    def RSub(self,fl,pa,p1,p2,p3,mu):
        s12 = 2.*p1*p2
        s13 = 2.*p1*p3
        s23 = 2.*p2*p3
        y132 = s13/(s12+s13+s23)
        y231 = s23/(s12+s13+s23)
        if y132 < self.amin or y231 < self.amin:
            return [ ]
        z1 = s12/(s12+s23)
        z2 = s12/(s12+s13)
        p13t = p1+p3-y132/(1.-y132)*p2
        p1t = 1./(1.-y231)*p1
        D132 = 1./y132 * (2./(1.-z1*(1.-y132))-(1.+z1))
        D231 = 1./y231 * (2./(1.-z2*(1.-y231))-(1.+z2))
        D132 *= self.ME2(fl,s12+s13+s23,(pa+p13t).M2())
        D231 *= self.ME2(fl,s12+s13+s23,(pa+p1t).M2())
        cpl = 8.*m.pi*self.alphas(mu*mu)*4./3.
        return [ cpl*D132, cpl*D231 ]

    def GenerateHPoint(self,lo,mu):
        y = pow(r.random(),1./(1.-self.ye));
        w = pow(y,self.ye)/(1.-self.ye);
        z = pow(r.random(),1./(1.-self.ze));
        w *= pow(z,self.ze)/(1.-self.ze);
        phi = 2.*m.pi*r.random()
        w *= (1.-y)/(16.*m.pi*m.pi);
        ij = r.randint(0,1)
        pe = self.MakeKinematics(z,y,phi,lo[0][2+ij].mom,lo[0][3-ij].mom)
        Dijk = self.RSub(lo[0][2].pid,lo[0][0].mom,pe[2*ij],pe[2-2*ij],pe[1],mu)
        if len(Dijk) == 0:
            return ( lo[0], 0. )
        R = self.Real(lo[2],pe[2*ij],pe[2-2*ij],pe[1],mu)
        return ( [
            lo[0][0],lo[0][1],
            Particle(21,pe[1],[2,1]),
            Particle(lo[0][2].pid,pe[2*ij],[1,0]),
            Particle(lo[0][3].pid,pe[2-2*ij],[0,2])
        ], lo[1]/lo[2]*w*(R-Dijk[0]-Dijk[1]) )

    def GenerateSPoint(self,lo,mu):
        l = m.log(mu*mu/(lo[0][2].mom+lo[0][3].mom).M2())
        V = [ -2., -3. - 2.*l, -8.+m.pi*m.pi - 3.*l - l*l ]
        I132 = [ 1., 3./2. + l, 5.-m.pi*m.pi/2. + 3./2.*l + l*l/2. ]
        I231 = [ 1., 3./2. + l, 5.-m.pi*m.pi/2. + 3./2.*l + l*l/2. ]
        if V[0]+I132[0]+I231[0] != 0.:
            print "Pole check failed (double pole)"
        if V[1]+I132[1]+I231[1] != 0.:
            print "Pole check failed (single pole)"
        K = self.alphas(mu*mu)/(2.*m.pi)*4./3.*(V[2]+I132[2]+I231[2])
        return ( lo[0], lo[1]*(1.+K) )

    def GenerateMCNLOPoint(self):
        lo = self.GeneratePoint()
        if r.random() < self.ws:
            nlo = self.GenerateHPoint(lo,self.ecms)
            return ( nlo[0], nlo[1]/self.ws )
        nlo = self.GenerateSPoint(lo,self.ecms)
        return ( nlo[0], nlo[1]/(1.-self.ws) )
