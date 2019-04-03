
import math as m

from vector import Vec4
from particle import Particle

class Algorithm:

    def Yij(self,p,q):
        pq = p.px*q.px+p.py*q.py+p.pz*q.pz
        return 2.*pow(min(p.E,q.E),2)*(1.0-min(max(pq/m.sqrt(p.P2()*q.P2()),-1.0),1.0))/self.ecm2

    def Cluster(self,event):
        self.ecm2 = (event[0].mom+event[1].mom).M2()
        p = [ i.mom for i in event[2:] ]
        kt2 = []
        n = len(p)
        imap = range(n)
        kt2ij = [ [ 0 for i in range(n) ] for j in range(n) ]
        dmin = 1
        for i in range(n):
            for j in range(i):
                dij = kt2ij[i][j] = self.Yij(p[i],p[j])
                if dij < dmin: dmin = dij; ii = i; jj = j
        while n>2:
            n -= 1
            kt2.append(dmin);
            jjx = imap[jj]
            p[jjx] += p[imap[ii]]
            for i in range(ii,n): imap[i] = imap[i+1]
            for j in range(jj): kt2ij[jjx][imap[j]] = self.Yij(p[jjx],p[imap[j]])
            for i in range(jj+1,n): kt2ij[imap[i]][jjx] = self.Yij(p[jjx],p[imap[i]])
            dmin = 1
            for i in range(n):
                for j in range(i):
                    dij = kt2ij[imap[i]][imap[j]]
                    if dij < dmin: dmin = dij; ii = i; jj = j
        return kt2

from histogram import Histo1D

class Analysis:

    def __init__(self,n=1):
        self.n = 0.
        self.ynm = [[ Histo1D(100,-4.3,-0.3,'/LL_JetRates/log10_y_{0}{1}\n'.format(i+2,i+3))
                      for i in range(4) ] for i in range(n+1) ]
        self.duralg = Algorithm()

    def Analyze(self,event,w,weights=[]):
        self.n += 1.
        kt2 = self.duralg.Cluster(event)
        for j in range(len(self.ynm[0])):
            self.ynm[0][j].Fill(m.log10(kt2[-1-j]) if len(kt2)>j else -5.,w)
            for i in range(len(weights)):
                self.ynm[i+1][j].Fill(m.log10(kt2[-1-j]) if len(kt2)>j else -5.,w*weights[i])

    def Finalize(self,name):
        for h in self.ynm[0]: h.ScaleW(1./self.n)
        file = open(name+".yoda","w")
        file.write("\n\n".join([ str(h) for h in self.ynm[0] ]))
        file.close() 
        for i in range(1,len(self.ynm)):
            for h in self.ynm[i]: h.ScaleW(1./self.n)
            file = open(name+".{0}.yoda".format(i),"w")
            file.write("\n\n".join([ str(h) for h in self.ynm[i] ]))
            file.close() 
