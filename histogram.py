
import sys

class Bin1D:

    def __init__(self,xmin,xmax):
        self.xmin = xmin
        self.xmax = xmax
        self.w = 0.
        self.w2 = 0.
        self.wx = 0.
        self.wx2 = 0.
        self.n = 0.

    def __repr__(self):
        return str(self)

    def __str__(self):
        return "{0:10.6e}\t{1:10.6e}\t{2:10.6e}\t{3:10.6e}\t{4:10.6e}\t{5:10.6e}\t{6}".format\
            (self.xmin,self.xmax,self.w,self.w2,self.wx,self.wx2,int(self.n))

    def Format(self,tag):
        return "{0}\t{0}\t{1:10.6e}\t{2:10.6e}\t{3:10.6e}\t{4:10.6e}\t{5}".format\
            (tag,self.w,self.w2,self.wx,self.wx2,int(self.n))

    def Fill(self,x,w):
        self.w += w
        self.w2 += w*w
        self.wx += w*x
        self.wx2 += w*w*x
        self.n += 1.

    def ScaleW(self,scale):
        self.w *= scale
        self.w2 *= scale*scale
        self.wx *= scale
        self.wx2 *= scale*scale

class Histo1D:

    def __init__(self,nbin,xmin,xmax,name):
        self.name = name
        self.bins = []
        self.uflow = Bin1D(-sys.float_info.max,xmin)
        width = (xmax-xmin)/nbin
        for i in range(nbin):
            self.bins.append(Bin1D(xmin+i*width,xmin+(i+1)*width))
        self.oflow = Bin1D(xmax,sys.float_info.max)
        self.total = Bin1D(-sys.float_info.max,sys.float_info.max)
        self.scale = 1.

    def __repr__(self):
        return str(self)

    def __str__(self):
        s = "BEGIN YODA_HISTO1D {0}\n\n".format(self.name)
        s += "Path={0}\n\n".format(self.name)
        s += "ScaledBy={0}\n".format(self.scale)
        s += "Title=\nType=Histo1D\n"
        s += "# ID\tID\tsumw\tsumw2\tsumwx\tsumwx2\tnumEntries\n"
        s += self.total.Format("Total")+"\n"
        s += self.uflow.Format("Underflow")+"\n"
        s += self.oflow.Format("Overflow")+"\n"
        s += "# xlow\txhigh\tsumw\tsumw2\tsumwx\tsumwx2\tnumEntries\n"
        for i in range(len(self.bins)):
            s += "{0}\n".format(self.bins[i])
        s += "END YODA_HISTO1D\n"
        return s

    def Fill(self,x,w):
        l = 0
        r = len(self.bins)-1
        c = int((l+r)/2)
        a = self.bins[c].xmin
        while r-l > 1:
            if x < a: r = c
            else: l = c
            c = int((l+r)/2)
            a = self.bins[c].xmin
        if x > self.bins[r].xmin:
            self.oflow.Fill(x,w)
        elif x < self.bins[l].xmin:
            self.uflow.Fill(x,w)
        else:
            self.bins[l].Fill(x,w)
        self.total.Fill(x,w)

    def ScaleW(self,scale):
        self.total.ScaleW(scale)
        self.uflow.ScaleW(scale)
        self.oflow.ScaleW(scale)
        for i in range(len(self.bins)):
            self.bins[i].ScaleW(scale)
        self.scale = scale


