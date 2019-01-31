import numpy as np

from vector import Vec4

def u_2(r):
    """
    Integration variable u_2 --- solution for r = 2u^1 - 1u^2
    """
    return 1. -  np.sqrt(1.-r)

def u_3(r):
    """
    Integration variable u_3 --- solution for r = 3u^2 - 2u^3
    """
    x = pow(1.-2.*r+2.*np.sqrt(r*(r-1.)+0.j),1./3.)
    y = (2.-(1.-1.j*np.sqrt(3.))/x-(1.+1.j*np.sqrt(3.))*x)/4.
    return y.real

def u_4(r):
    """
    Integration variable u_4 --- solution for r = 4u^3 - 3u^4
    """
    y = pow(r+np.sqrt(r*r*(1-r)+0.j),1./3.)
    x = 3./2.*(r/y+y)
    y = np.sqrt(1.+x)
    z = (1.+y-np.sqrt(2.-x+2./y))/3.
    return z.real

def f(x, a, r):
    """
    The equation ax^(a-1) - (a-1)x^a - r = 0
    To be used as argument in solver
    """
    return a*x**(a-1) - (a-1)*x**a - r

def fp(x, a):
    """
    First derivative of f
    """
    return a*(a-1)*(x**(a-2) - x**(a-1))

def fpp(x, a):
    """
    Second derivative of
    """
    return a*(a-1)*((a-2)*x**(a-3) - (a-1)*x**(a-2))

def get_u(a, r):
    """
    Solve f for u

    a = n + 1 -i in Simon's notation

    The lowest order case is n=3 and i = 2, i.e. a = 2
    """
    if a < 2 : raise Exception("a = {} not implemented".format(a))

    from scipy import optimize
    if a == 2: return u_2(r)
    elif a == 3: return u_3(r)
    elif a == 4: return u_4(r)
    else:
        return optimize.newton(lambda x : f(x, a, r), r, fprime=lambda x: fp(x,a), fprime2=lambda x: fpp(x,a))

def rho(Min, Mout, mext=0.0):
    """
    Helper function for mass term eq (5)
    """
    M2 = Min*Min
    return 0.125 * np.sqrt( (M2 - (Mout+mext)*(Mout+mext)) * (M2 - (Mout-mext)*(Mout-mext))) / M2

def rho_massless(Min, Mout):
    """
    Helper function for mass term eq (5)
    """
    M2  = Min*Min
    M22 = Mout*Mout
    return 0.125 * np.sqrt( M2*M2 - 2*M2*M22 + M22*M22) / M2

def generate_point(pa,pb,rans):

    # The final momenta
    MOM = [ -rans[-1]*pa, -rans[-2]*pb ]
    _Q  = -MOM[0]-MOM[1]

    # Storage of intermediate Masses, Qs
    M = [_Q.M()]
    ALLQ =[_Q]

    U, R = [], [] # Store the u and random numbers r
    for i in range(2, NP+1):
        # print("now i = {}".format(i))
        if i < NP:
            # print("Solving for u_{}, M_{}".format(i, i))
            r = rans[3*(i-2)+2]
            u = get_u(NP+1-i, r)
            U.append(u)
            R.append(r)
            # Simon's paper must be wrong here, check
            _M = np.sqrt(u*_Q.M2()) # M_i^2
        else:
            _M = 0
        # print("Got M_{}={}".format(i, _M))
        M.append(_M)

        q = 4*M[-2] * rho_massless(M[-2], M[-1])
        # Random numbers for costheta and phi
        costheta = 2*rans[3*(i-2)] - 1
        phi = 2.*np.pi*rans[3*(i-2)+1]

        # Generated 4 Vectors
        # p_(i-1)
        sintheta = np.sqrt(1. - costheta*costheta)
        p = q*Vec4(1, np.cos(phi)*sintheta, np.sin(phi)*sintheta, costheta)
        # print("p_{} = {} {}".format(i+1, p, np.sqrt(abs(p.M2()))))

        # now Q_i
        _Q = Vec4(np.sqrt(q*q + M[-1]*M[-1]), -p.px, -p.py, -p.pz)
        # print("Q_{} = {} {}".format(i, _Q, np.sqrt(abs(_Q.M2()))))

        p = ALLQ[i-2].BoostBack(p)
        _Q = ALLQ[i-2].BoostBack(_Q)
        # print ALLQ[i-2]-_Q-p
        # print "p boosted ",p,p.M2()
        # print "Q boosted ",_Q,np.sqrt(abs(_Q.M2()))
        # print "sum p+Q   ",(p+_Q),(p+_Q).M()
        MOM.append(p)
        ALLQ.append(_Q)
    MOM.append(_Q)
    return MOM

def generate_weight(pa,pb,mom):
    Q = -mom[0]-mom[1]
    rans = []
    for i in range(2, NP+1):
        # print("now i = {}".format(i))
        p = Q.Boost(mom[i])
        # print 'p = ',p
        costh = p[3]/p.P()
        phi = p.Phi()
        if phi < 0: phi += 2.*np.pi
        # print "phi = ",phi
        rans.append((1+costh)/2.)
        rans.append(phi/(2.*np.pi))
        if i < NP:
            m = (Q-mom[i]).M2() / Q.M2()
            u = f(m, NP+1-i, 0)
            # print Q.M2(),(Q-mom[i]).M2(),(mom[3]+mom[4]).M2(),m,u
            # print Q
            Q -= mom[i]
            # print Q
            rans.append(u)
        else:
            _M = 0
    rans.append(-(mom[1]*pa)/(pa*pb))
    rans.append(-(mom[0]*pb)/(pa*pb))
    return rans

def ME_ESW(P):
    """
    Calculate the matrix element for g(p1) g(p2) --> g(p3) g(p4) g(p5)

    Using eq (7.51) in QCD for collider physics.

    P ... list of 4 momentum vectors
    """
    from itertools import permutations
    permutations=list(permutations([0,1,2,3,4])) # All 120 permutations

    # M = const * A * B / C

    # A = sum_permutations {1 2} ^ 4
    A = 0
    B = 0
    for i in permutations:
        A+= (P[i[0]] * P[i[1]])**4
        B+= (P[i[0]] * P[i[1]]) * (P[i[1]] * P[i[2]]) * (P[i[2]] * P[i[3]]) * (P[i[3]] * P[i[4]]) * (P[i[4]] * P[i[0]])

    C = 1
    for i in range(5):
        for j in range(5):
            if i <j:
                # print("i={}, j={}: {} * {} = {}".format(i, j, P[i], P[j], P[i]*P[j]))
                C *= P[i]*P[j]

    return A*B/C

def ME_PLB(P):
    """
    Calculate the matrix element for g(p1) g(p2) --> g(p3) g(p4) g(p5)

    Using eq (18) in Berends et al, Phys Let B 103 (1981) p 124 ff.

    P ... list of 4 momentum vectors
    """
    from itertools import permutations, combinations
    permutations= [
            (0,1,2,3,4),
            (0,1,2,4,3),
            (0,1,3,2,4),
            (0,1,3,4,2),
            (0,1,4,2,3),
            (0,1,4,3,2),
            (0,2,1,3,4),
            (0,2,1,4,3),
            (0,2,3,1,4),
            (0,2,4,1,3),
            (0,3,1,2,4),
            (0,3,2,1,4),
            ]

    kpermutations = list(combinations([0,1,2,3,4], 2))

    # M = const * A * B / C

    # A = sum_permutations {1 2} ^ 4
    A = 0
    for i in kpermutations:
        A+= (P[i[0]] * P[i[1]])**4

    B = 0
    for i in permutations:
        # print("(k{} * k{})^4".format(i[0]+1, i[1]+1))
        B+= (P[i[0]] * P[i[1]]) * (P[i[1]] * P[i[2]]) * (P[i[2]] * P[i[3]]) * (P[i[3]] * P[i[4]]) * (P[i[4]] * P[i[0]])

    C = 1
    for i in range(5):
        for j in range(5):
            if i <j:
                # print("i={}, j={}: {} * {} = {}".format(i, j, P[i], P[j], P[i]*P[j]))
                C *= P[i]*P[j]

    return A*B/C
if __name__ == "__main__":

    import sys

    np.random.seed(1)
    pa = Vec4(7000,0,0,7000)
    pb = Vec4(7000,0,0,-7000)

    if len(sys.argv) <2:
        print("Please specify the number of external particles, exiting")
        sys.exit(1)

    NP = int(sys.argv[1]) # Number of external particles
    if NP<3:
        print("NP should be >=3 for the whole thing to make sense, exiting")
        sys.exit(1)

    rans = [ np.random.rand() for i in range(0,3*NP-4+2) ]

    moms = generate_point(pa,pb,rans)

    msum = Vec4()
    for num, mom in enumerate(moms):
        msum += mom
        print("p_{} = {} {}".format(num+1, mom, mom.M2()))
    print("Mom sum {}".format(msum))

    ranc = generate_weight(pa,pb,moms)

    for r in range(0,len(rans)):
        print("r_{} = {} -> dev. {}".format(r, ranc[r], ranc[r]/rans[r]-1))


    print("120*Berends: {:.20f}".format(120*ME_PLB(moms)))
    print("Ellis:       {:.20f}".format(ME_ESW(moms)))

