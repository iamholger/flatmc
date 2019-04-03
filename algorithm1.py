import numpy as np

from vector import Vec4

def u_2(r):
    """
    Integration variable u_2 --- solution for r = 2u^1 - 1u^2
    """
    return 1. -  np.sqrt(1.-r)

# def u3(r):
    # kappa = np.cbrt(  2*np.sqrt(r*(r-1) ) -2*r +1  )
    # return 0.5*(kappa + 1./kappa + 1)

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

    pa = Vec4(45.,0,0, 45)
    pb = Vec4(45.,0,0,-45)

    p1 = Vec4(30, 4, 5, 1)
    p2 = Vec4(30, 4, 5, -1)
    p3 = Vec4(30, 4, -5, 1)

    TEST = [pa,pb,p1,p2,p3]
    print("120*Berends: {}".format(120*ME_PLB(TEST)))
    print("Ellis:       {}".format(ME_ESW(TEST)))

    import time
    NEVAL=10000

    t1=time.time()
    for _ in range(NEVAL):ME_PLB(TEST)
    t2=time.time()

    print("{} evaluations of PLB took {} seconds".format(NEVAL, t2-t1))
    t1=time.time()
    for _ in range(NEVAL):ME_ESW(TEST)
    t2=time.time()

    print("{} evaluations of ESW took {} seconds".format(NEVAL, t2-t1))
    from IPython import embed
    embed()


    import sys

    if len(sys.argv) <2:
        print("Please specify the number of external particles, exiting")
        sys.exit(1)

    NP = int(sys.argv[1]) # Number of external particles
    if NP<3:
        print("NP should be >=3 for the whole thing to make sense, exiting")
        sys.exit(1)

    # This is the 4-Vector of the initial blob of energy
    Q  = Vec4(100,0,0,10)
    Q1 = Q
    Qabs = np.sqrt(Q*Q)
    M1 = Qabs

    # Storage of intermediate Masses, Qs
    M = [M1]
    ALLQ =[Q1]

    # The final momenta
    MOM = []

    U, R = [], [] # Store the u and random numbers r
    for i in range(2, NP+1):
        # print("now i = {}".format(i))
        if i < NP:
            # print("Solving for u_{}, M_{}".format(i, i))
            r = np.random.rand()
            u = get_u(NP+1-i, r)
            U.append(u)
            R.append(r)
            _M = np.prod(U)*Qabs # M_i
        else:
            _M = 0
        # print("Got M_{}={}".format(i, _M))
        M.append(_M)

        q = 4*M[-2] * rho_massless(M[-2], M[-1])
        # Random numbers for costheta and phi
        costheta = 2*np.random.rand() -1
        phi = np.pi*np.random.rand()

        # Generated 4 Vectors
        # p_(i-1)
        p = q*Vec4(1, np.cos(phi)*np.sqrt(1. - costheta*costheta), np.sin(phi)*np.sqrt(1. - costheta*costheta), costheta)
        # print("p_{} = {}".format(i-1, p))

        # now Q_i
        _Q = Vec4(np.sqrt(q*q + M[-1]*M[-1]), -p.px, -p.py, -p.pz)
        ALLQ.append(_Q)
        # print("Q_{} = {}".format(i, _Q))

        # TODO ask Stefan how boost is supposed to work with his class
        # boostVector = Vec4(ALLQ[-1][0]/ALLQ[-1].E, ALLQ[-1][1]/ALLQ[-1].E, ALLQ[-1][2]/ALLQ[-1].E, ALLQ[-1][3]/ALLQ[-2].E)
        boostVector = ALLQ[0]-ALLQ[-1]

        p.BoostBack(boostVector)
        MOM.append(p)
        _Q.BoostBack(boostVector)
    MOM.append(_Q)


    for num, mom in enumerate(MOM):
        print("p_{} = {}".format(num+1, mom))



    from IPython import embed
    embed()
