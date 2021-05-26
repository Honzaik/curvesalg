from sage.all import *
import sys

sys.setrecursionlimit(2000)


selectedPrime = 149

a = 1
b = 1

R.<x,y> = PolynomialRing(GF(selectedPrime),2,order='invlex')
Q.<x,y> = QuotientRing(R,R.ideal(y**2 - x**3 -a*x -b))


def expMod(poly, exp, modu):
    S.<s,t> = QuotientRing(Q, Q.ideal(modu)) #exp
    return Q(S(poly)**exp)

fPolyCache = {}
RPolyCache = {}
SPolyCache = {}
CPolyCache = {}
DPolyCache = {}

def getFPoly(m):
    toReturn = None
    if m in fPolyCache:
        return fPolyCache[m]

    if m == 0:
        toReturn = 0
    elif m == 1 or m == 2:
        toReturn = 1
    elif m == 3:
        toReturn = 3*x**4+6*a*x**2+12*b*x-a**2
    elif m == 4:
        toReturn = 2*(x**6 + 5*a*x**4 + 20*b*x**3 - 5*a*a*x**2 - 4*a*b*x - 8*b**2 - a**3)
    else:
        if m % 2 == 0:
            mHalf = m // 2
            toReturn = getFPoly(mHalf) * (getFPoly(mHalf + 2) * getFPoly(mHalf - 1)**2 - getFPoly(mHalf - 2) * getFPoly(mHalf + 1)**2)
        else:
            mHalf = (m-1)//2
            if mHalf % 2 == 0:
                toReturn = (16*(x**3+a*x+b)**2)*getFPoly(mHalf+2)*getFPoly(mHalf)**3 - getFPoly(mHalf-1)*getFPoly(mHalf+1)**3
            else:
                toReturn = getFPoly(mHalf+2)*getFPoly(mHalf)**3 -(16*(x**3+a*x+b)**2)*getFPoly(mHalf-1)*getFPoly(mHalf+1)**3
    fPolyCache[m] = toReturn
    return fPolyCache[m]

def getSBarPoly(m):
    ql = selectedPrime % m
    if ql == 1:
        return x**(selectedPrime**2)-x
    if ql % 2 == 0:
        return 4*(x**(selectedPrime**2)-x)*(x**3+a*x+b)*getFPoly(ql)**2+getFPoly(ql-1)*getFPoly(ql+1)
    else:
        return (x**(selectedPrime**2)-x)*getFPoly(ql)**2+(4*(x**3+a*x+b))*getFPoly(ql-1)*getFPoly(ql+1)
    

def getRPoly(m):
    if m in RPolyCache:
        return RPolyCache[m]

    if m == 1:
        toReturn = 1
    elif m % 2 == 0:
        toReturn = getFPoly(m+2)*getFPoly(m-1)**2 - getFPoly(m-2)*getFPoly(m+1)**2
    else:
        toReturn = getFPoly(m+2)*getFPoly(m-1)**2 - getFPoly(m-2)*getFPoly(m+1)**2

    RPolyCache[m] = toReturn
    return RPolyCache[m]

def getSPoly(m):
    if m in SPolyCache:
        return SPolyCache[m]
    if m == 1:
        toReturn = 1
    elif m % 2 == 0:
        toReturn = 16*(x**3 + a*x + b)**2 * getFPoly(m)**3
    else:
        toReturn = getFPoly(m)**3

    SPolyCache[m] = toReturn
    return SPolyCache[m]

def getCPoly(m):
    if m in CPolyCache:
        return CPolyCache[m]

    if m == 1:
        toReturn = 0
    elif m % 2 == 0:
        toReturn = getFPoly(m-1)*getFPoly(m-1)
    else:
        toReturn = 4*(x**3 + a*x + b)*getFPoly(m-1)*getFPoly(m-1)

    CPolyCache[m] = toReturn
    return CPolyCache[m]

def getDPoly(m):
    if m in DPolyCache:
        return DPolyCache[m]

    if m == 1:
        toReturn = 1
    elif m % 2 == 0:
        toReturn = 4*(x**3 + a*x + b)*getFPoly(m)**2
    else:
        toReturn = getFPoly(m)**2

    DPolyCache[m] = toReturn
    return DPolyCache[m]

def tyzero(l, m):
    fl = getFPoly(l)
    exponent = (selectedPrime**2-1) / 2
    polyToTest = getSPoly(m) * expMod((x**3+a*x+b), Integer(exponent), fl) + getRPoly(m)
    res = R(polyToTest) % R(fl)
    return res == 0
    
def eigen(l, gamma, gl):
    if gl == 1:
        return False

    exponent = (selectedPrime**2-1) / 2
    polyToTest = getSPoly(m) * expMod((x**3+a*x+b), Integer(exponent), fl) - getRPoly(m)
    res = R(polyToTest) % R(gl)
    print('eigen res:', res)
    return res == 0
    
def equalx(l, gl):
    K = GF(l) 
    ql = selectedPrime % l
    print(K(ql).is_square())
    if not K(ql).is_square() or tyzero(l, ql):
        return 0

    tau = (R(4*ql).sqrt()) % l
    gamma = (R(2*ql) * R(tau)**(-1)) % l
    if eigen(l, gamma, gl):
        return tau
    return -tau
    
def nonequalx(l, tau):

    m = selectedPrime % l

    fl = getFPoly(l)
    cm = R(getCPoly(m)) % R(fl)
    dm = R(getDPoly(m)) % R(fl)
    rm = R(getRPoly(m)) % R(fl)
    sm = R(getSPoly(m)) % R(fl)

    ctau = R(getCPoly(tau)) % R(fl)
    dtau = R(getDPoly(tau)) % R(fl)
    rtau = R(getRPoly(tau)) % R(fl)
    stau = R(getSPoly(tau)) % R(fl)

    F = Q.fraction_field()

    exponent = (selectedPrime**2 - 1)/2

    lam = y * (F(dm)/F(sm)) * (F(expMod(x**3 + a*x + b, Integer(exponent), fl) * sm - rm))/(F(dm * (expMod(x,Integer(selectedPrime**2), fl) - x) + cm))

    firstX = F(lam**2 - expMod(x,Integer(selectedPrime**2), fl) - x) + (F(cm)/F(dm))
    secondX = F(x**selectedPrime) - (F(ctau(x=x**selectedPrime))/F(dtau(x=x**selectedPrime)))


    hx = R((firstX - secondX).numerator()) % R(fl)
    if R(fl).gcd(R(hx)) == 1:
        return 0


    return 0
    hy = 0 % R(fl)
    
    
    if R(fl).gcd(R(hy)) == 1:
        return -1
    return 1

def schoff():
    B = 2
    l = 2
    curvePoly = x**3+a*x+b
    cycloPoly = x**selectedPrime - x
    if R(curvePoly).gcd(R(cycloPoly)) == 1:
        tau = 1
    else:
        tau = 0
    residues = [tau]
    moduli = [2]
    
    while B < 4*sqrt(selectedPrime):
        l = l.next_prime()
        B = B*l
        fl = getFPoly(l)
        sl = R(getSBarPoly(l))

        sl = sl % R(fl)
        gl = R(fl).gcd(sl)

        if gl == 1:
            print('equalx')
            tau = equalx(l, gl)
        else:
            print('not equal x')
            r = None
            tau = 0
            while r != 0:
                tau += 1
                r = nonequalx(l, tau)
            if r == -1:
                tau = -tau
        residues.append(tau)
        moduli.append(l)

    trace = crt(residues, moduli)
    return selectedPrime+1-trace

print(schoff())

print(tyzero(3, 2))