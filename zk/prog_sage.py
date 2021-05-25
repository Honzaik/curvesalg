from sympy import *

selectedPrime = 149
K = GF(selectedPrime)
x, y = symbols('x y')
a = K(1)
b = K(1)


f = poly(4*x**2 + 1, domain=K)
g = poly(8*x**3 + 1, domain=K)

gc = gcd(f,g)
print(rem(f,g,domain=K))
fPolyCache = {}

def getFPoly(m):
    toReturn = None
    if m in fPolyCache:
        return fPolyCache[m]

    if m == 0:
        toReturn = poly(0, x, domain=K)
    elif m == 1 or m == 2:
        toReturn = poly(1, x, domain=K)
    elif m == 3:
        toReturn = poly(3*x**4+6*a*x**2+12*b*x-a**2, domain=K)
    elif m == 4:
        toReturn = poly(2*(x**6 + 5*a*x**4 + 20*b*x**3 - 5*a*a*x**2 - 4*a*b*x - 8*b**2 - a**3), domain=K)
    else:
        if m % 2 == 0:
            mHalf = m // 2
            toReturn = getFPoly(mHalf) * (getFPoly(mHalf + 2) * getFPoly(mHalf - 1)**2 - getFPoly(mHalf - 2) * getFPoly(mHalf + 1)**2)
        else:
            mHalf = (m-1)//2
            if mHalf % 2 == 0:
                toReturn = poly(16*(x**3+a*x+b)**2, domain=K)*getFPoly(mHalf+2)*getFPoly(mHalf)**3 - getFPoly(mHalf-1)*getFPoly(mHalf+1)**3
            else:
                toReturn = getFPoly(mHalf+2)*getFPoly(mHalf)**3 - poly(16*(x**3+a*x+b)**2, domain=K)*getFPoly(mHalf-1)*getFPoly(mHalf+1)**3
    fPolyCache[m] = toReturn
    return fPolyCache[m]

def getSPoly(m):
    ql = selectedPrime % m
    if ql == 1:
        return poly(x**(selectedPrime**2)-x, domain=K)
    if ql % 2 == 0:
        return poly(4*(x**(selectedPrime**2)-x)*(x**3+a*x+b), domain=K)*getFPoly(ql)**2+getFPoly(ql-1)*getFPoly(ql+1)
    else:
        return poly(x**(selectedPrime**2)-x, domain=K)*getFPoly(ql)**2+poly(4*(x**3+a*x+b), domain=K)*getFPoly(ql-1)*getFPoly(ql+1)

def schoff():
    B = 2
    l = 2
    curvePoly = poly(x**3+a*x+b, domain=K)
    cycloPoly = poly(x**selectedPrime - x, domain=K)
    if (gcd(curvePoly, cycloPoly) == 1) :
        tau = 1
    else:
        tau = 0
    M = [(2,tau)]
    
    while B < 4*sqrt(selectedPrime):
        l = nextprime(l)
        B = B*l
        fl = getFPoly(l)
        sl = getSPoly(l)
        print(fl)
        

schoff()