import sys
from sage.all import *

selectedPrime = 149

a = 2
b = 1

'''
if len(sys.argv) != 4:
    print('need 4 arguments [a,b,prime]')
    exit()

a = Integer(sys.argv[1])
b = Integer(sys.argv[2])
selectedPrime = Integer(sys.argv[3])
'''


R.<x,y> = PolynomialRing(GF(selectedPrime),2,order='invlex')
Q.<x,y> = QuotientRing(R,R.ideal(y**2 - x**3 -a*x -b))
S.<o> = PolynomialRing(GF(selectedPrime)) #okruh pro testovani ireducibility

def expMod(poly, exp, modu):
    S.<s,t> = QuotientRing(Q, Q.ideal(modu)) #vytvor okruh K[x,y]/(y**2-x**3-a*x-b, modu)
    return Q(S(poly)**exp) #spocitej exp v něm

fPolyCache = {}
RPolyCache = {}
SPolyCache = {}
CPolyCache = {}
DPolyCache = {}

#funkce pocita polynomy \bar{f}_m
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

#funkce pocita polynomy \bar{s}_m 
def getSBarPoly(m):
    ql = selectedPrime % m
    fl = getFPoly(m)
    if ql == 1:
        return expMod(x, Integer(selectedPrime**2), fl)-x
    if ql % 2 == 0:
        return 4*(expMod(x, Integer(selectedPrime**2), fl)-x)*(x**3+a*x+b)*getFPoly(ql)**2+getFPoly(ql-1)*getFPoly(ql+1)
    else:
        return (expMod(x, Integer(selectedPrime**2), fl)-x)*getFPoly(ql)**2+(4*(x**3+a*x+b))*getFPoly(ql-1)*getFPoly(ql+1)
    

#funkce pocita polynomy r_m
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

#funkce pocita polynomy s_m
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

#funkce pocita polynomy c_m
def getCPoly(m):
    if m in CPolyCache:
        return CPolyCache[m]

    if m == 1:
        toReturn = 0
    elif m % 2 == 0:
        toReturn = getFPoly(m-1)*getFPoly(m+1)
    else:
        toReturn = 4*(x**3 + a*x + b)*getFPoly(m-1)*getFPoly(m+1)

    CPolyCache[m] = toReturn
    return CPolyCache[m]

#funkce pocita polynomy d_m
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

#kontrola jestli polyToTest je násobek f_l
def tyzero(l, m):
    fl = getFPoly(l)
    exponent = (selectedPrime**2-1) / 2
    polyToTest = getSPoly(m) * expMod((x**3+a*x+b), Integer(exponent), fl) + getRPoly(m)
    res = R(polyToTest) % R(fl)
    return res == 0

#gl = gcd(sl, fl)
def eigen(l, gamma, gl):
    if gl == 1:
        return False

    fl = getFPoly(l)


    exponent = (selectedPrime-1) / 2
    polyToTest = getSPoly(Integer(gamma)) * expMod((x**3 + a*x + b), Integer(exponent), fl) - getRPoly(Integer(gamma))
    res = R(polyToTest) % R(gl)

    return res == 0
    
def equalx(l, gl):
    K = GF(l) 
    ql = selectedPrime % l

    if tyzero(l, ql): #true jestlize phi^2(P) = [-ql](P) pro P z E[l]*, poté t_l = 0
        return 0

    #nyní phi^2(P) = [ql](P) pro nějaké P z E[l]*

    tau = (K(4*ql).sqrt())

    #2 kandidáti na t_l => +-tau

    gamma = (K(2*ql) * K(tau)**(-1))

    if eigen(l, gamma, gl): #kontrola jestli phi(P) = [gamma](P) pro nějaké P z E[l]*, když false musí být phi(P) = [-gamma](P)
        return tau
    return -tau
    
def nonequalx(l, tau):

    m = selectedPrime % l

    fl = getFPoly(l)

    if m > 0: #pokud je m=0 (nastane v pripade kdy q=selectedPrime je p^e, e > 1)
        cm = R(getCPoly(m))
        dm = R(getDPoly(m))
        rm = R(getRPoly(m))
        sm = R(getSPoly(m))

    ctau = R(getCPoly(tau))
    dtau = R(getDPoly(tau))
    rtau = R(getRPoly(tau))
    stau = R(getSPoly(tau))

    F = Q.fraction_field()

    exponent = (selectedPrime**2 - 1)/2

    if m > 0:
        lam = y * (F(dm)/F(sm)) * (F(expMod(x**3 + a*x + b, Integer(exponent), fl) * sm - rm))/(F(dm * (expMod(x,Integer(selectedPrime**2), fl) - x) + cm))

        firstX = F(lam**2 - expMod(x,Integer(selectedPrime**2), fl) - x) + (F(cm)/F(dm))
    else:
        firstX = expMod(x,Integer(selectedPrime**2), fl)

    secondX = F(expMod(x, Integer(selectedPrime), fl)) - (F(ctau(x=expMod(x, Integer(selectedPrime), fl)))/F(dtau(x=expMod(x, Integer(selectedPrime), fl))))

    hx = R((firstX - secondX).numerator()) % R(fl)

    firstGCD =  R(fl).gcd(hx)

    if firstGCD == 1:
        return 0

    if m > 0:
        firstY = lam * (2 * expMod(x, Integer(selectedPrime**2), fl) - lam**2 + x - F(cm)/F(dm)) - expMod(y, Integer(selectedPrime**2), fl)
    else:
        firstY = expMod(y,Integer(selectedPrime**2), fl)

    secondY = expMod(y, Integer(selectedPrime), fl) * F(rtau(x=expMod(x, Integer(selectedPrime), fl)))/F(stau(x=expMod(x, Integer(selectedPrime), fl)))

    hy = R((firstY - secondY).numerator())(y=1) % R(fl)

    secondGCD = R(fl).gcd(hy)

    if secondGCD == 1:
        return -1
    return 1

def schoff():
    B = 2
    l = 2
    curvePoly = o**3+a*o+b
    if curvePoly.is_irreducible() == 1: #t_l = 0 <=> E[F_q] má involuci <=> x**3+a*x+b má kořen v F_q
        tau = 1
    else:
        tau = 0
    residues = [Integer(tau)]
    moduli = [Integer(l)]
    totalMod = l
    
    while B < 4*sqrt(selectedPrime):
        l = l.next_prime()
        totalMod *= l
        B = B*l
        fl = getFPoly(l)

        if selectedPrime % l == 0: #ql = 0, phi^2 nikdy není zero isogeny.
            sl = 1
        else:
            sl = R(getSBarPoly(l))

        sl = sl % R(fl)
        gl = R(fl).gcd(sl) #gl != 1 <==> existuje P v E[l]* t.ž. phi^2(P) = [+-ql](P)

        if gl != 1: #víme, že phi^2(P) = [+-ql](P) pro nějaké P v E[l]*
            tau = equalx(l, gl)
        else: #phi^2(P) != [+-ql](P) pro každé P z E[l]*, můžeme tedy použít vzorec pro součet různých bodů
            r = 0
            tau = 0
            while r == 0: #dokud r není +-1
                tau += 1
                r = nonequalx(l, tau) #testuje phi^2 + [ql] = [tau] pro různá tau > 0, tl != 0 jelikož phi^2(P) != [+-ql](P)
            if r == -1:
                tau = -tau
        
        residues.append(Integer(tau))
        moduli.append(Integer(l))

    trace = crt(residues, moduli)
    if trace >= (totalMod/2): #stopa je v Z_totalMod, chceme ji převést do Z
        trace -= totalMod

    return selectedPrime+1-trace

print(schoff())