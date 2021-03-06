

# This file was *autogenerated* from the file prog.sage
from sage.all_cmdline import *   # import sage library

_sage_const_149 = Integer(149); _sage_const_2 = Integer(2); _sage_const_1 = Integer(1); _sage_const_3 = Integer(3); _sage_const_0 = Integer(0); _sage_const_4 = Integer(4); _sage_const_6 = Integer(6); _sage_const_12 = Integer(12); _sage_const_5 = Integer(5); _sage_const_20 = Integer(20); _sage_const_8 = Integer(8); _sage_const_16 = Integer(16)
import sys
from sage.all import *

selectedPrime = _sage_const_149 

a = _sage_const_2 
b = _sage_const_1 


if len(sys.argv) != 4:
    print('need 4 arguments [a,b,prime]')
    exit()

a = Integer(sys.argv[1])
b = Integer(sys.argv[2])
selectedPrime = Integer(sys.argv[3])



R = PolynomialRing(GF(selectedPrime),_sage_const_2 ,order='invlex', names=('x', 'y',)); (x, y,) = R._first_ngens(2)
Q = QuotientRing(R,R.ideal(y**_sage_const_2  - x**_sage_const_3  -a*x -b), names=('x', 'y',)); (x, y,) = Q._first_ngens(2)
S = PolynomialRing(GF(selectedPrime), names=('o',)); (o,) = S._first_ngens(1)#okruh pro testovani ireducibility

def expMod(poly, exp, modu):
    S = QuotientRing(Q, Q.ideal(modu), names=('s', 't',)); (s, t,) = S._first_ngens(2)#vytvor okruh K[x,y]/(y**2-x**3-a*x-b, modu)
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

    if m == _sage_const_0 :
        toReturn = _sage_const_0 
    elif m == _sage_const_1  or m == _sage_const_2 :
        toReturn = _sage_const_1 
    elif m == _sage_const_3 :
        toReturn = _sage_const_3 *x**_sage_const_4 +_sage_const_6 *a*x**_sage_const_2 +_sage_const_12 *b*x-a**_sage_const_2 
    elif m == _sage_const_4 :
        toReturn = _sage_const_2 *(x**_sage_const_6  + _sage_const_5 *a*x**_sage_const_4  + _sage_const_20 *b*x**_sage_const_3  - _sage_const_5 *a*a*x**_sage_const_2  - _sage_const_4 *a*b*x - _sage_const_8 *b**_sage_const_2  - a**_sage_const_3 )
    else:
        if m % _sage_const_2  == _sage_const_0 :
            mHalf = m // _sage_const_2 
            toReturn = getFPoly(mHalf) * (getFPoly(mHalf + _sage_const_2 ) * getFPoly(mHalf - _sage_const_1 )**_sage_const_2  - getFPoly(mHalf - _sage_const_2 ) * getFPoly(mHalf + _sage_const_1 )**_sage_const_2 )
        else:
            mHalf = (m-_sage_const_1 )//_sage_const_2 
            if mHalf % _sage_const_2  == _sage_const_0 :
                toReturn = (_sage_const_16 *(x**_sage_const_3 +a*x+b)**_sage_const_2 )*getFPoly(mHalf+_sage_const_2 )*getFPoly(mHalf)**_sage_const_3  - getFPoly(mHalf-_sage_const_1 )*getFPoly(mHalf+_sage_const_1 )**_sage_const_3 
            else:
                toReturn = getFPoly(mHalf+_sage_const_2 )*getFPoly(mHalf)**_sage_const_3  -(_sage_const_16 *(x**_sage_const_3 +a*x+b)**_sage_const_2 )*getFPoly(mHalf-_sage_const_1 )*getFPoly(mHalf+_sage_const_1 )**_sage_const_3 
    fPolyCache[m] = toReturn
    return fPolyCache[m]

#funkce pocita polynomy \bar{s}_m 
def getSBarPoly(m):
    ql = selectedPrime % m
    fl = getFPoly(m)
    if ql == _sage_const_1 :
        return expMod(x, Integer(selectedPrime**_sage_const_2 ), fl)-x
    if ql % _sage_const_2  == _sage_const_0 :
        return _sage_const_4 *(expMod(x, Integer(selectedPrime**_sage_const_2 ), fl)-x)*(x**_sage_const_3 +a*x+b)*getFPoly(ql)**_sage_const_2 +getFPoly(ql-_sage_const_1 )*getFPoly(ql+_sage_const_1 )
    else:
        return (expMod(x, Integer(selectedPrime**_sage_const_2 ), fl)-x)*getFPoly(ql)**_sage_const_2 +(_sage_const_4 *(x**_sage_const_3 +a*x+b))*getFPoly(ql-_sage_const_1 )*getFPoly(ql+_sage_const_1 )
    

#funkce pocita polynomy r_m
def getRPoly(m):
    if m in RPolyCache:
        return RPolyCache[m]

    if m == _sage_const_1 :
        toReturn = _sage_const_1 
    elif m % _sage_const_2  == _sage_const_0 :
        toReturn = getFPoly(m+_sage_const_2 )*getFPoly(m-_sage_const_1 )**_sage_const_2  - getFPoly(m-_sage_const_2 )*getFPoly(m+_sage_const_1 )**_sage_const_2 
    else:
        toReturn = getFPoly(m+_sage_const_2 )*getFPoly(m-_sage_const_1 )**_sage_const_2  - getFPoly(m-_sage_const_2 )*getFPoly(m+_sage_const_1 )**_sage_const_2 

    RPolyCache[m] = toReturn
    return RPolyCache[m]

#funkce pocita polynomy s_m
def getSPoly(m):
    if m in SPolyCache:
        return SPolyCache[m]
    if m == _sage_const_1 :
        toReturn = _sage_const_1 
    elif m % _sage_const_2  == _sage_const_0 :
        toReturn = _sage_const_16 *(x**_sage_const_3  + a*x + b)**_sage_const_2  * getFPoly(m)**_sage_const_3 
    else:
        toReturn = getFPoly(m)**_sage_const_3 

    SPolyCache[m] = toReturn
    return SPolyCache[m]

#funkce pocita polynomy c_m
def getCPoly(m):
    if m in CPolyCache:
        return CPolyCache[m]

    if m == _sage_const_1 :
        toReturn = _sage_const_0 
    elif m % _sage_const_2  == _sage_const_0 :
        toReturn = getFPoly(m-_sage_const_1 )*getFPoly(m+_sage_const_1 )
    else:
        toReturn = _sage_const_4 *(x**_sage_const_3  + a*x + b)*getFPoly(m-_sage_const_1 )*getFPoly(m+_sage_const_1 )

    CPolyCache[m] = toReturn
    return CPolyCache[m]

#funkce pocita polynomy d_m
def getDPoly(m):
    if m in DPolyCache:
        return DPolyCache[m]

    if m == _sage_const_1 :
        toReturn = _sage_const_1 
    elif m % _sage_const_2  == _sage_const_0 :
        toReturn = _sage_const_4 *(x**_sage_const_3  + a*x + b)*getFPoly(m)**_sage_const_2 
    else:
        toReturn = getFPoly(m)**_sage_const_2 

    DPolyCache[m] = toReturn
    return DPolyCache[m]

#kontrola jestli polyToTest je násobek f_l
def tyzero(l, m):
    fl = getFPoly(l)
    exponent = (selectedPrime**_sage_const_2 -_sage_const_1 ) / _sage_const_2 
    polyToTest = getSPoly(m) * expMod((x**_sage_const_3 +a*x+b), Integer(exponent), fl) + getRPoly(m)
    res = R(polyToTest) % R(fl)
    return res == _sage_const_0 

#gl = gcd(sl, fl)
def eigen(l, gamma, gl):
    if gl == _sage_const_1 :
        return False

    fl = getFPoly(l)


    exponent = (selectedPrime-_sage_const_1 ) / _sage_const_2 
    polyToTest = getSPoly(Integer(gamma)) * expMod((x**_sage_const_3  + a*x + b), Integer(exponent), fl) - getRPoly(Integer(gamma))
    res = R(polyToTest) % R(gl)

    return res == _sage_const_0 
    
def equalx(l, gl):
    K = GF(l) 
    ql = selectedPrime % l

    if tyzero(l, ql): #true jestlize phi^2(P) = [-ql](P) pro P z E[l]*, poté t_l = 0
        return _sage_const_0 

    #nyní phi^2(P) = [ql](P) pro nějaké P z E[l]*

    tau = (K(_sage_const_4 *ql).sqrt())

    #2 kandidáti na t_l => +-tau

    gamma = (K(_sage_const_2 *ql) * K(tau)**(-_sage_const_1 ))

    if eigen(l, gamma, gl): #kontrola jestli phi(P) = [gamma](P) pro nějaké P z E[l]*, když false musí být phi(P) = [-gamma](P)
        return tau
    return -tau
    
def nonequalx(l, tau):

    m = selectedPrime % l

    fl = getFPoly(l)

    if m > _sage_const_0 : #pokud je m=0 (nastane v pripade kdy q=selectedPrime je p^e, e > 1)
        cm = R(getCPoly(m))
        dm = R(getDPoly(m))
        rm = R(getRPoly(m))
        sm = R(getSPoly(m))

    ctau = R(getCPoly(tau))
    dtau = R(getDPoly(tau))
    rtau = R(getRPoly(tau))
    stau = R(getSPoly(tau))

    F = Q.fraction_field()

    exponent = (selectedPrime**_sage_const_2  - _sage_const_1 )/_sage_const_2 

    if m > _sage_const_0 :
        lam = y * (F(dm)/F(sm)) * (F(expMod(x**_sage_const_3  + a*x + b, Integer(exponent), fl) * sm - rm))/(F(dm * (expMod(x,Integer(selectedPrime**_sage_const_2 ), fl) - x) + cm))

        firstX = F(lam**_sage_const_2  - expMod(x,Integer(selectedPrime**_sage_const_2 ), fl) - x) + (F(cm)/F(dm))
    else:
        firstX = expMod(x,Integer(selectedPrime**_sage_const_2 ), fl)

    secondX = F(expMod(x, Integer(selectedPrime), fl)) - (F(ctau(x=expMod(x, Integer(selectedPrime), fl)))/F(dtau(x=expMod(x, Integer(selectedPrime), fl))))

    hx = R((firstX - secondX).numerator()) % R(fl)

    firstGCD =  R(fl).gcd(hx)

    if firstGCD == _sage_const_1 :
        return _sage_const_0 

    if m > _sage_const_0 :
        firstY = lam * (_sage_const_2  * expMod(x, Integer(selectedPrime**_sage_const_2 ), fl) - lam**_sage_const_2  + x - F(cm)/F(dm)) - expMod(y, Integer(selectedPrime**_sage_const_2 ), fl)
    else:
        firstY = expMod(y,Integer(selectedPrime**_sage_const_2 ), fl)

    secondY = expMod(y, Integer(selectedPrime), fl) * F(rtau(x=expMod(x, Integer(selectedPrime), fl)))/F(stau(x=expMod(x, Integer(selectedPrime), fl)))

    hy = R((firstY - secondY).numerator())(y=_sage_const_1 ) % R(fl)

    secondGCD = R(fl).gcd(hy)

    if secondGCD == _sage_const_1 :
        return -_sage_const_1 
    return _sage_const_1 

def schoff():
    B = _sage_const_2 
    l = _sage_const_2 
    curvePoly = o**_sage_const_3 +a*o+b
    if curvePoly.is_irreducible() == _sage_const_1 : #t_l = 0 <=> E[F_q] má involuci <=> x**3+a*x+b má kořen v F_q
        tau = _sage_const_1 
    else:
        tau = _sage_const_0 
    residues = [Integer(tau)]
    moduli = [Integer(l)]
    totalMod = l
    
    while B < _sage_const_4 *sqrt(selectedPrime):
        l = l.next_prime()
        totalMod *= l
        B = B*l
        fl = getFPoly(l)

        if selectedPrime % l == _sage_const_0 : #ql = 0, phi^2 nikdy není zero isogeny.
            sl = _sage_const_1 
        else:
            sl = R(getSBarPoly(l))

        sl = sl % R(fl)
        gl = R(fl).gcd(sl) #gl != 1 <==> existuje P v E[l]* t.ž. phi^2(P) = [+-ql](P)

        if gl != _sage_const_1 : #víme, že phi^2(P) = [+-ql](P) pro nějaké P v E[l]*
            tau = equalx(l, gl)
        else: #phi^2(P) != [+-ql](P) pro každé P z E[l]*, můžeme tedy použít vzorec pro součet různých bodů
            r = _sage_const_0 
            tau = _sage_const_0 
            while r == _sage_const_0 : #dokud r není +-1
                tau += _sage_const_1 
                r = nonequalx(l, tau) #testuje phi^2 + [ql] = [tau] pro různá tau > 0, tl != 0 jelikož phi^2(P) != [+-ql](P)
            if r == -_sage_const_1 :
                tau = -tau
        
        residues.append(Integer(tau))
        moduli.append(Integer(l))

    trace = crt(residues, moduli)
    if trace >= (totalMod/_sage_const_2 ): #stopa je v Z_totalMod, chceme ji převést do Z
        trace -= totalMod

    return selectedPrime+_sage_const_1 -trace

print(schoff())

