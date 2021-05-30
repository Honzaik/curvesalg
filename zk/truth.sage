from sage.all import *

a = 1
b = 1
p = 191


curve = EllipticCurve(GF(p),[a,b])
print(curve.cardinality())