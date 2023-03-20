import math
from scipy.optimize import root

R = 0.083144
T = 328.15
Tc = 374.3
Pc = 40.65
w = 0.3268

a = 0.45724 * R**2 * Tc**2 / Pc * (1 + (0.37464 + 1.54226 * w - 0.26992 * w**2) * (1 - math.sqrt(T/Tc)))**2
b = 0.0778 * R * Tc / Pc
c = b * (1 + math.sqrt(2))
d = b * (1 - math.sqrt(2))

def eos(v, P):
    return R * T / (v - b) - a / ((v + c) * (v + d)) - P

def func(P):
    V = root(eos, [b+0.01, Pc-1], args=(P,), method='hybr', tol=1e-12).x
    VL, VG = min(V), max(V)
    ZG = P * VG / R / T
    ZL = P * VL / R / T
    F = (-1 + ZG - ((-c+d) * R * T * math.log(VG) + (c-d) * R * T * math.log(-b+VG) + a * (math.log(c+VG) - math.log(d+VG))) / ((c-d) * R * T) - math.log(ZG)) - (-1 + ZL - ((-c+d) * R * T * math.log(VL) + (c-d) * R * T * math.log(-b+VL) + a * (math.log(c+VL) - math.log(d+VL))) / ((c-d) * R * T) - math.log(ZL))
    return F

res = root(func, [1/10], method='hybr', tol=1e-12)
P = res.x[0]

V = root(eos, [b+0.01, Pc-1], args=(P,), method='hybr', tol=1e-12).x
VL, VG = min(V), max(V)

print("F = ", func(P))
print("P, VL, VG = ", P, VL, VG)
