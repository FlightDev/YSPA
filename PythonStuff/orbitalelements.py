from visual import *
from math import *

r = vector(0.59574, -1.70604, 0.30176)
r_dot = vector(-0.477734, -0.509794, 0.450681)
epsilon = 23.43687   # obliquity of the Ecliptic

r = r.rotate(-epsilon, vector(1, 0, 0))
r_dot = r_dot.rotate(-epsilon, vector(1, 0, 0))

h = cross(r, r_dot)
e = cross(r_dot, h) - norm(r)
q = mag(h) ** 2 / (1 + mag(e))
a = q/(1 - mag(e))
zhat = vector(0, 0, 1)
xhat = vector(1, 0, 0)
i = asin(h.z/mag(h))
N = cross(zhat, h)
Omega = acos(dot(xhat, N)/mag(N))
W = acos(dot(e, N)/ (mag(e) * mag(N)))

print "e", e
print "a", a
print "i", degrees(i)
print "OMEGA", degrees(Omega)
print "W", degrees(W)
