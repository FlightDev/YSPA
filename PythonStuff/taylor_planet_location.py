from math import *
from visual import *

r = vector(-0.00879000, -1.1700000, 0.12900000)
rog = r
rdot = vector(1.0600000, 0.9200000, -0.12400000)
AU_to_meters = 1.496E11

start_solar_date = 0.00000
end_solar_date = 5.00000
k = 0.01720209895
tau = start_solar_date * k
tau_step = 1.0E-6
end_tau = end_solar_date * k
d_tau = end_tau - tau
print d_tau

def f(r, rdot, d_tau):
    fact_1 = 1
    fact_2 = (d_tau ** 2) / (2 * mag(r) ** 3)
    fact_3 = (dot(r, rdot) * d_tau ** 3)/(2 * mag(r) ** 5)
    return fact_1 - fact_2 + fact_3

def g(r, rdot, d_tau):
    fact_1 = d_tau
    fact_2 = (d_tau ** 3)/(6 * mag(r) ** 3)
    fact_3 = (dot(r, rdot) * d_tau ** 4)/(4 * mag(r) ** 5)
    return fact_1 - fact_2 + fact_3

def accel(r):
    return -1.0 * 1.0 * r /(mag(r) ** 3)

taylor_pos = f(r, rdot, d_tau) * r + g(r, rdot, d_tau) * rdot

rnew = r + rdot * tau_step + 0.5 * accel(r) * (tau_step ** 2)

while tau <= end_tau:
    rtemp = rnew
    rnew = 2.0 * rnew - r + accel(rnew) * tau_step ** 2
    r = rtemp
    tau += tau_step

last_vect = r - taylor_pos

print taylor_pos
print r
print mag(last_vect)/mag(rog)*AU_to_meters
