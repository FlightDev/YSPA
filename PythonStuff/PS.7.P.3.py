import matplotlib.pyplot as plt
from visual import *
from math import *
import de421
from jplephem import Ephemeris

def find_rho(RA, Dec):
    RA = radians(RA)
    Dec = radians(Dec)
    z = sin(Dec)
    y = sin(RA) * cos(Dec)
    x = cos(RA) * cos(Dec)
    return vector(x, y, z)

def get_earth_to_sun(JD):
    eph = Ephemeris(de421)
    AU = 149597870.7 # km/AU
    barycenter = eph.position('earthmoon', JD)
    moonvector = eph.position('moon', JD)
    earth = barycenter - moonvector * eph.earth_share
    return -1.0 * vector(earth)/AU

def get_f_and_g(tau, r):
    #r = mag(r)
    f = 1 - ((tau ** 2) / (2 * (r ** 3)))
    g = tau - ((tau ** 3) / (6 * (r ** 3)))
    return f, g

def get_triple_products(rho1, rho2, rho3, a1, a3, R1, R2, R3):
    rho1 = norm(rho1)
    rho2 = norm(rho2)
    rho3 = norm(rho3)
    product1 = (a1 * (dot(cross(R1, rho2), rho3)) - dot(cross(R2, rho2), rho3) + a3 * (dot(cross(R3, rho2), rho3))) / (a1 * (dot(cross(rho1, rho2), rho3)))
    product2 = (a1 * (dot(cross(rho1, R1), rho3)) - dot(cross(rho1, R2), rho3) + a3 * (dot(cross(rho1, R3), rho3))) / (-1 * (dot(cross(rho1, rho2), rho3)))
    product3 = (a1 * (dot(cross(rho2, R1), rho1)) - dot(cross(rho2, R2), rho1) + a3 * (dot(cross(rho2, R3), rho1))) / (a3 * (dot(cross(rho2, rho3), rho1)))
    return product1, product2, product3

def accel(pos, mu):
    return -1.0 * mu * pos / (mag(pos) ** 3)

def verlet_position(current_pos, last_pos, d_tau, mu):
    r1 = 2.0 * current_pos - last_pos + accel(current_pos, mu) * (d_tau ** 2)
    return [r1, current_pos]

def location_to_angles(location):
    dec = asin(location.z)
    #print dec
    ra = acos(location.x/cos(dec))
    if location.y < 0:
        ra = 2 * pi - ra
    return degrees(ra), degrees(dec)

JDs = [2458309.750000,2458312.750000,2458315.750000,2458318.750000,2458321.750000,2458324.750000]
RAs =[331.7975,334.0454167,336.465,339.08667,341.94375,345.0720833]
Decs=[2.0456,3.5653,5.1900,6.9264,8.7783,10.7469]
Mags=[16.4,16.2,16.1,15.9,15.8,15.7]
day_in_july = [10, 13, 16, 19, 22, 25]
Rs = []

for JD in JDs:
    Rs.append(get_earth_to_sun(JD))

plt.scatter(Decs, RAs)
plt.xlabel("Declination (Degrees)")
plt.ylabel("Right Ascension (Degrees)")
plt.title("Right Ascension vs. Declination of 2061 Anza")
#plt.show() #TODO Uncomment Broski

k = 0.01720209895
epsilon = 23.43687   # obliquity of the Ecliptic


r2_estimate = 1.5

index1 = 1
index2 = 2
index3 = 5

#date in July, 2018(?)
t1 = day_in_july[index1]
t2 = day_in_july[index2]
t3 = day_in_july[index3]
tau1 = k * (t1 - t2)
tau2 = k * (t3 - t2)


rho1 = find_rho(RAs[index1], Decs[index1])
rho2 = find_rho(RAs[index2], Decs[index2])
rho3 = find_rho(RAs[index3], Decs[index3])
"""
rho1 = find_rho(RAs[index1], Decs[index1]).rotate(-epsilon, vector(1, 0, 0))
rho2 = find_rho(RAs[index2], Decs[index2]).rotate(-epsilon, vector(1, 0, 0))
rho3 = find_rho(RAs[index3], Decs[index3]).rotate(-epsilon, vector(1, 0, 0))

for i in range(len(Rs)):
    Rs[i] = Rs[i].rotate(-epsilon, vector(1, 0, 0))
"""
f1, g1 = get_f_and_g(tau1, r2_estimate)
f3, g3 = get_f_and_g(tau2, r2_estimate)

a1 = g3 / ((f1 * g3) - (f3 * g1))
a3 = -g1 / ((f1 * g3) - (f3 * g1))


product1, product2, product3 = get_triple_products(rho1, rho2, rho3, a1, a3, Rs[index1], Rs[index2], Rs[index3])
r2 = norm(rho2) * product2 - Rs[index2]
print "rs", Rs[index2]
print "HI, ", r2
print product1, product2, product3

r2 = vector(1.5, 0, 0)
for i in range(10):
    f1, g1 = get_f_and_g(tau1, mag(r2))
    f3, g3 = get_f_and_g(tau2, mag(r2))
    a1 = g3 / (f1 * g3 - f3 * g1)
    a3 = -g1 / (f1 * g3 - f3 * g1)
    product1, product2, product3 = get_triple_products(rho1, rho2, rho3, a1, a3, Rs[index1], Rs[index2], Rs[index3])
    r2 = norm(rho2) * product2 - Rs[index2]

"""
f1, g1 = get_f_and_g(tau1, mag(r2))
f3, g3 = get_f_and_g(tau2, mag(r2))
a1 = g3 / (f1 * g3 - f3 * g1)
a3 = -g1 / (f1 * g3 - f3 * g1)
product1, product2, product3 = get_triple_products(rho1, rho2, rho3, a1, a3, Rs[index1], Rs[index2], Rs[index3])
"""
r1 = norm(rho1) * product1 - Rs[index1]
r2 = norm(rho2) * product2 - Rs[index2]
r3 = norm(rho3) * product3 - Rs[index3]
r2_dot = ((f3) / (g1 * f3 - g3 * f1)) * r1 - ((f1)/(g1 * f3 - g3 * f1)) * r3
print "R2 ", r2
print "R2_DOT ", r2_dot


h = cross(r2, r2_dot)
e = cross(r2_dot, h) - norm(r2)
q = mag(h) ** 2 / (1 + mag(e))
a = q/(1-mag(e))
print "H:", h
print "E:", mag(e)
print "q:", q
print 'a:', a

find_index = 4
t = 0.00
tau = 0
end_t = day_in_july[find_index] - day_in_july[index2]
d_tau = 0.000001
k = 0.01720209895
end_tau = k * (end_t - t)
mu = 1.0

r_last = r2
r_dot = r2_dot

r_now = r_last + r_dot * d_tau + 0.5 * accel(r_last, mu) * (d_tau ** 2)
tau += d_tau

while tau <= end_tau:
    tau += d_tau
    temp = verlet_position(r_now, r_last, d_tau, mu)
    r_last = temp[1]
    r_now = temp[0]

AU = 149597870.7 # km/AU
jd = JDs[find_index]
eph = Ephemeris(de421)


R_geocentric = get_earth_to_sun(jd)
sun_to_asteroid = r_now
#sun_to_asteroid = sun_to_asteroid.rotate(radians(epsilon), vector(1, 0, 0))
print sun_to_asteroid
earth_to_asteroid = norm(R_geocentric + sun_to_asteroid) # Gets unit vector from Earth to Asteroid
print earth_to_asteroid
angles = location_to_angles(earth_to_asteroid)

ra, dec = angles[0], angles[1]
print angles[0], angles[1]
print "Residuals:"
print "RA: ", abs(ra - RAs[find_index])/RAs[find_index] * 100
print "DEC: ", abs(dec - Decs[find_index])/Decs[find_index] * 100
