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

def is_negative(num):
    if num < 0:
        return -1.0
    return 1.0




RAs = [  19.820671, 19.82062667, 19.82057944, 19.8207975, 19.82036417, 19.57390194, 18.161944, 18.16140639, 18.12, 18.119722, 18.11944  ,17.83268167, 17.83251944, 17.83234972, 17.83219167, 17.83204833, 17.831835]
Decs = [ 28.62544, 28.6259833, 28.6265722, 28.62676111, 28.62928056, 31.40670556, 35.92722, 35.92431944,  36.7675, 35.766, 35.7655 ,33.76921667, 33.76738056, 33.7656611, 33.76386944, 35.76268056, 33.76027222]
JDs = [ 2458281.856678, 2458281.857743, 2458281.858808, 2458281.859606, 2458281.863056, 2458286.836458, 2458315.567743, 2458315.569248, 2458316.711586, 2458316.718773, 2458316.725613 , 2458326.554572, 2458326.561562, 2458326.568299, 2458326.575336, 2458326.582002, 2458326.590035]


"""
JDs = [2458309.750000,2458312.750000,2458315.750000,2458318.750000,2458321.750000,2458324.750000]
RAs =[331.7975,334.0454167,336.465,339.08667,341.94375,345.0720833]
Decs=[2.0456,3.5653,5.1900,6.9264,8.7783,10.7469]
"""
#JDs = [10, 13, 16, 19, 22, 25]
Rs = []

#JDs = [2458281.856331, 2458315.564005, 2458306.711586, 2458325.554225]
#RAs =[297.3100, 272.4291, 274.7500, 267.4902]
#Decs=[28.6254, 35.9272, 35.7675, 33.7692]
Rs = []

for JD in JDs:
    Rs.append(get_earth_to_sun(JD))
"""
plt.scatter(Decs, RAs)
plt.xlabel("Declination (Degrees)")
plt.ylabel("Right Ascension (Degrees)")
plt.title("Right Ascension vs. Declination of 2061 Anza")
#plt.show() #TODO Uncomment Broski
"""

k = 0.01720209895
epsilon = 23.43687   # obliquity of the Ecliptic
r2_estimate = 1.4

index1 = 0
index2 = 7
index3 = 15

#date in July, 2018(?)
t1 = JDs[index1]
t2 = JDs[index2]
t3 = JDs[index3]
tau1 = k * (t1 - t2)
tau2 = k * (t3 - t2)

rho1 = find_rho(RAs[index1], Decs[index1])
rho2 = find_rho(RAs[index2], Decs[index2])
rho3 = find_rho(RAs[index3], Decs[index3])
f1, g1 = get_f_and_g(tau1, r2_estimate)
f3, g3 = get_f_and_g(tau2, r2_estimate)

a1 = g3 / ((f1 * g3) - (f3 * g1))
a3 = -g1 / ((f1 * g3) - (f3 * g1))


product1, product2, product3 = get_triple_products(rho1, rho2, rho3, a1, a3, Rs[index1], Rs[index2], Rs[index3])
r2 = norm(rho2) * product2 - Rs[index2]

r2 = vector(1.5, 0, 0)

for i in range(10):
    f1, g1 = get_f_and_g(tau1, mag(r2))
    f3, g3 = get_f_and_g(tau2, mag(r2))
    a1 = g3 / (f1 * g3 - f3 * g1)
    a3 = -g1 / (f1 * g3 - f3 * g1)
    product1, product2, product3 = get_triple_products(rho1, rho2, rho3, a1, a3, Rs[index1], Rs[index2], Rs[index3])
    r2 = norm(rho2) * product2 - Rs[index2]

#End of method of Gauss


r1 = norm(rho1) * product1 - Rs[index1]
r2 = norm(rho2) * product2 - Rs[index2]
r3 = norm(rho3) * product3 - Rs[index3]
r2_dot = ((f3) / (g1 * f3 - g3 * f1)) * r1 - ((f1)/(g1 * f3 - g3 * f1)) * r3

print 1234567890
print r2
print r2_dot
print 1234567890

find_index = 10
t = 0.00
tau = 0
end_t = JDs[find_index] - JDs[index2]
print is_negative(end_t)
d_tau = 0.000001 * is_negative(end_t)
k = 0.01720209895
end_tau = k * (end_t - t)
mu = 1.0

r_last = r2
r_dot = r2_dot

r_now = r_last + r_dot * d_tau + 0.5 * accel(r_last, mu) * (d_tau ** 2)
tau += d_tau

while abs(tau) <= abs(end_tau):
    tau += d_tau
    temp = verlet_position(r_now, r_last, d_tau, mu)
    r_last = temp[1]
    r_now = temp[0]

print r_now

AU = 149597870.7 # km/AU
jd = JDs[find_index]
eph = Ephemeris(de421)

R_geocentric = get_earth_to_sun(jd)
sun_to_asteroid = r_now
earth_to_asteroid = norm(R_geocentric + sun_to_asteroid) # Gets unit vector from Earth to Asteroid
angles = location_to_angles(earth_to_asteroid)

ra, dec = angles[0], angles[1]
print ra, dec
print RAs[find_index], Decs[find_index]
print "Residuals:"
print "RA: ", abs(ra - RAs[find_index])/RAs[find_index] * 100
print "DEC: ", abs(dec - Decs[find_index])/Decs[find_index] * 100
