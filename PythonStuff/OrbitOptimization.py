import matplotlib.pyplot as plt
from visual import *
from math import *
import de421
from jplephem import Ephemeris
import numpy as np

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

def find_residuals(check_index, decs, ras, JDs, r, r_dot):
    dec_res_total = 0
    ra_res_total = 0
    total_ra = 0
    total_dec = 0
    for i in range(len(decs)):
        time_difference = JDs[i] - JDs[check_index]
        tau = 0
        d_tau = 0.0001 * is_negative(time_difference)
        k = 0.01720209895
        mu = 1.0
        end_tau = time_difference * k
        r_last = r
        r_now = r_last + r_dot * d_tau + 0.5 * accel(r_last, mu) * (d_tau ** 2)
        tau += d_tau

        while abs(tau) <= abs(end_tau):
            tau += d_tau
            temp = verlet_position(r_now, r_last, d_tau, mu)
            r_last = temp[1]
            r_now = temp[0]
        AU = 149597870.7 # km/AU
        jd = JDs[i]
        eph = Ephemeris(de421)

        R_geocentric = get_earth_to_sun(JD)
        sun_to_asteroid = r_now
        earth_to_asteroid = norm(R_geocentric + sun_to_asteroid) # Gets unit vector from Earth to Asteroid
        RA, Dec = location_to_angles(earth_to_asteroid)
        RA_Residual = abs(ras[i] - RA)
        Dec_Residual = abs(decs[i] - Dec)
        ra_res_total += RA_Residual ** 2
        dec_res_total += Dec_Residual ** 2

        tau = 0
        d_tau = 0.000001
        k = 0.01720209895
        mu = 1.0

    sigma = sqrt((ra_res_total + dec_res_total))
    return sigma








#REALISH DATA
RAs = [  19.820671, 19.82062667, 19.82057944, 19.8207975, 19.82036417, 19.57390194, 18.161944, 18.16140639, 18.12, 18.119722, 18.11944  ,17.83268167, 17.83251944, 17.83234972, 17.83219167, 17.83204833, 17.831835]
Decs = [ 28.62544, 28.6259833, 28.6265722, 28.62676111, 28.62928056, 31.40670556, 35.92722, 35.92431944,  36.7675, 35.766, 35.7655 ,33.76921667, 33.76738056, 33.7656611, 33.76386944, 35.76268056, 33.76027222]
JDs = [ 2458281.856678, 2458281.857743, 2458281.858808, 2458281.859606, 2458281.863056, 2458286.836458, 2458315.567743, 2458315.569248, 2458316.711586, 2458316.718773, 2458316.725613 , 2458326.554572, 2458326.561562, 2458326.568299, 2458326.575336, 2458326.582002, 2458326.590035]

Rs = []
#TEST DATA
#JDs = [2458309.750000,2458312.750000,2458315.750000,2458318.750000,2458321.750000,2458324.750000]
#RAs =[331.7975,334.0454167,336.465,339.08667,341.94375,345.0720833]
#Decs=[2.0456,3.5653,5.1900,6.9264,8.7783,10.7469]

#JDs = [2458309.750000,2458312.750000,2458315.750000,2458318.750000]
#RAs = [331.7975,334.0454167,336.465,339.08667]
#Decs = [2.0456,3.5653,5.1900,6.9264]


for JD in JDs:
    Rs.append(get_earth_to_sun(JD))

k = 0.01720209895
epsilon = 23.43687   # obliquity of the Ecliptic
r2_estimate = 1.5

index1 = 0
index2 = 2
index3 = 5

#date in July, 2018(?)
t1 = JDs[index1]
t2 = JDs[index2]
t3 = JDs[index3]
tau1 = k * (t1 - t2)
tau2 = k * (t3 - t2)

rho1 = find_rho(RAs[index1], Decs[index1])
rho2 = find_rho(RAs[index2], Decs[index2])
rho3 = find_rho(RAs[index3], Decs[index3])

#rho1 = rho1.rotate(radians(epsilon), vector(1, 0, 0))
#rho2 = rho2.rotate(radians(epsilon), vector(1, 0, 0))
#rho3 = rho3.rotate(radians(epsilon), vector(1, 0, 0))

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
r2 = r2.rotate(-epsilon, vector(1, 0, 0))
r2_dot = r2_dot.rotate(-epsilon, vector(1, 0, 0))
print r2
print r2_dot

old_sigma = find_residuals(index2, Decs, RAs, JDs, r2, r2_dot)
r2_temp = r2
r2_dot_temp = r2_dot
i = 1
#(check_index, decs, ras, JDs, r, r_dot)
while True:
    print "ITERATION: ", i
    r2 += vector(np.random.normal(0.0, 0.0001, 3))
    r2_dot += vector(np.random.normal(0.0, 0.0001, 3))
    sigma = find_residuals(index2, Decs, RAs, JDs, r2, r2_dot)

    if sigma > old_sigma:
        r2 = r2_temp
        r2_dot = r2_dot_temp
        sigma = old_sigma
        print "FAILED"

    r2_temp = r2
    r2_dot_temp = r2_dot
    old_sigma = sigma

    print "POSITION: ", r2
    print "VELOCITY: ", r2_dot
    print "SIGMA: ", sigma
    i += 1

print "Residuals:"
print "RA: ", ra
print "DEC: ", dec

f,(ax1, ax2) = plt.subplots(2, sharex = True)

ax1.plot(times, ras)
ax2.plot(times, decs)
ax1.scatter(JDs, RAs)
ax2.scatter(JDs, Decs)
plt.show()

"""
0
<0.613031, -1.04491, -0.370265>
1
<0.647705, -1.00807, -0.35579>
2
<0.681496, -0.969855, -0.34083>
3
<0.714319, -0.930266, -0.325387>
4
<0.746086, -0.889305, -0.309463>
5
<0.776707, -0.846977, -0.293064>
"""
