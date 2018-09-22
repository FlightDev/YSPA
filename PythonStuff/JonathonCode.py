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
    return -1.0 * vector(earth[0], earth[1], earth[2])/AU

def get_f_and_g(tau, r):
    #r = mag(r)
    f = 1 - ((tau * 2) / (2 * (mag(r) * 3)))
    g = tau - ((tau * 3) / (6 * (mag(r) * 3)))
    return f, g

def best_f_and_g(tau, r, rdot):
    r_mag = mag(r)
    f_1 = 1
    f_2 = (tau * 2) / (2 * (r_mag * 3))
    f_3 = dot(r, rdot) * (tau * 3) / (2 * (r_mag * 5))
    f_4 = ((tau * 4) / 24) * ((3 / (r_mag * 3)) * (((dot(rdot, rdot)) / (r_mag * 2)) - ((1) / (r_mag * 3))) - ((15 * (dot(r, rdot) * 2)) / (r_mag * 7)) + (1 / (r_mag ** 6)))

    g_1 = tau
    g_2 = (tau * 3) / (6 * (r_mag * 3))
    g_3 = dot(r, rdot) * (tau * 4) / (4 * (r_mag * 5))

    f = f_1 - f_2 + f_3 + f_4
    g = g_1 - g_2 + g_3
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

def find_endpoints(check_index, decs, ras, JDs, r, r_dot):
    dec_res_total = 0
    ra_res_total = 0
    total_ra = 0
    total_dec = 0
    times = []
    times.append(JDs[0])
    times.append(JDs[len(JDs) - 1])
    new_RAs = []
    new_decs = []
    for i in range(len(JDs)):
        k = 0.01720209895
        time_difference = JDs[i] - JDs[check_index]
        tau = time_difference * k

        f, g = best_f_and_g(tau, r, r_dot)
        r_now = f * r + g * r_dot

        AU = 149597870.7 #km/AU
        jd = JDs[i]
        eph = Ephemeris(de421)

        R_geocentric = get_earth_to_sun(JD)
        sun_to_asteroid = r_now
        earth_to_asteroid = norm(R_geocentric + sun_to_asteroid) #Gets unit vector from Earth to Asteroid
        RA, Dec = location_to_angles(earth_to_asteroid)
        new_RAs.append(RA)
        new_decs.append(Dec)
        """
        tau = 0
        d_tau = 0.000001
        k = 0.01720209895
        mu = 1.0
        """

    return new_RAs, new_decs

#REAL DATA

JDs = [2458313.643, 2458318.611, 2458318.615, 2458318.628, 2458318.631, 2458318.638, 2458318.644, 2458320.617, 2458320.629, 2458321.711, 2458321.724, 2458321.732, 2458323.939, 2458323.942, 2458323.953, 2458323.966, 2458323.975, 2458326.246, 2458326.68, 2458329.043, 2458329.06, 2458332.555, 2458332.562, 2458332.568, 2458332.571, 2458332.575, 2458333.139, 2458333.139, 2458333.148, 2458333.15]
RAs = [276.0694167, 273.062205, 273.060045, 273.050958, 273.049455, 273.045345, 273.041835, 271.906755, 271.908045, 271.295835, 271.288005, 271.283667, 270.098683, 270.097245, 270.092458, 270.08945, 270.080505, 268.71679, 268.716417, 267.614742, 267.606691, 266.189425, 266.137095, 266.1351, 266.132575, 266.132654, 265.909604, 265.908121, 265.906492, 265.905608]
Decs = [-12.340383, -15.204692, -15.207186, -15.203397, -15.217278, -15.220786, -15.224583, -16.342956, -16.341075, -16.952889, -16.960417, -16.964778, -18.167397, -18.168828, -18.173417, -18.176403, -18.186472, -19.616666, -19.616944, -20.8106, -20.818869, -22.507581, -22.507433, -22.509272, -22.510661, -22.510661, -22.770333, -22.774228, -22.775256, -22.775256]

#EPHEMERIS DATA
#JDs = [2458281.711585995, 2458286.711585995, 2458291.711585995, 2458296.711585995, 2458301.711585995, 2458306.711585995, 2458311.711585995, 2458316.711585995, 2458321.711585995, 2458326.711585995, 2458331.711585995]
#Decs = [270.85951, 273.58042, 276.21501, 278.77361, 281.26552, 283.69926, 286.08260, 288.42270, 290.72621, 292.99935, 295.24799]
#RAs = [-126.43995000000001, -110.90895, -95.5608, -80.4009, -65.43075, -50.6493, -36.0534, -21.6384, -7.398899999999999, 6.67125, 20.5785]

Rs = []
for JD in JDs:
    Rs.append(get_earth_to_sun(JD))

k = 0.01720209895
epsilon = 23.43687   # obliquity of the Ecliptic
r2_estimate = vector(1.5, 0, 0)

index1 = 0
index2 = len(JDs) / 2
index3 = len(JDs) - 1
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

for i in range(1000):
    f1, g1 = get_f_and_g(tau1, r2)
    f3, g3 = get_f_and_g(tau2, r2)
    a1 = g3 / (f1 * g3 - f3 * g1)
    a3 = -g1 / (f1 * g3 - f3 * g1)
    product1, product2, product3 = get_triple_products(rho1, rho2, rho3, a1, a3, Rs[index1], Rs[index2], Rs[index3])
    r2 = norm(rho2) * product2 - Rs[index2]

#End of method of Gauss



r1 = norm(rho1) * product1 - Rs[index1]
r2 = norm(rho2) * product2 - Rs[index2]
r3 = norm(rho3) * product3 - Rs[index3]

r2_dotGauss = ((f3) / (g1 * f3 - g3 * f1)) * r1 - ((f1)/(g1 * f3 - g3 * f1)) * r3
r2Gauss = r2

#ITERATION:  20
#R2 <0.657047, -1.56004, -0.587581>
#R2_DOT <-0.433102, -0.252873, -0.597922>
Gauss_Decs = []
Gauss_RAs = []

r2 = vector(0.657047, -1.56004, -0.587581)
r2_dot = vector(-0.433102, -0.252873, -0.597922)

#<0.528811, -1.49139, -0.572998>
#<0.466179, 0.316431, -0.270101>

r2Gauss = vector(0.528811, -1.49139, -0.572998)
r2_dotGauss = vector(0.466179, 0.316431, -0.270101)

Gauss_RAs, Gauss_Decs = find_endpoints(index2, Decs, RAs, JDs, r2Gauss, r2_dotGauss)
optimized_RAs, optimized_Decs = find_endpoints(index2, Decs, RAs, JDs, r2, r2_dot)

f,(ax1, ax2, ax3) = plt.subplots(3, sharex = False)

ax1.set_title("RA vs JD")
ax2.set_title("DEC vs JD")
ax3.set_title("DEC vs RA")

ax1.set_xlabel("Julian Date")
ax1.set_ylabel("RA (degrees)")

ax2.set_xlabel("Julian Date")
ax2.set_ylabel("Declination (degrees)")

ax3.set_xlabel("RA (degrees)")
ax3.set_ylabel("Declination (degrees)")


ax1.plot(JDs, Gauss_RAs, color = 'green')
ax1.plot(JDs, optimized_RAs, color = 'red')
ax1.scatter(JDs, RAs, color = 'blue')

ax2.plot(JDs, Gauss_Decs, color = 'green')
ax2.plot(JDs, optimized_Decs, color = 'red')
ax2.scatter(JDs, Decs, color = 'blue')

ax3.plot(Gauss_RAs, Gauss_Decs, color = 'green')
ax3.plot(optimized_RAs, optimized_Decs, color = 'red')
ax3.scatter(RAs, Decs, color = 'blue')

plt.gca().invert_xaxis()


plt.tight_layout(h_pad = 1.50)
plt.show()