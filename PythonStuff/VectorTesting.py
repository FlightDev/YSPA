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

def best_f_and_g(tau, r, rdot):
    r_mag = mag(r)
    f_1 = 1
    f_2 = (tau ** 2) / (2 * (r_mag ** 3))
    f_3 = dot(r, rdot) * (tau ** 3) / (2 * (r_mag ** 5))
    f_4 = ((tau ** 4) / 24) * ((3 / (r_mag ** 3)) * (((dot(rdot, rdot)) / (r_mag ** 2)) - ((1) / (r_mag ** 3))) - ((15 * (dot(r, rdot) ** 2)) / (r_mag ** 7)) + (1 / (r_mag ** 6)))

    g_1 = tau
    g_2 = (tau ** 3) / (6 * (r_mag ** 3))
    g_3 = dot(r, rdot) * (tau ** 4) / (4 * (r_mag ** 5))

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
JDs = []
RAs = []
Decs = []
file = open("data.txt", 'r')
data = []
for line in file:
    datum = line.split("\t")
    JDs.append(float(datum[0]))
    RAs.append(float(datum[1]))
    Decs.append(float(datum[2][:-1]))

#EPHEMERIS DATA
#JDs = [2458281.711585995, 2458286.711585995, 2458291.711585995, 2458296.711585995, 2458301.711585995, 2458306.711585995, 2458311.711585995, 2458316.711585995, 2458321.711585995, 2458326.711585995, 2458331.711585995]
#Decs = [270.85951, 273.58042, 276.21501, 278.77361, 281.26552, 283.69926, 286.08260, 288.42270, 290.72621, 292.99935, 295.24799]
#RAs = [-126.43995000000001, -110.90895, -95.5608, -80.4009, -65.43075, -50.6493, -36.0534, -21.6384, -7.398899999999999, 6.67125, 20.5785]

Rs = []
for JD in JDs:
    Rs.append(get_earth_to_sun(JD))

k = 0.01720209895
epsilon = 23.43687   # obliquity of the Ecliptic
r2_estimate = 1.5

index1 = 0
index2 = len(JDs) / 2
index3 = len(JDs) - 1
t1 = JDs[index1]
t2 = JDs[index2]
t3 = JDs[index3]
tau1 = k * (t1 - t2)
tau2 = k * (t3 - t2)

for i in range(50):
    f1, g1 = get_f_and_g(tau1, r2_estimate)
    f3, g3 = get_f_and_g(tau2, r2_estimate)
    a1 = g3 / (f1 * g3 - f3 * g1)
    a3 = -g1 / (f1 * g3 - f3 * g1)
    product1, product2, product3 = get_triple_products(rho1, rho2, rho3, a1, a3, Rs[index1], Rs[index2], Rs[index3])
    r2Gauss = norm(rho2) * product2 - Rs[index2]

#End of method of Gauss

r1 = norm(rho1) * product1 - Rs[index1]
r3 = norm(rho3) * product3 - Rs[index3]

r2_dotGauss = (f3 / (g1 * f3 - g3 * f1)) * r1 - (f1 / (g1 * f3 - g3 * f1)) * r3


#ITERATION:  20
#R2 <0.657047, -1.56004, -0.587581>
#R2_DOT <-0.433102, -0.252873, -0.597922>
Gauss_Decs = [-17.78225341312213, -17.782630768645543, -17.783856329813037, -17.78413896830974, -17.784798190833357, -17.785362940968643, -17.956520038070842, -17.957474637880008, -18.03945829321004, -18.040395100162538, -18.040971039497123, -18.184375917881816, -18.184550873825795, -18.185191939382907, -18.185948672292774, -18.186472000000002, -18.30464367939497, -18.324316704867236, -18.417531896381373, -18.41812571548247, -18.521895046716985, -18.522071891880376, -18.522223391728723, -18.522299113388055, -18.52240004631402, -18.53630706994539, -18.53630706994539, -18.536523937211772, -18.536572109095037]
Gauss_RAs = [276.3177997847802, 276.3128380717267, 276.296715334289, 276.2929953172769, 276.28431617513525, 276.2768779102753, 273.8836990497294, 273.8694809555755, 272.60532729840924, 272.5903575099203, 272.5811479607188, 270.1193234306484, 270.11608688264846, 270.10422215056406, 270.09020548613177, 270.080505, 267.72214376646593, 267.29205933452965, 265.0686258010063, 265.05335692846404, 262.13397372293934, 262.12856164425335, 262.1239240872669, 262.12160578221324, 262.11851519963545, 261.6883460025626, 261.6883460025626, 261.68157156075625, 261.68006651211016]

r2 = vector(0.657047, -1.56004, -0.587581)
r2_dot = vector(-0.433102, -0.252873, -0.597922)



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
