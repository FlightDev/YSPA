import matplotlib.pyplot as plt
from visual import *
from math import *
import de421
from jplephem import Ephemeris
import numpy as np
import rebound

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

    date = "2018-7-16 01:52"


    for i in range(len(decs)):
        #print str(i + 1) + "/" + str(len(decs))
        k = 0.01720209895
        time_difference = JDs[i] - JDs[check_index]
        tau = 0
        d_tau = 0.001 * is_negative(time_difference)
        end_tau = time_difference * k

        #TEMP INTEGRATOR THINGY HERE (THE NEW ONE BRO)
        sim = rebound.Simulation()
        date = "2018-7-16 01:52"

        #TEST INFO#

        #SUN
        sim.add(m = 1.0)

        #ASTEROID
        sim.add(m = 0.0, x = r.x, y =r.y , z = r.z, vx = r_dot.x, vy = r_dot.y, vz = r_dot.z)

        #EARTH
        sim.add(m = 1.0/330000.0, x = 4.189675889417898E-01, y =-8.496104448760934E-01, z =-3.683119473905782E-01, vx = 1.539665060981763E-02/k, vy= 6.454046892673179E-03/k, vz = 2.797224517462366E-03/k)
        #sim.add('Earth')

        #JUPITER
        sim.add(m = 1.0/1000.0, x =-3.213494013743870E+00, y =-4.009890793133180E+00, z =-1.640521133996669E+00, vx = 5.974697377415830E-03/k, vy =-3.756443410778172E-03/k, vz =-1.755593829452442E-03/k)
        #sim.add('Jupiter')

        #SATURN
        sim.add(m = 1.0/3300.0,  x = 1.084878715668408E+00, y =-9.231862109791130E+00, z =-3.860012135382278E+00, vx= 5.246951399009119E-03/k, vy= 6.191309101028240E-04, vz= 3.018304002418789E-05/k)
        #sim.add('Saturn')

        sim.move_to_com()

        ps = sim.particles

        sim.integrate(end_tau)

        r_now = vector(ps[1].x, ps[1].y, ps[1].z)

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

    sigma = sqrt(ra_res_total + dec_res_total)
    return sigma

def RAs_and_DECs(check_index, JDs, r, r_dot):
    ras = []
    decs = []

    for i in range(len(JDs)):
        #print str(i + 1) + "/" + str(len(decs))
        k = 0.01720209895
        time_difference = JDs[i] - JDs[check_index]
        tau = 0
        d_tau = 0.001 * is_negative(time_difference)
        end_tau = time_difference * k

        #TEMP INTEGRATOR THINGY HERE (THE NEW ONE BRO)
        sim = rebound.Simulation()

        #TEST INFO#

        #SUN
        sim.add(m = 1.0)

        #ASTEROID
        sim.add(m = 0.0, x = r.x, y =r.y , z = r.z, vx = r_dot.x, vy = r_dot.y, vz = r_dot.z)

        #EARTH
        sim.add(m = 1.0/330000.0, x = 4.189675889417898E-01, y =-8.496104448760934E-01, z =-3.683119473905782E-01, vx = 1.539665060981763E-02/k, vy= 6.454046892673179E-03/k, vz = 2.797224517462366E-03/k)
        #sim.add('Earth')

        #JUPITER
        sim.add(m = 1.0/1000.0, x =-3.213494013743870E+00, y =-4.009890793133180E+00, z =-1.640521133996669E+00, vx = 5.974697377415830E-03/k, vy =-3.756443410778172E-03/k, vz =-1.755593829452442E-03/k)
        #sim.add('Jupiter')

        #SATURN
        sim.add(m = 1.0/3300.0,  x = 1.084878715668408E+00, y =-9.231862109791130E+00, z =-3.860012135382278E+00, vx= 5.246951399009119E-03/k, vy= 6.191309101028240E-04, vz= 3.018304002418789E-05/k)
        #sim.add('Saturn')

        sim.move_to_com()

        ps = sim.particles

        sim.integrate(end_tau)

        r_now = vector(ps[1].x, ps[1].y, ps[1].z)

        AU = 149597870.7 # km/AU
        jd = JDs[i]
        eph = Ephemeris(de421)

        R_geocentric = get_earth_to_sun(JD)
        sun_to_asteroid = r_now
        earth_to_asteroid = norm(R_geocentric + sun_to_asteroid) # Gets unit vector from Earth to Asteroid
        RA, Dec = location_to_angles(earth_to_asteroid)
        print "HI"
        ras.append(RA)
        decs.append(Dec)

    return ras, decs





#REALISH DATA
#RAs = [  19.820671, 19.82062667, 19.82057944, 19.8207975, 19.82036417, 19.57390194, 18.161944, 18.16140639, 18.12, 18.119722, 18.11944  ,17.83268167, 17.83251944, 17.83234972, 17.83219167, 17.83204833, 17.831835]
#Decs = [ 28.62544, 28.6259833, 28.6265722, 28.62676111, 28.62928056, 31.40670556, 35.92722, 35.92431944,  36.7675, 35.766, 35.7655 ,33.76921667, 33.76738056, 33.7656611, 33.76386944, 35.76268056, 33.76027222]
#June 11, 2018:July
2458281.856678, 2458315.711586
JDs = [2458281.856678, 2458281.857743, 2458281.858808, 2458281.859606, 2458281.863056, 2458286.836458, 2458315.564236, 2458315.569248, 2458315.578218, 2458315.711586, 2458316.718773, 2458316.725613, 2458316.727662, 2458326.554572, 2458326.561562, 2458326.568299, 2458326.575336, 2458326.582002, 2458326.590035]
RAs = [297.31006665, 297.30940004999997, 297.3086916, 297.3077958, 297.30546255, 293.6085291, 272.429166, 272.42109585, 272.4214917, 271.8, 271.7958333, 271.7919291, 271.7908584, 267.49022505, 267.4877916, 267.4852458, 267.48287505, 267.4872495, 267.477525]
Decs = [28.625444, 28.62598333, 28.62657222, 28.62676111, 28.62928056, 31.40670556, 35.9272222, 35.92431944, 35.92481944, 35.7675, 35.7666, 35.76640278, 35.76515278, 33.76921667, 33.76738056, 33.76566111, 33.76386944, 33.76268056, 33.760272222]

#EPHEM DATA
#RAs = [2458281.711585995, 2458286.711585995, 2458291.711585995, 2458296.711585995, 2458301.711585995, 2458306.711585995, 2458311.711585995, 2458316.711585995, 2458321.711585995, 2458326.711585995, 2458331.711585995]
#Decs = [270.85951, 273.58042, 276.21501, 278.77361, 281.26552, 283.69926, 286.08260, 288.42270, 290.72621, 292.99935, 295.24799]
#JDs = [-8.42933, -7.39393, -6.37072, -5.36006, -4.36205, -3.37662, -2.40356, -1.44256, -0.49326, 0.44475, 1.37190]


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
index2 = 9
index3 = 17

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
#r2 = r2.rotate(-epsilon, vector(1, 0, 0))
#r2_dot = r2_dot.rotate(-epsilon, vector(1, 0, 0))
print r2
print r2_dot
"""
 <0.608031, -1.59416, 0.226923>
VELOCITY:  <-0.580291, -0.0752849, 0.15723>

k = 0.01720209895
(0.4334311, -1.30895, -0.0372353)
(0.01124734, 0.00144496115, 0.004523545)
"""
r2 = vector(0.608031, -1.59416, 0.226923)
r2_dot = vector(-0.580291, -0.0752849, 0.15723)
"""
old_sigma = find_residuals(index2, Decs, RAs, JDs, r2, r2_dot)
r2_temp = vector(r2.x, r2.y, r2.z)
r2_dot_temp = vector(r2_dot.x, r2_dot.y, r2_dot.z)
i = 1
pass_counter = 0.0
last_changed = 0.0
#(check_index, decs, ras, JDs, r, r_dot)
while True:
    print "ITERATION: ", i
    change_stuff = 0.01
    r2 += vector(np.random.normal(0.0, change_stuff, 3))
    r2_dot += vector(np.random.normal(0.0, change_stuff, 3))
    if mag(r2) > 5 or mag(r2_dot) > 5:
        r2 = r2_temp
        r2_dot = r2_dot_temp
        continue
    sigma = find_residuals(index2, Decs, RAs, JDs, r2, r2_dot)

    if sigma > old_sigma:
        r2 = vector(r2_temp.x, r2_temp.y, r2_temp.z)
        r2_dot = vector(r2_dot_temp.x, r2_dot_temp.y, r2_dot_temp.z)
        sigma = old_sigma
        print "FAILED"
    else:
        last_changed = i
        pass_counter += 1.0
        print "PASSED"
    r2_temp = vector(r2.x, r2.y, r2.z)
    r2_dot_temp = vector(r2_dot.x, r2_dot.y, r2_dot.z)
    old_sigma = sigma
    print "POSITION: ", r2
    print "VELOCITY: ", r2_dot
    print "CHI SQUARED: ", sigma
    print str(pass_counter * 1.0 / i * 100) + "% PASSED"
    i += 1
    if i - last_changed == 32:
        change_stuff /= 10
        print "CHANGE NOW SET TO ", change_stuff
"""

ras, decs = RAs_and_DECs(index2, JDs, r2, r2_dot)
print len(ras), len(decs), len(JDs)

f,(ax1, ax2) = plt.subplots(2, sharex = True)

ax1.plot(JDs, ras)
ax2.plot(JDs, decs)
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
