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

def find_residuals(check_index, decs, ras, JDs, r, r_dot):
    dec_res_total = 0
    ra_res_total = 0
    total_ra = 0
    total_dec = 0

    for i in range(len(decs)):
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

for i in range(1000):
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
print r2
print r2_dot

#ITERATION:  20
#R2 <0.657047, -1.56004, -0.587581>
#R2_DOT <-0.433102, -0.252873, -0.597922>
#r2 = vector(0.657047, -1.56004, -0.587581)
#r2_dot = vector(-0.433102, -0.252873, -0.597922)

old_sigma = find_residuals(index2, Decs, RAs, JDs, r2, r2_dot)
r2_temp = vector(r2.x, r2.y, r2.z)
r2_dot_temp = vector(r2_dot.x, r2_dot.y, r2_dot.z)

print "R2:", r2_temp
print "R2_DOT:", r2_dot_temp

#FINDS X DIRECTION
temp_x = vector(r2_temp.x + 0.01, r2_temp.y, r2_temp.z)
res_x_1 = find_residuals(index2, Decs, RAs, JDs, r2_temp, r2_dot_temp)
res_x_2 = find_residuals(index2, Decs, RAs, JDs, temp_x, r2_dot)
x_slope = is_negative(res_x_2 - res_x_1)

#FINDS Y DIRECTION
temp_y = vector(r2_temp.x, r2_temp.y + 0.01, r2_temp.z)
res_y_1 = find_residuals(index2, Decs, RAs, JDs, r2_temp, r2_dot_temp)
res_y_2 = find_residuals(index2, Decs, RAs, JDs, temp_y, r2_dot)
y_slope = is_negative(res_y_2 - res_y_1)

#FINDS Z DIRECTION
temp_z = vector(r2_temp.x, r2_temp.y, r2_temp.z + 0.01)
res_z_1 = find_residuals(index2, Decs, RAs, JDs, r2_temp, r2_dot_temp)
res_z_2 = find_residuals(index2, Decs, RAs, JDs, temp_z, r2_dot)
z_slope = is_negative(res_z_2 - res_z_1)

#FINDS X_DOT DIRECTION
temp_x = vector(r2_dot_temp.x + 0.01, r2_dot_temp.y, r2_dot_temp.z)
res_x_1 = find_residuals(index2, Decs, RAs, JDs, r2_temp, r2_dot_temp)
res_x_2 = find_residuals(index2, Decs, RAs, JDs, r2_temp, temp_x)
x_dot_slope = is_negative(res_x_2 - res_x_1)

#FINDS Y_DOT DIRECTION
temp_y = vector(r2_dot_temp.x, r2_dot_temp.y + 0.01, r2_dot_temp.z)
res_y_1 = find_residuals(index2, Decs, RAs, JDs, r2_temp, r2_dot_temp)
res_y_2 = find_residuals(index2, Decs, RAs, JDs, r2_temp, temp_y)
y_dot_slope = is_negative(res_y_2 - res_y_1)

#FINDS Z_DOT DIRECTION
temp_z = vector(r2_dot_temp.x, r2_dot_temp.y, r2_dot_temp.z + 0.01)
res_z_1 = find_residuals(index2, Decs, RAs, JDs, r2_temp, r2_dot_temp)
res_z_2 = find_residuals(index2, Decs, RAs, JDs, r2_temp, temp_z)
z_dot_slope = is_negative(res_z_2 - res_z_1)

print "X_SLOPE:", x_slope
print "Y_SLOPE:", y_slope
print "Z_SLOPE:", z_slope
print "X_DOT_SLOPE:", x_dot_slope
print "Y_DOT_SLOPE:", y_dot_slope
print "Z_DOT_SLOPE:", z_dot_slope
lower_limit = 0.0001
it = 0


try:
    while it < 500:
        it += 1
        print "ITERATION: ", it
        print "R2", r2_temp
        print "R2_DOT", r2_dot_temp
        change_basic = 0.1
        #OPTIMIZES X
        last_res = find_residuals(index2, Decs, RAs, JDs, r2_temp, r2_dot_temp)
        change_interval = change_basic
        print last_res
        while True:
            new_vector = vector(r2_temp.x + change_interval * x_slope, r2_temp.y, r2_temp.z)
            new_res = find_residuals(index2, Decs, RAs, JDs, new_vector, r2_dot_temp)
            if new_res >= last_res:
                if change_interval <= lower_limit:
                    last_res = new_res
                    break
                change_interval *= 0.85
                x_slope *= -1.0
            else:
                last_res = new_res
                r2_temp = new_vector
            print last_res
        print "X DONE"
        #OPTIMIZES Y
        last_res = find_residuals(index2, Decs, RAs, JDs, r2_temp, r2_dot_temp)
        change_interval = change_basic
        print last_res
        while True:
            new_vector = vector(r2_temp.x, r2_temp.y + change_interval * y_slope, r2_temp.z)
            new_res = find_residuals(index2, Decs, RAs, JDs, new_vector, r2_dot_temp)
            if new_res >= last_res:
                if change_interval <= lower_limit:
                    last_res = new_res
                    break
                change_interval *= 0.85
                y_slope *= -1.0
            else:
                last_res = new_res
                r2_temp = new_vector
            print last_res
        print "Y DONE"
        #OPTIMIZES Z
        last_res = find_residuals(index2, Decs, RAs, JDs, r2_temp, r2_dot_temp)
        change_interval = change_basic
        print last_res
        while True:
            new_vector = vector(r2_temp.x, r2_temp.y, r2_temp.z + change_interval * z_slope)
            new_res = find_residuals(index2, Decs, RAs, JDs, new_vector, r2_dot_temp)
            if new_res >= last_res:
                if change_interval <= lower_limit:
                    last_res = new_res
                    break
                change_interval *= 0.85
                z_slope *= -1.0
            else:
                last_res = new_res
                r2_temp = new_vector
            print last_res
        print "Z DONE"
        #OPTIMIZES X_DOT
        last_res = find_residuals(index2, Decs, RAs, JDs, r2_temp, r2_dot_temp)
        change_interval = change_basic
        print last_res
        while True:
            new_vector = vector(r2_dot_temp.x + change_interval * x_dot_slope, r2_dot_temp.y, r2_dot_temp.z)
            new_res = find_residuals(index2, Decs, RAs, JDs, r2_temp, new_vector)
            if new_res >= last_res:
                if change_interval <= lower_limit:
                    last_res = new_res
                    break
                change_interval *= 0.85
                x_dot_slope *= -1.0
            else:
                last_res = new_res
                r2_dot_temp = new_vector
            print last_res
        print "X DOT DONE"
        #OPTIMIZES Y_DOT
        last_res = find_residuals(index2, Decs, RAs, JDs, r2_temp, r2_dot_temp)
        change_interval = change_basic
        print last_res
        while True:
            new_vector = vector(r2_dot_temp.x, r2_dot_temp.y + change_interval * y_dot_slope, r2_dot_temp.z)
            new_res = find_residuals(index2, Decs, RAs, JDs, r2_temp, new_vector)
            if new_res >= last_res:
                if change_interval <= lower_limit:
                    last_res = new_res
                    break
                change_interval *= 0.85
                y_dot_slope *= -1.0
            else:
                last_res = new_res
                r2_dot_temp = new_vector
            print last_res
        print "Y DOT DONE"
        #OPTIMIZES Z_DOT
        last_res = find_residuals(index2, Decs, RAs, JDs, r2_temp, r2_dot_temp)
        change_interval = change_basic
        print last_res
        while True:
            new_vector = vector(r2_dot_temp.x, r2_dot_temp.y, r2_dot_temp.z + change_interval * z_dot_slope)
            new_res = find_residuals(index2, Decs, RAs, JDs, r2_temp, new_vector)
            if new_res >= last_res:
                if change_interval <= lower_limit:
                    last_res = new_res
                    break
                change_interval *= 0.85
                z_dot_slope *= -1.0
            else:
                last_res = new_res
                r2_dot_temp = new_vector
            print last_res
        print "Z DOT DONE"
        change_basic /= 100.0
        lower_limit /= 100.0
finally:
    print "ITERATION:", it
    print "R2", r2_temp
    print "R2_DOT", r2_dot_temp
    print "FINAL RES:", last_res