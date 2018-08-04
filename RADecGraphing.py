import matplotlib.pyplot as plt
from visual import *
from math import *
import de421
from jplephem import Ephemeris

def get_earth_to_sun(JD):
    eph = Ephemeris(de421)
    AU = 149597870.7 # km/AU
    barycenter = eph.position('earthmoon', JD)
    moonvector = eph.position('moon', JD)
    earth = barycenter - moonvector * eph.earth_share
    return -1.0 * vector(earth)/AU

def accel(pos, mu):
    return -1.0 * mu * pos / (mag(pos) ** 3)

def verlet_position(current_pos, last_pos, d_tau, mu):
    r1 = 2.0 * current_pos - last_pos + accel(current_pos, mu) * (d_tau ** 2)
    return [r1, current_pos]
"""
data = open('radecdata.txt', 'r')
s = ""
for line in data:
    s += line
lines = s.split('\n')
lines = lines[:len(lines) - 1]
data = []
for line in lines:
    data.append(line.split(' '))
jd = []
ra = []
dec = []
for i in range(len(data)):
    jd.append(data[i][0])
    ra.append(data[i][1])
    dec.append(data[i][2])
    <0.675905, -1.54554, 0.170722>
    VELOCITY:  <0.0580764, 0.0165273, 0.307812>
"""

k = 0.01720209895

def neg(num):
    if num < 0:
        return -1.0
    else:
        return 1.0

def find_final_vector(jd):
    base_day = 2458316.711586
    r_last = vector(0.675905, -1.54554, 0.170722)
    r_dot = vector(0.0580764, 0.0165273, 0.307812)
    tau = 0.0
    t = 0.00
    if jd == base_day:
        return r_last
    end_t = jd - base_day
    d_tau = 0.001 * neg(end_t)
    k = 0.01720209895
    end_tau = k * (end_t - t)
    mu = 1.0
    r_dot /= k
    r_now = r_last + r_dot * d_tau + 0.5 * accel(r_last, mu) * (d_tau ** 2)
    tau += d_tau

    while abs(tau) <= abs(end_tau):
        tau += d_tau
        temp = verlet_position(r_now, r_last, d_tau, mu)
        r_last = temp[1]
        r_now = temp[0]

    return r_now

AU = 149597870.7 # km/AU
epsilon = 23.43687   # obliquity of the Ecliptic
eph = Ephemeris(de421)

def location_to_angles(location):
    dec = asin(location.z)
    #print dec
    ra = acos(location.x/cos(dec))
    if location.y < 0:
        ra = 2 * pi - ra
    return degrees(ra), degrees(dec)

def get_ra_dec(r, julian_date):
    barycenter = eph.position('earthmoon', julian_date)
    moonvector = eph.position('moon', julian_date)
    earth = barycenter - moonvector * eph.earth_share
    R0 = vector(earth)/AU # This is the Sun to Earth center vector in Equatorial system

    R_geocentric = -1*R0
    sun_to_asteroid = r_now
    sun_to_asteroid = sun_to_asteroid.rotate(radians(epsilon), vector(1, 0, 0))
    earth_to_asteroid = norm(R_geocentric + sun_to_asteroid)
    #earth_to_asteroid = vector(0.938685, -0.326344, 0.111221)
    angles = location_to_angles(earth_to_asteroid)
    return angles[0], angles[1]

def find_residuals(check_index, decs, ras, JDs, r, r_dot):
    dec = []
    ra = []

    for i in range(len(decs)):
        k = 0.01720209895
        time_difference = JDs[i] - JDs[check_index]
        tau = 0.0
        d_tau = 0.0001 * neg(time_difference)
        end_tau = time_difference * k
        mu = 1.0
        r_last = r
        r_now = r_last + r_dot * d_tau + 0.5 * accel(r_last, mu) * (d_tau ** 2)
        tau += d_tau

        while abs(tau) <= abs(end_tau):
            tau += d_tau
            temp = verlet_position(r_now, r_last, d_tau, mu)
            r_last = temp[1]
            r_now = temp[0]


        AU = 149597870.7# km/AU
        jd = JDs[i]
        eph = Ephemeris(de421)

        R_geocentric = get_earth_to_sun(jd)
        sun_to_asteroid = r_now
        earth_to_asteroid = norm(R_geocentric + sun_to_asteroid) # Gets unit vector from Earth to Asteroid
        RA, Dec = location_to_angles(earth_to_asteroid)
        ra.append(RA)
        dec.append(Dec)
        tau = 0
        d_tau = 0.000001
        k = 0.01720209895
        mu = 1.0

    return ra, dec

jd = [2458281.856678, 2458281.857743, 2458281.858808, 2458281.859606, 2458281.863056, 2458286.836458, 2458315.564236, 2458315.569248, 2458315.578218, 2458316.711586, 2458316.718773, 2458316.725613, 2458316.727662, 2458326.554572, 2458326.561562, 2458326.568299, 2458326.575336, 2458326.582002, 2458326.590035]
ra = [297.31006665, 297.30940004999997, 297.3086916, 297.3077958, 297.30546255, 293.6085291, 272.429166, 272.42109585, 272.4214917, 271.8, 271.7958333, 271.7919291, 271.7908584, 267.49022505, 267.4877916, 267.4852458, 267.48287505, 267.4872495, 267.477525]
dec = [28.625444, 28.62598333, 28.62657222, 28.62676111, 28.62928056, 31.40670556, 35.9272222, 35.92431944, 35.92481944, 35.7675, 35.7666, 35.76640278, 35.76515278, 33.76921667, 33.76738056, 33.76566111, 33.76386944, 33.76268056, 33.760272222]
ras = []
decs = []
"""
<0.59378, -1.68061, 0.283395>
VELOCITY:  <-0.489096, -0.451879, 0.411871>

"""
#r = vector(0.604909, -1.59533, 0.23126)
#rdot = vector(-0.513917, -0.241117, 0.273125)
r = vector(0.59378, -1.68061, 0.283395)
rdot = vector(-0.489096, -0.451879, 0.411871)
#r = r.rotate(-epsilon, vector(1, 0, 0))
#r = rdot.rotate(-epsilon, vector(1, 0, 0))
ras, decs = find_residuals(9, dec, ra, jd, r, rdot)

f ,(ax1, ax2) = plt.subplots(2)

ax1.plot(jd, ras)
ax2.plot(jd, decs)



ax1.scatter(jd, ra, color='r')
ax2.scatter(jd, dec, color='b')
ax1.set_xlabel("Julian Date")
ax2.set_xlabel("Julian Date")
ax1.set_ylabel("Right Ascension")
ax2.set_ylabel("Declination")
plt.show()
