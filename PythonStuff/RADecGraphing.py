import matplotlib.pyplot as plt
from visual import *
from math import *
import de421
from jplephem import Ephemeris


def accel(pos, mu):
    return -1.0 * mu * pos / (mag(pos) ** 3)

def verlet_position(current_pos, last_pos, d_tau, mu):
    r1 = 2.0 * current_pos - last_pos + accel(current_pos, mu) * (d_tau ** 2)
    return [r1, current_pos]

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



k = 0.01720209895

def neg(num):
    if num < 0:
        return -1.0
    else:
        return 1.0

def find_final_vector(jd):
    base_day = 2458315.7500000
    r_last = vector(0.751289875,  -1.0397039,  0.083614846)
    r_dot = vector(1.058954097E-02, 1.4886272E-02, -5.53132485E-04)
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

    if end_t > 0:
        while tau <= end_tau:
            tau += d_tau
            temp = verlet_position(r_now, r_last, d_tau, mu)
            r_last = temp[1]
            r_now = temp[0]
    else:
        while tau >= end_tau:
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


j_d = 2458315.750
time = [-3, 3]
times = []
for i in range(len(time)):
    times.append(time[i] + j_d)
ras = []
decs = []

for t in time:
    r_now = find_final_vector(j_d + t)
    angles = get_ra_dec(r_now, j_d + t)
    ras.append(angles[0])
    decs.append(angles[1])
f, (ax1, ax2) = plt.subplots(2, sharex = False, sharey = False)

ax1.plot(times, ras)
ax2.plot(times, decs)



ax1.scatter(jd, ra, color='r')
ax2.scatter(jd, dec, color='b')
ax1.set_xlabel("Julian Date")
ax2.set_xlabel("Julian Date")
ax1.set_ylabel("Right Ascension")
ax2.set_ylabel("Declination")
plt.show()
