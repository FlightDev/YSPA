from visual import *
from math import *
import de421
from jplephem import Ephemeris


def accel(pos, mu):
    return -1.0 * mu * pos / (mag(pos) ** 3)

def verlet_position(current_pos, last_pos, d_tau, mu):
    r1 = 2.0 * current_pos - last_pos + accel(current_pos, mu) * (d_tau ** 2)
    return [r1, current_pos]

tau = 0.0
t = 0.00
end_t = 5.00
d_tau = 0.000001
k = 0.01720209895
end_tau = k * (end_t - t)
mu = 1.0

r_last = vector(0.244, 2.170, -0.445)
r_dot = vector(-0.731, -0.0041, 0.0502)
r_now = r_last + r_dot * d_tau + 0.5 * accel(r_last, mu) * (d_tau ** 2)
tau += d_tau

#sun = sphere(color = color.yellow, pos = (0, 0, 0), radius = 0.3)
#o2 = sphere(color = color.blue, pos = (r_now.x, r_now.y, r_now.z), radius = 0.1)
#o2_t = curve(color = color.blue)

while tau <= end_tau:
    #rate(50)
    #print tau, end_tau
    tau += d_tau
    temp = verlet_position(r_now, r_last, d_tau, mu)
    r_last = temp[1]
    r_now = temp[0]
    #o2.pos = (r_now.x, r_now.y, r_now.z)
    #o2_t.append((r_now.x, r_now.y, r_now.z))

print r_now


AU = 149597870.7 # km/AU
jd = 2458305.666667 # July 1, 2018, 0:00 EDT
epsilon = 23.43687   # obliquity of the Ecliptic
eph = Ephemeris(de421)


def location_to_angles(location):
    dec = asin(location.z)
    ra = acos(location.x/cos(dec))
    return degrees(ra), degrees(dec)

barycenter = eph.position('earthmoon', jd)
moonvector = eph.position('moon', jd)
earth = barycenter - moonvector * eph.earth_share
R0 = vector(earth)/AU # This is the Sun to Earth center vector in Equatorial system

R_geocentric = 1*R0
print R_geocentric
sun_to_asteroid = r_now
sun_to_asteroid = sun_to_asteroid.rotate(radians(epsilon), vector(1, 0, 0))
print sun_to_asteroid
earth_to_asteroid = norm(-R_geocentric + sun_to_asteroid) # Gets unit vector from Earth to Asteroid
print earth_to_asteroid
angles = location_to_angles(earth_to_asteroid)


print angles[0], angles[1]


print 12345
