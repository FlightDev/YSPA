import de421
from jplephem import Ephemeris
from visual import *
from  math import *

AU = 149597870.7 # km/AU
jd = 2458315.75 # July 1, 2018, 0:00 EDT
epsilon = 23.43687   # obliquity of the Ecliptic
eph = Ephemeris(de421)

def location_to_angles(location):
    dec = asin(location.z)
    ra = acos(location.x/cos(dec))
    if location.y < 0:
        ra = 2 * pi - ra
    return degrees(ra), degrees(dec)

barycenter = eph.position('earthmoon', jd)
moonvector = eph.position('moon', jd)
earth = barycenter - moonvector * eph.earth_share
R0 = vector(earth)/AU # This is the Sun to Earth center vector in Equatorial system

R_geocentric = 1*R0
print R_geocentric
sun_to_asteroid = vector(0.751289875, -1.0397039, 0.08361484)
#sun_to_asteroid = vector(0.244, 2.17, -0.445)
sun_to_asteroid = sun_to_asteroid.rotate(radians(epsilon), vector(1, 0, 0))
print sun_to_asteroid
earth_to_asteroid = norm(-R_geocentric + sun_to_asteroid) # Gets unit vector from Earth to Asteroid
print earth_to_asteroid
angles = location_to_angles(earth_to_asteroid)


print angles[0], angles[1]
