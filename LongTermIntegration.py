import rebound
from math import *
from visual import *
import matplotlib.pyplot as plt

k = 0.01720209895

sim = rebound.Simulation()
"""
0.59574, -1.70604, 0.30176
-0.477734, -0.509794, 0.450681
"""

r = vector(0.511, -1.5195, -0.5674)
r_dot = vector(0.4615, 0.2939, -0.2884)
epsilon = 23.43687   # obliquity of the Ecliptic
#r = r.rotate(epsilon, vector(1, 0, 0))
#r_dot = r_dot.rotate(epsilon, vector(1, 0, 0))
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

ecc = []
semimaj = []
times = []
I = []

t0 = 2458315.578218
time_step = 0.01000
t = 0

ps = sim.particles

#earth = vector(ps[2].x, ps[2].y, ps[2].z)
#asteroid = vector(ps[1].x, ps[1].y, ps[1].z)
#min_dist = 1.06056208892

sun = sphere(pos = (ps[0].x, ps[0].y, ps[0].z), radius = 0.1, color = color.yellow)
earth = sphere(pos = (ps[2].x, ps[2].y, ps[2].z), radius = 0.02, material = materials.earth)
jupiter = sphere(pos = (ps[3].x, ps[3].y, ps[3].z), radius = 0.04, color = color.orange)
saturn = sphere(pos = (ps[4].x, ps[4].y, ps[4].z), radius = 0.03, color = (255, 200, 200))
asteroid = sphere(pos = (ps[1].x, ps[1].y, ps[1].z), radius = 0.005, color = (255, 0, 0))
sun_trail = curve(color = color.yellow)
earth_trail = curve(color = color.blue)
jupiter_trail = curve(color = color.orange)
saturn_trail = curve(color = (255, 200, 200))
asteroid_trail = curve(color = color.red)

min_dist = 10000000
while t < 25000/k:
    #rate(100)
    sim.integrate(t)
    t += time_step
    earthpos = vector(ps[2].x, ps[2].y, ps[2].z)
    asteroidpos = vector(ps[1].x, ps[1].y, ps[1].z)
    dist = mag(earthpos - asteroidpos)
    if dist < min_dist:
        min_dist = dist
        print "NEW MIN DIST:", dist, " AT TIME:", t
    #print t
    sun.pos = (ps[0].x, ps[0].y, ps[0].z)
    earth.pos = (ps[2].x, ps[2].y, ps[2].z)
    jupiter.pos = (ps[3].x, ps[3].y, ps[3].z)
    saturn.pos = (ps[4].x, ps[4].y, ps[4].z)
    asteroid.pos = (ps[1].x, ps[1].y, ps[1].z)
    sun_trail.append(sun.pos)
    earth_trail.append(earth.pos)
    jupiter_trail.append(jupiter.pos)
    saturn_trail.append(saturn.pos)
    asteroid_trail.append(asteroid.pos)
print min_dist


"""
while t < 2 * pi * 100.0:
    sim.integrate(t)
    t += time_step
    rdot = vector(ps[1].vx, ps[1].vy, ps[1].vz)
    r = vector(ps[1].x, ps[1].y, ps[1].z)
    h = cross(r, rdot)
    e = mag(cross(rdot, h) - norm(r))
    i = acos(dot(norm(h), vector(0, 0, 1)))
    #asteroid = vector(ps[1].x, ps[1].y, ps[1].z)
    #current_dist = mag(r - asteroid)
    #if current_dist < min_dist:
    #    min_dist = current_dist
    a = mag(h) ** 2 / (1 - e ** 2)
    ecc.append(e)
    semimaj.append(a)
    I.append(i)
    times.append(t)
    print t
"""
f, (ax1, ax2, ax3) = plt.subplots(3, sharex = True, sharey = False)
plt.title("Orbital Elements vs. Time of Asteroid")
ax1.set_ylabel("Eccentricity")
ax2.set_ylabel("Semi-Major Axis")
ax3.set_ylabel("Inclination")
ax1.set_xlabel("Days after JD " + str(t0))
ax2.set_xlabel("Days after JD " + str(t0))
ax3.set_xlabel("Days after JD " + str(t0))
ax1.plot(times, ecc, color='g')
ax2.plot(times, semimaj, color='r')
ax3.plot(times, I, color='b')
plt.show()

"""
while t < 2 * pi:
    sim.integrate(t)
    earth = vector(ps[2].x, ps[2].y, ps[2].z)
    asteroid = vector(ps[1].x, ps[1].y, ps[1].z)
    current_dist = mag(earth - asteroid)
    if current_dist < min_dist:
        min_dist = current_dist

print min_dist
"""
