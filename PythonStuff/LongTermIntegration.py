import rebound
from math import *
from visual import *
import matplotlib.pyplot as plt

k = 0.01720209895

sim = rebound.Simulation()



#SUN
sim.add(m = 1.0)
#ASTEROID
sim.add(m = 0.0, x = 0.0, y =0.0 , z = 1.5, vx = 0.7, vy = 0.7, vz = 0.0)
#EARTH
sim.add(m = 1.0/330000.0, x = 0.60507, y = -0.80851, z = -0.00005299, vx = 0.013529/k, vy = 0.010179/k, vz = 0.0000002635/k)
#JUPITER
sim.add(m = 1.0/1000.0, x = -3.13617, y = -4.37686, z = 0.088304, vx = 0.00604369/k, vy = -0.0040367/k, vz = -0.00011849/k)
#SATURN
sim.add(m = 1.0/3300.0, x = 1.1523, y = -9.9908, z = 0.127841, vx = 0.0052346/k, vy = 0.00062367/k, vz = -0.00021891/k)

sim.move_to_com()

ecc = []
semimaj = []
times = []
I = []

t0 = 2458329.50
time_step = 0.01000
t = 0

ps = sim.particles

earth = vector(ps[2].x, ps[2].y, ps[2].z)
asteroid = vector(ps[1].x, ps[1].y, ps[1].z)
min_dist = 1.06056208892

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
