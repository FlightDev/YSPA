import rebound
import math
from visual import *
import numpy as np
import matplotlib.pyplot as plt

k  =  0.01720209895
t0 = 2458329.50

sim = rebound.Simulation()

sim.add(m = 1.0) # Sun
# Earth
sim.add(m = 0.000003003, x = 6.050702901916951E-01, y = -8.085113449604454E-01, z = -5.299403058075317E-05,
vx = 1.352973714877966E-02/k, vy = 1.017946114599288E-02/k, vz = 2.635007516883264E-07/k)
# Jupiter
sim.add(m = 0.0009543, x = -3.136171885349048, y = -4.376867815031186, z = 8.830406151426499E-02,
vx = 6.043693520370460E-03/k, vy = -4.036729465948762E-03/k, vz = -1.184985444663097E-04/k)
# Saturn
sim.add(m = 0.0002857, x = 1.152369370228355, y = -9.990801814597656, z = 1.278418069572467E-01,
vx = 5.234609262074102E-03/k, vy = 6.209478126028669E-04/k, vz = -2.189111161287214E-04/k)
# MFFA
sim.add(m = 0, x = 0, y = 0, z = 1.5, vx = 0.7, vy = 0.7, vz = 0)


sim.move_to_com()
ps = sim.particles
time = 0
end_time = 2. * math.pi * 50

e_values = []
a_values = []
i_values = []
times = []

while time < end_time:
    # Earth
    sim.integrate(time)

    rdot = vector(ps[3].vx, ps[3].vy, ps[3].vz)
    r = vector(ps[3].x, ps[3].y, ps[3].z)
    h = vector(cross(r, rdot))
    h_mag = mag(h)
    rhat = norm(r)
    e = mag(vector(cross(rdot, h) - rhat))
    a = ((h_mag)**2)/(1-e**2)
    z_hat = vector(1,0,0)
    i = acos(dot(norm(h), z_hat))
    e_values.append(e)
    a_values.append(a)
    i_values.append(i)
    times.append(time)
    print time
    time += 0.01


f,(ax1, ax2, ax3) = plt.subplots(3)
ax1.plot(times, i_values)
ax2.plot(times, e_values)
ax3.plot(times, a_values)
ax1.set_xlabel("Julian Day")
ax2.set_xlabel("Julian Day")
ax3.set_xlabel("Julian Day")
ax1.set_ylabel("Inclination")
ax2.set_ylabel("Eccentricity")
ax3.set_ylabel("Semi-major Axis")
plt.show()

"""
f,(ax1, ax2) = plt.subplots(2)
ax1.plot(times, e_values)
ax2.plot(times, a_values)
ax1.set_xlabel("Julian Day")
ax1.set_ylabel("Eccentricity")
ax2.set_xlabel("Julian Day")
ax2.set_ylabel("Semi-major Axis")
plt.show()
"""
