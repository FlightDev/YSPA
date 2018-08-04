import rebound
import math
from visual import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

k = 0.01720209895
how_long = 600
start_pos = vector(0, 0, 1.5)
start_v = vector(0.7, 0.7, 0)

sim = rebound.Simulation()

#sun:
sim.add(m = 1.)

#asteroid:
sim.add(m = 0, x = start_pos.x, y = start_pos.y, z = start_pos.z, vx = start_v.x/k, vy = start_v.y/k, vz = start_v.z/k)

#earth:
sim.add(m = 0.000003003, x = 6.050702901916951E-01, y = -8.085113449604454E-01, z = -5.299403058075317E-05, vx = (1.352973714877966E-02)/k, vy = (1.017946114599288E-02)/k, vz = (2.635007516883264E-07)/k )

#jupiter:
sim.add(m = 0.0009543, x = -3.136171264149830E+00, y = -4.376868856434548E+00, z = 8.830403590272071E-02, vx = 6.044270575179947E-03/k, vy =-4.035730426004453E-03/k, vz = -1.184535381952951E-04/k)

#saturn:
sim.add(m =.0002857,  x = 1.152370623788473E+00, y =-9.990803088412557E+00, z = 1.278423486688079E-01, vx = 5.235192499208867E-03/k, vy = 6.213724626462464E-04/k, vz = -2.191864499860967E-04/k )

sim.dt = 0.01
sim.move_to_com()

time = 0
end_time = 2.*math.pi*how_long

ps = sim.particles

#sim.integrate(end_time)
#earth position vector
r = vector (ps[1].x, ps[1].y, ps[1].z )
earth = vector (ps[2].x, ps[2].y, ps[2].z )
e_a_distance = mag(r) - mag(earth)
closest_distance = abs( mag(r) - mag(earth) )
closest_time = time#/k

while time < end_time:
    sim.integrate(time)
    r = vector(ps[1].x, ps[1].y, ps[1].z )
    #see if asteroid hits the earth
    earth = vector (ps[2].x, ps[2].y, ps[2].z )
    e_a_distance = mag(r) - mag(earth)
    if abs(e_a_distance) < closest_distance:
        closest_distance = mag(e_a_distance)
        closest_time = time#/k

    rdot = vector(ps[1].vx, ps[1].vy, ps[1].vz )

    time = time + 0.01

print "closest distance = ", closest_distance
print "time = ", closest_time
