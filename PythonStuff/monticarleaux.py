from math import *
import matplotlib.pyplot as plt
import numpy as np
#from visual import *
import de421
from jplephem import Ephemeris
import rebound

r = vector(0.59574, -1.70604, 0.30176)
r_dot = vector(-0.477734, -0.509794, 0.450681)
epsilon = 23.43687   # obliquity of the Ecliptic
r = r.rotate(-epsilon, vector(1, 0, 0))
rdot = r_dot.rotate(-epsilon, vector(1, 0, 0))
rlist = []
rdotlist = []

for i in range(1000):
    rlist.append(r+np.random.normal(0,0.001, 3))
    rdotlist.append(rdot + np.random.normal(0,0.001, 3))

sim = rebound.Simulation()
k  =  0.01720209895
JD =  2458329.50
sim.add(m = 1.)
sim.add(m = 0.000003003, x = 4.652709809216394E-01, y =-8.966566096285266E-01, z =-5.459720950324835E-05, vx= (1.500916294305504E-02)/k, vy= (7.813599276848059E-03)/k, vz= (-2.948047283284480E-07)/k)
sim.add(m=0, x=0.89917,  y=-1.56911,  z=-0.490915,  vx=0.507158,  vy=-0.0278803,  vz=0.463388)
sim.add(m = 0.0009543, x = -3.194977267469506, y = -4.337020438246942, z = 8.945405865419359E-02, vx= (5.986341370125854E-03)/k, vy=(-4.115361763361076E-03)/k, vz=(-1.167705200709614E-04)/k)
sim.add(m = 0.0002857166,x = 1.101173501013428, y =-9.996740706242282, z = 1.299843189539670E-01, vx = ( 5.239197851566815E-03)/k, vy= (5.929373084668421E-04)/k, vz=(-2.188911200273129E-04)/k)
for i in range(10000):
    sim.add(m=0, x=rlist[i].x,y=rlist[i].y,z=rlist[i].z,vx=rdotlist[i].x,vy=rdotlist[i].y,vz=rdotlist[i].z)
sim.status()
sim.dt=0.1
sim.move_to_com()
time=0.0
end_time=2*pi*1000.
ps = sim.particles
a = []
for i in range(10000):
    a.append(10000)
print a , len(a)
while(time<end_time):
    sim.integrate(time)
    earth = vector(ps[1].x,ps[1].y,ps[1].z)

    for i in range(10000):
        rtest = vector(ps[5+i].x,ps[5+i].y,ps[5+i].z)
        rho = mag(rtest - earth)
        if rho<a[i]:
            a[i]=rho
    time = time + sim.dt
print a
count = 0
for i in range(len(a)):
    if a[i]<0.001:
        count+=1
print count
