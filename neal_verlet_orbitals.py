from visual import *
from math import *
from random import *
"""
import de421
from jplephem import Ephemeris
"""

class orbital_object:
    def __init__(self, position, mass, d_tau, velocity, radius):
        self.r1 = position
        self.r2 = position
        self.m = mass
        self.a = -1
        self.vi = velocity
        self.rad = radius

    def update_position(self, d_tau, objects):
        self.update_acceleration(objects)
        r_temp = self.r1
        self.r1 = 2 * self.r1 - self.r2 + self.a * (d_tau ** 2)
        self.r2 = r_temp

    def update_acceleration(self, objects):
        net_accel = vector(0, 0, 0)
        for object in objects:
            rad = object.pos() - self.r1
            net_accel += object.get_mass()*rad/(mag(rad)**3)
        self.a = net_accel

    def pos(self):
        return self.r1

    def get_mass(self):
        return self.m

    def get_radius(self):
        return self.rad

    def set_initial(self, d_tau, objects):
        self.update_acceleration(objects)
        self.r1 = self.r2 + self.vi * d_tau + 0.5 * self.a * (d_tau ** 2)

def update_all_objects(d_tau, objects):
    if len(objects) < 2:
        return -1
    for i in range(len(objects)):
        temp_objects = []
        for j in range(len(objects)):
            if i != j:
                temp_objects.append(objects[j])
        objects[i].update_position(d_tau, temp_objects)

def initialize_objects(d_tau, objects):
    if len(objects) < 2:
        return -1
    for i in range(len(objects)):
        temp_objects = []
        for j in range(len(objects)):
            if i != j:
                temp_objects.append(objects[j])
        objects[i].set_initial(d_tau, temp_objects)

celest_obj = []
d_tau = 0.0001
sun = orbital_object(vector(0, 0, 0), 1, d_tau, vector(0, 0, 0), 0.3)
earth = orbital_object(vector(1, 0, 0), 1.0 / 330000, d_tau, vector(0, 1, 0), 0.05)
mars = orbital_object(vector(1.52, 0, 0), 1.0 / 3300000, d_tau, vector(0, 0.805, 0), 0.04)
celest_obj.append(sun)
celest_obj.append(earth)
celest_obj.append(mars)
print celest_obj[1].get_radius()
o = []
o1 = sphere(pos=celest_obj[0].pos(), radius = celest_obj[0].get_radius(), color = color.yellow, material = materials.emissive)
o2 = sphere(pos=celest_obj[1].pos(), radius = celest_obj[1].get_radius(), material = materials.earth)
o3 = sphere(pos=celest_obj[2].pos(), radius = celest_obj[2].get_radius(), color = color.red)
o.append(o1)
o.append(o2)
o.append(o3)

"""
for i in range(10):
    initial_position = vector(randint(500, 1751) / 1000.0, randint(500, 1751) / 1000.0, randint(500, 1751) / 1000.0)
    initial_velocity = vector(randint(100, 2501) / 1000.0, randint(100, 2501) / 1000.0, randint(100, 2501) / 1000.0)
    mass = 1.0 / randint(10000000, 1000000000)
    rad = 1.0 / randint(100, 1000)
    i_en = 2
    #o_temp = sphere(pos=celest_obj[i + i_en].pos(), radius = 0.3, color = color.green)
    o_temp = sphere(pos=celest_obj[i+i_en].pos(), radius = celest_obj[i+i_en].get_radius(), color = color.green)
    o.append(o_temp)
initialize_objects(d_tau, celest_obj)
"""

while True:
    rate(10000)
    update_all_objects(d_tau, celest_obj)
    for i in range(len(celest_obj)):
        o[i].pos = celest_obj[i].pos()
