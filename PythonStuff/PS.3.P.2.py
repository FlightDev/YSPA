from math import *

print "hi"

def to_unit_vector(RA, declination):
    s = radians(RA / 24.0 * 360)
    h = radians(declination)
    z = sin(h)
    y = cos(h) * sin(s)
    x = cos(h) * cos(s)
    return x, y, z

def rotate_vector(vector):
    x = vector[0]
    radius = sqrt((vector[1] ** 2) + (vector[2] ** 2))
    angle = arctan(vector[2] / vector[1])
    angle -= radians(23.5)
    y = cos(angle) * radius
    z = sin(angle) * radius
    return x, y, z

def into

print to_unit_vector(6, 0)
print to_unit_vector(0, 0)
print to_unit_vector(6, 180)
print rotate_vector(to_unit_vector(6, 10))
