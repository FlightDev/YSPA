from math import *

print "hi"

def to_unit_vector(RA, dec):
    RA = radians(RA * 15)
    dec = radians(dec)
    z = sin(dec)
    x = cos(RA) * cos(dec)
    y = sin(RA) * cos(dec)
    return x, y, z

def rotate_vector(vector):
    x = vector[0]
    radius = sqrt((vector[1] ** 2) + (vector[2] ** 2))
    angle = atan(vector[2] / vector[1])
    angle += radians(23.5)
    y = cos(angle) * radius
    z = sin(angle) * radius
    return x, y, z

def into_ecliptical(vector):

    return 1

print to_unit_vector(12, 90 )
print to_unit_vector(0, 0)
print to_unit_vector(6, 180)
