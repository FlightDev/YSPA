from math import *

def to_unit_vector(RA, dec):
    RA = radians(RA * 15)
    dec = radians(dec)
    z = sin(dec)
    x = cos(RA) * cos(dec)
    y = sin(RA) * cos(dec)
    return x, y, z

def rotate_vector(RA, dec):
    epsilon = radians(23.5)
    dec = radians(dec)
    RA = radians(RA * 15)
    x = cos(RA) * cos(dec)
    y = cos(epsilon) * sin(RA) * cos(dec) + sin(epsilon) * sin (dec)
    z = cos(epsilon) * sin(dec) - sin(epsilon) * sin(RA) * cos(dec)
    return x, y, z

def into_ecliptical(RA, dec):
    dec = radians(dec)
    RA = radians(RA * 15)
    epsilon = radians(23.5)
    lat = asin(cos(epsilon) * sin(dec) - sin(epsilon) * sin(RA) * cos(dec))
    long = acos((cos(RA) * cos(dec))/cos(lat))
    long = degrees(long)
    lat = degrees(lat)
    return decimal_to_sexagesimal(long), decimal_to_sexagesimal(lat)


def decimal_to_sexagesimal(angle):
    time = []
    negative = 1.0
    sexagesimal = ""
    if angle < 0:
        negative = -1.0
        angle *= -1.0
    time.append(int(angle)) #finds the degree of the angle
    angle -= int(angle)
    time.append(int(angle*60.0)) #find the arcminutes of the angle
    angle -= time[1] /  60.0
    time.append(angle*3600.0) #find the arcseconds of the angle
    time[0] *= negative
    time[0] = int(time[0])
    return time #returns the angle in the correct format

print into_ecliptical(18.28, 36.255)
print into_ecliptical(0, 0)
