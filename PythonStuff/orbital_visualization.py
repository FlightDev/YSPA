from visual import *
from math import *

x = arrow(pos = (0, 0, 0), axis = vector(10, 0, 0), shaftwidth = 0.1, color = color.red)
y = arrow(pos = (0, 0, 0), axis = vector(0, 10, 0), shaftwidth = 0.1, color = color.green)
z = arrow(pos = (0, 0, 0), axis = vector(0, 0, 10), shaftwidth = 0.1, color = color.blue)

def graph_orbital(q, e):
    a = q / (1.0 - e)
    return [a, q, e]

def rotate_vector(axis, angle, v):
    x_temp = v[0]
    y_temp = v[1]
    z_temp = v[2]
    angle = -1.0 * radians(angle)
    print angle
    if axis == 'z':
        x = x_temp * cos(angle) - y_temp * sin(angle)
        y = x_temp * sin(angle) + y_temp * cos(angle)
        z = z_temp
        print 1
    elif axis == 'y':
        y = y_temp * cos(angle) - z_temp * sin(angle)
        z = z_temp * sin(angle) + z_temp * cos(angle)
        x = x_temp
        print 2
    elif axis == 'x':
        z = z_temp * cos(angle) - x_temp * sin(angle)
        x = z_temp * sin(angle) + x_temp * cos(angle)
        y = y_temp
        print 3
    else:
        print 4
        return -1
    print (x, y, z)
    return (x, y, z)

q = [1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
e = [0.0, 0.2, 0.4, 0.6, 0.8, 0.9999]
Omega = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
omega = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
inclination = [45.0, 30.0, 15.0, 5.0, 0.0, -10.0]
colors = [color.red, color.orange, color.yellow, color.green, color.blue, color.magenta]
orbits = []

#for i in range(1):

sun = sphere(radius = 0.1, material = materials.emissive, color = color.yellow)

for i in range(len(q)):
    orbits.append(graph_orbital(q[i], e[i]))
    orbit = curve(color = colors[i])
    for j in range(257):
        angle = j / 128. * pi
        r = orbits[i][0] * (1 - orbits[i][2] ** 2) / (1 + orbits[i][2] * cos(angle))
        v = (cos(angle) * r, sin(angle) * r, 0.0)
        print v
        v = rotate_vector('z', omega[i], v)
        v = rotate_vector('x', inclination[i], v)
        v = rotate_vector('z', Omega[i], v)
        orbit.append(pos=(v[0] * -1.0, v[1], v[2]))
