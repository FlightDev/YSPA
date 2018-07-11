from visual import *
import math as m
from Tkinter import *
"""
master = Tk()

Label(master, text="Right Ascension(Hour Angle)").grid(row=0)
Label(master, text="Declination (Degrees)").grid(row=1)

e1 = Entry(master)
e2 = Entry(master)

e1.grid(row=0, column=1)
e2.grid(row=1, column=1)

print e1
"""
RA = 0
dec = 0

earth = sphere(pos = (0, 0, 0), material=materials.earth, radius = 3)
earth.rotate(angle = m.pi/4, axis = vector(0, 1, 0))
celest = sphere(pos = (0, 0, 0), opacity = 0.1, color = color.white, radius = 10)
#point_of_aries = arrow(pos = (0, 0, 0), axis = vector(10, 0, 0), shaftwidth = 0.4, color = color.white)

equator = ring(pos = (0, 0, 0), axis = (0, 1, 0), thickness = 0.1, radius = 10, color = color.green)
eq_label = label(pos = (0, 0, 10), text = "Celestial Equator", xoffset = 15, yoffset = 15, zoffset = 15, height=10, font = 'sans')
ecliptic = ring(pos = (0, 0, 0), axis = (0, m.cos(radians(23.5)), m.sin(radians(23.5))), thickness = 0.1, radius = 10, color = color.magenta)
ecliptic_label = label(pos = (0, m.sin(radians(23.5)) * 10, -1 * m.cos(radians(23.5)) * 10), text = "Ecliptic", xoffset = 15, yoffset = 15, height=10, font = 'sans')
prime_meridian = ring(pos = (0, 0, 0), axis = (0, 0, 1), thickness = 0.1, radius = 10, color = color.green)
prime_label = label(pos = (m.cos(radians(45)) * 10, m.cos(radians(45)) * 10, 0), text = "Prime Meridian", xoffset = 15, yoffset = 15, height=10, font = 'sans')


north_pole = box(pos = (0, 0, 0), height = 25, width = 0.1, length = 0.1)
pole_label = label(pos = (0, 7, 0), text = "North Celestial Pole", xoffset = 15, yoffset = 15, height=10, font = 'sans')

def add_arrow(RA, dec):
    RA = (RA) % 24
    RA = m.radians(RA * 15)
    dec = m.radians(dec)
    print RA, dec
    y = m.sin(dec) * 10
    x = m.cos(RA) * m.cos(dec) * 10
    z = m.sin(RA) * m.cos(dec) * -10
    print x, y, z
    point = arrow(pos = (0, 0, 0), axis = vector(x, y, z), shaftwidth = 0.1, color = color.white)

add_arrow(RA, dec)

#mainloop()
