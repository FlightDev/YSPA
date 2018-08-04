from visual import *
import math as m

scene = display(forward=(-1, -1, 1))

def is_float(value): #determines if a given string can be cast to a float
    try:
        float(value)
        return True
    except ValueError:
        return False
RA = 0
dec = 0
prose = label(text = "[RA]")
while True:
    if scene.kb.keys: # event waiting to be processed?
        s = scene.kb.getkey() # get keyboard info
        if s == '\n' and is_float(prose.text):
            RA = float(prose.text)
            break
        elif len(s) == 1:
            prose.text += s # append new character
        elif ((s == 'backspace' or s == 'delete') and
                len(prose.text)) > 0:
            prose.text = prose.text[:-1] # erase letter
        elif s == 'shift+delete':
            prose.text = '' # erase all text
prose.visible = False
prose_one = label(text = "[Declination]")
while True:
    if scene.kb.keys: # event waiting to be processed?
        s = scene.kb.getkey() # get keyboard info
        if s == '\n' and is_float(prose_one.text):
            dec = float(prose_one.text)
            break
        elif len(s) == 1:
            prose_one.text += s # append new character
        elif ((s == 'backspace' or s == 'delete') and
                len(prose_one.text)) > 0:
            prose_one.text = prose_one.text[:-1] # erase letter
        elif s == 'shift+delete':
            prose_one.text = '' # erase all text
prose_one.visible = False

earth = sphere(pos = (0, 0, 0), material=materials.earth, radius = 3, opacity = 0.5)
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

while True:
    rate(10)
    earth.rotate(angle = m.radians(3), axis = vector(0, 1, 0))
