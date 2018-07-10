from visual import *
import math as m

celest = sphere(pos = (0, 0, 0), opacity = 0.1, color = color.white, radius = 10)
point_of_aries = minute_hand = arrow(pos = (0, 0, 0), axis = vector(10, 0, 0), shaftwidth = 0.4, color = color.white)

equator = ring(pos = (0, 0, 0), axis = (0, 1, 0), thickness = 0.1, radius = 10, color = color.green)
eq_label = label(pos = (0, 0, 10), text = "Celestial Equator", xoffset = 15, yoffset = 15, height=10, font = 'sans')
ecliptic = ring(pos = (0, 0, 0), axis = (0, m.cos(radians(23.5)), m.sin(radians(23.5))), thickness = 0.1, radius = 10, color = color.magenta)
ecliptic_label = label(pos = (0, 0, 10), text = "Celestial Equator", xoffset = 15, yoffset = 15, height=10, font = 'sans')
prime_meridian = ring(pos = (0, 0, 0), axis = (0, 0, 1), thickness = 0.1, radius = 10, color = color.green)
prime_label = label(pos = (0, 0, 10), text = "Celestial Equator", xoffset = 15, yoffset = 15, height=10, font = 'sans')

north_pole = box(pos = (0, 0, 0), height = 25, width = 0.1, length = 0.1)
