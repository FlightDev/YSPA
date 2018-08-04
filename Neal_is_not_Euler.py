from visual import *
from math import *


# Approximate masses of the earth and the sun
earth_mass = 1/330000
mars_mass = earth_mass/10


earth_position = vector(1.0, 0.0, 0.0)
earth_velocity = vector(0.0, 1.0, 0.0)
mars_position = vector(1.524, 0.0, 0.0)
mars_velocity = vector(0.0, 0.805, 0.0)

tau = 0.0
end_tau = 5.0
timestep = 0.0001
k = 0.01720209895

#scene.stereo = "redcyan"
sun = sphere(color = color.yellow, pos = (0, 0, 0), radius = 0.1)
earth = sphere(material = materials.earth, pos=earth_position, radius = 0.03)
mars = sphere(color = color.red, pos=mars_position, radius = 0.03)
earth_trail = curve(color = color.white)
mars_trail = curve(color = color.red)
x = arrow(pos = (0, 0, 0), axis = (0.2, 0, 0), color = color.red, shaftwidth = 0.01)
y = arrow(pos = (0, 0, 0), axis = (0, 0.2, 0), color = color.green, shaftwidth = 0.01)
z = arrow(pos = (0, 0, 0), axis = (0, 0, 0.2), color = color.blue, shaftwidth = 0.01)
parahelelion = arrow(pos = (0, 0, 0), axis = earth_position, color = color.white, shaftwidth = 0.01)
h = cross(earth_position, earth_velocity)
h_arrow = arrow(pos = (0, 0, 0), axis = h, color = color.white, shaftwidth = 0.001)

def find_period(pos):
    return sqrt((4 * pi ** 2 * mag(pos) ** 3)) / k


def update_acceleration(pos, planet, mass):
    mu = 1.0
    rad = pos - planet.pos
    return -1.0*mu*pos/(mag(pos)**3)- mass*rad/(mag(rad) ** 3)

earth_period = find_period(earth_position)
mars_period = find_period(mars_position)


while True:
#while tau <= end_tau:
    rate(10000)
    earth_trail.append(earth_position)
    earth.pos = earth_position
    earth_position += earth_velocity * timestep
    earth_velocity += update_acceleration(earth_position, mars, mars_mass) * timestep
    mars_trail.append(mars_position)
    mars.pos = mars_position
    mars_position += mars_velocity * timestep
    mars_velocity += update_acceleration(mars_position, earth, earth_mass) * timestep
    tau += timestep

print earth_period, mars_period
