from visual import *

mu = 1.0
k = 0.01720209895
asteroidPos = vector(0.9, 0.0, 0.0)
asteroidOriginalPos = vector(0.9, 0.0, 0.0)
asteroidVel = vector(0.0, 1.3, 0.0)

xscale = 0.2
vscale = 0.05
xAxis = arrow(pos = (0, 0, 0), axis = (1, 0, 0), shaftwidth = 0.01, color = color.red)
yAxis = arrow(pos = (0, 0, 0), axis = (0, 1, 0), shaftwidth = 0.01, color = color.yellow)
zAxis = arrow(pos = (0, 0, 0), axis = (0, 0, 1), shaftwidth = 0.01, color = color.green)
sun = sphere(pos = vector(0.0, 0.0, 0.0), radius = 0.1)
planet = sphere(pos = asteroidPos, radius = 0.03, color = color.blue)
planet.velocity = asteroidVel
varr = arrow(pos = xscale * planet.pos, axis = vscale * planet.velocity, length = 0.21, shaftwidth = 0.01, color = color.red)
planetTrail = curve(color = color.blue)


def acceleration(position):
    return -1.0 * mu * position / (mag(position)**3)

def eulerIntegrator(position, velocity):
    time = 0.0
    dt = 0.01
    oldVelocity = velocity
    didntMove = True
    while (position.y > asteroidOriginalPos.y or didntMove):
        rate(100)
        didntMove = False
        acc = acceleration(position)
        oldVelocity = velocity
        velocity += acc * dt
        position += oldVelocity * dt
        planet.pos = position
        planet.velocity = velocity
        varr.pos = planet.pos
        varr.axis = xscale * planet.velocity
        planetTrail.append(planet.pos)
        time += dt
    return position, velocity, time


def stormerVerlet(position, velocity):
    time = 0.0
    dTau = 0.01
    pos0 = position
    print pos0
    pos1 = position + velocity * dTau + 0.5 * acceleration(position) * dTau * dTau
    didntMove = True
    while(position.y > asteroidOriginalPos.y or didntMove):
        didntMove = False
        position = 2 * pos1 - pos0 + acceleration(pos1) * dTau * dTau
        #planet.pos = position
        #planetTrail.append(position)
        pos0 = pos1
        pos1 = position
        time += dTau
    return position, time

def stormerVerletTime(position, position1, dT):
    return 2.0 * (position1 - position) + acceleration(position1) * dT * dT

'''
eulerasteroidList = eulerIntegrator(asteroidPos, asteroidVel)
print 'Position = ', eulerasteroidList[0]
print 'Velocity = ', eulerasteroidList[1]
print 'Time = ', eulerasteroidList[2]

asteroidPos = vector(0.9, 0.0, 0.0)
asteroidVel = vector(0.0, 1.3, 0.0)

stormerasteroidList = stormerVerlet(asteroidPos, asteroidVel)
print 'Position = ', stormerasteroidList[0]
print 'Time = ', stormerasteroidList[1]

'''

#percent error
time = 0.00
dT = 0.01

asteroidPos = vector(0.9000, 0.0000, 0.0000)
asteroidVel = vector(0.0000, 1.3000, 0.0000)

eulerPosition = vector(0.9000, 0.0000, 0.0000)
eulerVelocity = vector(0.0000, 1.3000, 0.0000)

stormerPosition = vector(0.9000, 0.0000, 0.0000)

a = acceleration(eulerPosition)
#eulerPosition += eulerVelocity * dT
#eulerVelocity += a * dT
stormerPosition1 = stormerPosition + asteroidVel * dT + 0.5 * acceleration(stormerPosition) * (dT ** 2)

#stormerPosition2 = stormerPosition1
print eulerPosition, eulerVelocity, stormerPosition, stormerPosition1

while True:
    print eulerPosition, stormerPosition, time
    if abs(mag(eulerPosition - stormerPosition)) / mag(stormerPosition) > 0.01:
        print time
        break
    #euler

    eulerPosition += eulerVelocity * dT
    acc = acceleration(eulerPosition)
    eulerVelocity += acc * dT

    #stormer
    temp = stromerPosition1
    stormerPosition2 = 2.0 * stormerPosition1 - stormerPosition + acceleration(stormerPosition1) * (dT ** 2)
    stormerPosition = stormerPosition1
    stormerPosition1 = stormerPosition2



    #increment time
    time += dT
    #print eulerPosition, stormerPosition1
