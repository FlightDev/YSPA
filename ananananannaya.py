#Creates the Prelim orbit for our orbit
#Uses the method of Gauss and an angles only formula to numerically integrate

from math import *
import numpy as np
import de421
from jplephem import Ephemeris
from visual import *

def acceleration (position, position2):
    a = -1.0*position/((mag(position))**3) - 0.000003003*(position2-position)/((mag(position2-position))**3)
    return a

#Ephem generation code to generate R vector
def EphemR (JulianDate):
    #Import Ephem 421 to get vector from Earth to Sun (R)
    AU = 149597870.7 # km/AU

    eph = Ephemeris(de421)
    barycenter = eph.position('earthmoon', JulianDate)
    moonvector = eph.position('moon', JulianDate)
    earth = barycenter - moonvector * eph.earth_share
    R0 = vector(earth)/(AU)  #This is the Sun to Earth center vector in Equatorial system
    R_geocentric = -1.*R0  #Multiply by -1 to get Earth to Sun
    #gg =                #Our geocentric vector

    #R_topocentric =  R_geocentric - gg  #Subtract our topocentric vector
    return R_geocentric

#Defining the f series (Taylor Series)
def fSeries (r, timef):
    f1 = 1
    f2 = (timef**2/(2*(mag(r))**3))
    ftotal = f1-f2
    return ftotal

#Defining the g series (Taylor Series)
def gSeries (r, timeg):
    g1 = timeg
    g2 = (timeg**3)/(6*(mag(r)**3))
    gtotal = g1-g2
    return gtotal

def findRho (a1, a3, rho1, rho2, rho3, R1, R2, R3):
    Rho1hat = norm(rho1)
    Rho2hat = norm(rho2)
    Rho3hat = norm(rho3)
    Rho1 = ((a1*(dot(cross(R1, Rho2hat), Rho3hat))) - dot(cross(R2, Rho2hat), Rho3hat) + a3*dot(cross(R3, Rho2hat), Rho3hat))/(a1*dot(cross(Rho1hat, Rho2hat), Rho3hat))
    Rho2 = ((a1*(dot(cross(Rho1hat, R1), Rho3hat))) - dot(cross(Rho1hat, R2), Rho3hat) + a3*dot(cross(Rho1hat, R3), Rho3hat))/(-1*dot(cross(Rho1hat, Rho2hat), Rho3hat))
    Rho3 = ((a1*(dot(cross(Rho2hat, R1), Rho1hat))) - dot(cross(Rho2hat, R2), Rho1hat) + a3*dot(cross(Rho2hat, R3), Rho1hat))/(a3*dot(cross(Rho2hat, Rho3hat), Rho1hat))
    return Rho1, Rho2, Rho3

def find_residuals(r, rdot, inputJulian, dec_to_check, ra_to_check, timePass):

    #Verlet Integrator
    k = 0.01720209895 # Gaussian gravitational consgtant
    dtau = 0.001  # This is Gaussian Days
    dt = dtau / k     # This is in mean solar days
    tao = 0.0
    deltao = 0.001

    r = r2guess
    rdot = r2Dot

    r0Verlet = r
    r1Verlet = r0Verlet + (rdot * deltao) + (0.5*acceleration(r0Verlet, Rho2hat*Rho2)*(deltao**2))

    while tao <= timePass*k:
        r2Update = (2*r1Verlet) - r0Verlet + acceleration(r1Verlet, Rho2hat*Rho2)*(deltao**2)
        r0Verlet = r1Verlet
        r1Verlet = r2Update

        tao += deltao

    #print r1Verlet


    #Ephemeris Generation

    AU = 149597870.7 # km/AU

    eph = Ephemeris(de421)
    barycenter = eph.position('earthmoon', inputJulian)
    moonvector = eph.position('moon', inputJulian)
    earth = barycenter - moonvector * eph.earth_share
    R0 = vector(earth)/AU # This is the Sun to Earth center vector in Equatorial system
    #print R0
    R_geocentric = -1.0*R0

    #gg = ?

    rInput = r1Verlet

    rho = (rInput + R_geocentric)
    rhohat = norm(rho)
    declination = arcsin(rhohat.z)
    RA = arccos(rhohat.x/(cos(declination)))
    if rhohat.y < 0:
        RA = 2*math.pi - RA

    RA_residuals = abs(RA - ra_to_check)/ra_to_check
    Dec_residuals = abs(declination - dec_to_check)/dec_to_check

    return RA_residuals, Dec_residuals

    #print "RA = ", degrees(RA)
    #print "Dec = ", degrees(declination)



#Define Time intervals (t and tau)
k = 0.01720209895 # Gaussian gravitational consgtant
epsilon = radians(23.43688)

time1 = 2458312.750000 #7-19 5 combo
time2 = 2458315.750000 #7-24 3.7
time3 = 2458324.750000 #7-29 1comb

tau1 = k*(time1-time2)  #Modified time (Tau2 = middle reference)
tau3 = k*(time3-time2)

#Getting the vector from the Earth to the asteroid (rho)
Obs1 = [radians(334.0454167), radians(3.5653)]
Obs2 = [radians(336.465), radians(5.1900)]
Obs3 = [radians(345.0720833), radians(10.7469)]

#Puttting in observational data and returning rho vectors
Rho1hat = vector((np.cos(Obs1[0])*np.cos(Obs1[1])),(np.sin(Obs1[0])*np.cos(Obs1[1])), (np.sin(Obs1[1])))
Rho2hat = vector((np.cos(Obs2[0])*np.cos(Obs2[1])),(np.sin(Obs2[0])*np.cos(Obs2[1])), (np.sin(Obs2[1])))
Rho3hat = vector((np.cos(Obs3[0])*np.cos(Obs3[1])),(np.sin(Obs3[0])*np.cos(Obs3[1])), (np.sin(Obs3[1])))

R1 = (EphemR(time1))
R2 = (EphemR(time2))
R3 = (EphemR(time3))

#Rho2mag_new = vector(1.7,0,0)
#Rho2mag_old = vector(0,0,0)
r2guess = vector(1.5,0,0)

for i in range (0, 10):
    f1 = fSeries(r2guess, tau1)
    g1 = gSeries(r2guess, tau1)
    f3 = fSeries(r2guess, tau3)
    g3 = gSeries(r2guess, tau3)
    a1 = g3/(f1*g3 - f3*g1)
    a3 = (-g1)/(f1*g3 - f3*g1)
    Rho1, Rho2, Rho3 = findRho(a1, a3, Rho1hat, Rho2hat, Rho3hat, R1, R2, R3)
    r2guess = Rho2*Rho2hat - R2
    r1 = Rho1 * Rho1hat - R1
    r3 = Rho3 * Rho3hat - R3
    r2Dot = (f3/(g1*f3 - g3*f1))*r1 - (f1/(g1*f3 - g3*f1))*r3

#End of method of Gauss

"""
print Rho1, Rho2, Rho3
print "R2 = ", r2guess
print "R2mag = ", mag(r2guess)
print "R2 Dot = ", r2Dot
print "Rho2 = ", Rho2*Rho2hat
"""
last_vector = vector(1000, 1000, 1000)
last_total_resid = 100000000
i = 0
while (True):
    RAresid, Decresid = find_residuals(r2guess, r2guess, 2458321.750000, 8.7783, 341.94375, 6)
    i += 1
    print i, RAresid + Decresid
    resid = RAresid + Decresid
    if RAresid + Decresid < 1:
        break
    last_vector = r2guess
    r2guess = vector(np.random.lognormal(),np.random.lognormal(),np.random.lognormal()) + r2guess
    f1 = fSeries(r2guess, tau1)
    g1 = gSeries(r2guess, tau1)
    f3 = fSeries(r2guess, tau3)
    g3 = gSeries(r2guess, tau3)
    a1 = g3/(f1*g3 - f3*g1)
    a3 = (-g1)/(f1*g3 - f3*g1)
    Rho1, Rho2, Rho3 = findRho(a1, a3, Rho1hat, Rho2hat, Rho3hat, R1, R2, R3)
    r2guess = Rho2*Rho2hat - R2
    r1 = Rho1 * Rho1hat - R1
    r3 = Rho3 * Rho3hat - R3
    r2Dot = (f3/(g1*f3 - g3*f1))*r1 - (f1/(g1*f3 - g3*f1))*r3

print r2guess
print r2Dot
