from numpy import matrix
from numpy import linalg
from math import *


dec = [31.241614, 31.261435, 31.658777, 31.389536, 31.356873, 31.354845, 31.189566, 31.20734, 31.207386]
ra = [293.312693, 293.876543, 293.562800, 293.757655, 293.757094, 293.772857, 293.988409, 294.034602, 294.049422]
x = [337.306, 440.775, 753.332, 1056.967, 1055.110, 1080.333, 1422.026, 1496.856, 1520.755]
y = [175.383, 210.925, 951.820, 440.622, 379.151, 374.455, 58.266, 91.164, 90.815]
"""
dec = [42.693616667, 42.693608333, 42.713777778, 42.651097222, 42.626066667, 42.587472222, 42.664613889, 42.742877778]
ra =[168.286625, 168.25279167, 168.13083333, 168.21608333, 168.14645833, 168.12245833, 168.43845833, 168.25616667]
x = [312.40, 357.72, 516.36, 415.99, 514.57, 555.22, 115.79, 342.69]
y = [292.33, 287.08, 230.86, 358.72, 393.37, 459.89, 368.98, 197.45]
"""
def find_residuals(ra, dec, x, y, ra_stuff, dec_stuff):
    ras = []
    decs = []
    sigma_ra = 0
    sigma_dec = 0
    for i in range(len(ra)):
        ras.append((ra_stuff[0] + ra_stuff[1] * x[i] + ra_stuff[2] * y[i]) - ra[i])
        decs.append((dec_stuff[0] + dec_stuff[1] * x[i] + dec_stuff[2] * y[i]) - dec[i])
        sigma_ra += ((ra[i] - ra_stuff[0] - ra_stuff[1] * x[i] - ra_stuff[2] * y[i]) ** 2)
        sigma_dec += ((dec[i] - dec_stuff[0] - dec_stuff[1] * x[i] - dec_stuff[2] * y[i]) ** 2)
    sigma_ra = sqrt(sigma_ra / (len(ra) - 3))
    sigma_dec = sqrt(sigma_dec / (len(dec) - 3))
    print ras, decs, sigma_ra, sigma_dec

def find_needed_matricies(dec, ra, x, y):
    xi = 0.0
    yi = 0.0
    x_sq = 0.0
    xy = 0.0
    y_sq = 0.0
    a = 0.0
    ax = 0.0
    ay = 0.0
    declination = 0.0
    decx = 0.0
    decy = 0.0
    for i in range(len(dec)):
        xi += x[i]
        yi += y[i]
        x_sq += x[i] ** 2
        xy += x[i] * y[i]
        y_sq += y[i] ** 2
        a += ra[i]
        ax += ra[i] * x[i]
        ay += ra[i] * y[i]
        declination += dec[i]
        decx += dec[i] * x[i]
        decy += dec[i] * y[i]
    a1 = matrix([[a], [ax], [ay]])
    a2 = matrix([[len(ra), xi, yi],
                    [xi, x_sq, xy],
                    [yi, xy, y_sq]])
    d1 = matrix([[declination], [decx], [decy]])
    d2 = matrix([[len(dec), xi, yi],
                     [xi, x_sq, xy],
                     [yi, xy, y_sq]])
    return a1, a2, d1, d2

m = find_needed_matricies(dec, ra, x, y)
a1 = m[0]
a2 = m[1]
d1 = m[2]
d2 = m[3]

ra_c = linalg.solve(a2,a1)
dec_c = linalg.solve(d2,d1)

print ra_c
print dec_c
res = find_residuals(ra, dec, x, y, ra_c, dec_c)
print res
