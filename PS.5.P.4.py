from visual import *
from math import *


def euler_acceleration(pos, mu):
    return -1.0 * mu * pos / (mag(pos) ** 3)

def verlet_position(current_pos, last_pos, d_tau, mu):
    r1 = 2.0 * current_pos - last_pos + euler_acceleration(current_pos, mu) * (d_tau ** 2)
    return [r1, current_pos]


tau = 0.00
d_tau = 0.01
end_tau = 1000.0
mu = 1.0
k = 0.01720209895



r_last = vector(0.9000, 0.0000, 0.0000)
r_dot = vector(0.0000, 1.3000, 0.0000)
r_now = r_last + r_dot * d_tau + 0.5 * euler_acceleration(r_last, mu) * (d_tau ** 2)
r = vector(0.9000, 0.0000, 0.0000)
tau += d_tau
accel = euler_acceleration(r, mu)
r += r_dot * d_tau
r_dot += accel * d_tau

sun = sphere(color = color.yellow, pos = (0, 0, 0), radius = 0.3)
o1 = sphere(color = color.red, pos = (r.x, r.y, r.z), radius = 0.01)
o2 = sphere(color = color.blue, pos = (r_now.x, r_now.y, r_now.z), radius = 0.01)
o1_t = curve(color = color.red)
o2_t = curve(color = color.blue)

while True:
    rate(50)
    print r, r_now
    tau += d_tau
    accel = euler_acceleration(r, mu)
    r += r_dot * d_tau
    r_dot += accel * d_tau
    temp = verlet_position(r_now, r_last, d_tau, mu)
    r_last = temp[1]
    r_now = temp[0]
    o1.pos = (r.x, r.y, r.z)
    o2.pos = (r_now.x, r_now.y, r_now.z)
    o1_t.append((r.x, r.y, r.z))
    o2_t.append((r_now.x, r_now.y, r_now.z))
    if abs(mag(r) - mag(r_now))/ mag(r_now) >= 0.01:
        print 12345
        break
print tau/k, tau
