import matplotlib.pyplot as plt


def parse_file(s):
    ## TODO:
    dict = {}
    return dict

filename = "DATALOG6.TXT"

lines = []
text = open(filename, 'r')
for line in text:
    lines.append(line)

millisec = []
uSv_h = []
CPS = []
CPM = []
Temp = []
Alt = []
roll = []
pitch = []
heading = []
x_accel = []
y_accel = []
z_accel = []
balloonstate = []

for line in range(1110, 8100):
    s = lines[line].split(", ")
    if float(s[5]) < 0:
        continue
    millisec.append(int(s[0]))
    uSv_h.append(float(s[1]))
    CPS.append(float(s[2]))
    CPM.append(float(s[3]))
    Temp.append(float(s[4]))
    Alt.append(float(s[5]))
    roll.append(float(s[6]))
    pitch.append(float(s[7]))
    heading.append(float(s[8]))
    x_accel.append(float(s[9]))
    y_accel.append(float(s[10]))
    z_accel.append(float(s[11]))
    balloonstate.append(int(s[12][:1]))


#1025000-7900000
"""
millisec = millisec[1110:8100]
uSv_h = uSv_h[1110:8100]
Temp = Temp[1110:8100]
Alt = Alt[1110:8100]
"""
f, (x1, x2, x3) = plt.subplots(3)
x1.scatter(millisec, Alt)
x2.scatter(Alt, Temp)
x3.scatter(Alt, uSv_h)
plt.show()
