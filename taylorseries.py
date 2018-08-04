import matplotlib.pyplot as plt
from math import *
import matplotlib.patches as mpatches

def factorial(x):
    if x == 1:
        return 1
    else:
        return x * factorial(x - 1)

def find_taylor(x, num_terms):
    sum = 0.0
    counter = 1.0
    for i in range(num_terms):
        sum += counter * (x ** (i * 2 + 1)) / (factorial(i * 2 + 1))
        counter *= -1.0
    return sum


x = []
sinx = []
taylorx = []
colors = ['b', 'g', 'r', 'c', 'm', 'y']
lines = ['dashed', 'dashdot', 'dotted', ':', '--', '-.']
for i in range(0, 65):
    x.append((i - 32) / 16.0 * pi)
    sinx.append(sin(x[i]))

for i in range(2, 7):
    temp = []
    for j in range(len(x)):
        temp.append(find_taylor(x[j], i))
    taylorx.append(temp)


plt.style.use('dark_background')
plt.xlabel('x input')
plt.ylabel('y algorithm output')
plt.plot(x, sinx, 'w', linewidth = 2)
plt.ylim(-50, 50)
handle = []
handle.append(mpatches.Patch(color='w', label='math.sin'))
for i in range(len(taylorx)):
    plt.plot(x, taylorx[i], colors[i], linewidth = 0.33*(i+2), linestyle = lines[i])
    handle.append(mpatches.Patch(color=colors[i], label='Taylor ' + str(i + 2) + ' terms', linestyle = lines[i], linewidth = 0.33*(i+2)))
plt.legend(handles=handle)
plt.savefig('Neal_Taylor_Sin_Plot.png')
plt.show()
