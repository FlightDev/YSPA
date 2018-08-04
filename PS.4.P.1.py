import matplotlib.pyplot as plt
from math import *

x_coor = [1.1, 1.6, 2.0, 2.1, 2.9, 3.2, 3.3, 4.4, 4.9]
y_coor = [72.61, 72.91, 73.00, 73.11, 73.52, 73.70, 76.10, 74.26, 74.51]
plt.scatter(x_coor, y_coor)


tx = [2, 3, 5, 7, 9]
ty = [4, 5, 7, 10, 15]

def fit_least_squares(x_vals, y_vals):
    x = 0.0
    y = 0.0
    xy = 0.0
    x_sq = 0.0
    N = 0.0
    x_mean_x = 0.0
    x_mean_y = 0.0
    for i in range(len(x_vals)):
        x += x_vals[i]
        y += y_vals[i]
        xy += x_vals[i] * y_vals[i]
        x_sq += x_vals[i] ** 2
        N += 1
    mean_x = x / N
    mean_y = y / N
    for i in range(len(x_vals)):
        x_mean_y += mean_y * x_vals[i]
        x_mean_x += mean_x * x_vals[i]
    slope = (xy - x_mean_y)/(x_sq - x_mean_x)
    intercept = (y - slope * x)/(N)
    return slope, intercept

def find_val(x, slope, intercept):
    return slope * x + intercept

def standard_deviation(x_vals, y_vals, slope, intercept):
    diff_sum = 0
    for i in range(len(x_vals)):
        diff_sum += (y_vals[i] - find_val(x_vals[i], slope, intercept)) ** 2
    diff_sum /= (len(x_vals) - 2)
    return sqrt(diff_sum)

def throw_outliers(x_coor, y_coor, slope, intercept, standard_deviation):
    for i in range(len(x_coor) - 1, -1, -1):
        if y_coor[i] > find_val(x_coor[i], slope, intercept) + 3 * sd or y_coor[i] < find_val(x_coor[i], slope, intercept) - 3 * sd:
            del x_coor[i]
            del y_coor[i]
            print x_coor[i], y_coor[i]
    return x_coor, y_coor

def get_line_points(x_coor, slope, intercept):
    x1 = x_coor[0]
    x2 = x_coor[len(x_coor ) - 1]
    y1 = find_val(x1, slope, intercept)
    y2 = find_val(x2, slope, intercept)
    return x1, y1, x2, y2

line = fit_least_squares(x_coor, y_coor)
sd = standard_deviation(x_coor, y_coor, line[0], line[1])
print sd
print find_val(3.5, line[0], line[1])
throw_outliers(x_coor, y_coor, line[0], line[1], sd)
lines = []
lines.append(get_line_points(x_coor, line[0], line[1]))
lines.append(get_line_points(x_coor, line[0], line[1] + 3 * sd))
lines.append(get_line_points(x_coor, line[0], line[1] - 3 * sd))
for line in lines:
    plt.plot([line[0], line[2]], [line[1], line[3]])
plt.show()
