import matplotlib.pyplot as plt
from math import *

#Calibrated, Catalog, Red
x_coor = [10.665, 9.997, 9.377, 8.882, 8.849, 12.525, 11.974, 10.660, 11.978, 14.517, 14.72]
y_coor = [11.603, 10.439, 10.160, 9.438, 9.554, 14.867, 13.740, 11.833, 13.708, 15.630, 15.730]

#Calibrated, Catalog, V
#x_coor = [14.289, 14.022, 12.167, 12.208, 17.850, 13.222, 17.249, 16.389, 16.170]
#y_coor = [15.140, 14.270, 12.810, 11.582, 17.810, 13.840, 17.520, 16.690, 16.040]
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
print line[0], line[1]
print find_val(16.038, line[0], line[1])
#throw_outliers(x_coor, y_coor, line[0], line[1], sd)
lines = []
lines.append(get_line_points(x_coor, line[0], line[1]))
lines.append(get_line_points(x_coor, line[0], line[1] + 3 * sd))
lines.append(get_line_points(x_coor, line[0], line[1] - 3 * sd))
for line in lines:
    plt.plot([line[0], line[2]], [line[1], line[3]])
plt.show()
