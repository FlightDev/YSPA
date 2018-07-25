from astropy.io import fits
import matplotlib.pyplot as plt
file = "/home/student/Desktop/Combine2.fit"
image = fits.open(file)
image = image[0].data
print type(image)
plt.imshow(image, cmap='gray')
plt.show

x = 1198
y = 850


def center_pixel(matrix, x, y, size):
    if size % 2 == 0:
        return -1
    mat = []
    offset = (size - 1) / 2
    for i in range(size):
        temp = []
        for j in range(size):
            temp.append(matrix[y - offset + i][x - offset + j])
        mat.append(temp)
    return mat

def center_weight(matrix):
    sum_x_weight = 0
    sum_y_weight = 0
    sum_mass = 0
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            sum_x_weight += j * matrix[i][j]
            sum_y_weight += i * matrix[i][j]
            sum_mass += matrix[i][j]
    center_x = sum_x_weight * 1.0 / sum_mass
    center_y = sum_y_weight * 1.0 / sum_mass
    return [center_x, center_y, len(matrix) * len(matrix[0])]

def find_count(matrix):
    sum_mass = 0
    sum_counter = 0
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            sum_counter += 1
            sum_mass += matrix[i][j]
    return [sum_mass, sum_counter]

def subtract_count(matrix, scale):
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            matrix[i][j] -= scale
    return matrix

def find_centroid(matrix, x, y):
    aperature = center_pixel(matrix, x, y, 7)
    ring = center_pixel(matrix, x, y, 9)
    factor = (find_count(ring[0]) - find_count(aperature[0])) / (ring[1] - aperature[1])
    aperature = subtract_count(aperature[0])
    aperature = find_count(aperature)
    return aperature

print find_centroid(image, x, y)


warm_up = [[0, 33, 21, 33, 8],
           [0, 56, 51, 53, 26],
           [23, 120, 149, 73, 18],
           [55, 101, 116, 50, 16],
           [11, 78, 26, 2, 10]]
