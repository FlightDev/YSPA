from astropy.io import fits
import matplotlib.pyplot as plt
file = "/home/student/Desktop/Combine2.fit"
image = fits.open(file)
image = image[0].data
print type(image)
plt.imshow(image, cmap='gray')
plt.show

x = 784
y = 75
def center_pixel(matrix, x, y):
    mat = []
    for i in range(7):
        temp = []
        for j in range(7):
            temp.append(matrix[y - 3 + i][x - 3 + j])
        mat.append(temp)
    return mat


def center_weight(matrix):
    print matrix
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
    return center_x, center_y

warm_up = [[0, 33, 21, 33, 8],
           [0, 56, 51, 53, 26],
           [23, 120, 149, 73, 18],
           [55, 101, 116, 50, 16],
           [11, 78, 26, 2, 10]]


print len(image)
m = center_pixel(image, x, y)
print center_weight(m)
