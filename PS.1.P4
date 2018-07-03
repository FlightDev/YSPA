#Problem Set 1 - Problem 4

def determinant(matrix):
    dimension = len(matrix)
    for row in matrix: #checks if given a square matrix or not
        if dimension != len(row):
            return -1
    if len(matrix) == 1 and len(matrix[0]) == 1: #if base case, return the number
        return matrix[0][0]
    negative = -1
    sum = 0
    for column_number in range(0, len(matrix)):
        negative *= -1 #creates a switch to tell if a factor needs to be added or subtracted
        temp_matrix = []
        for row in range (1, len(matrix)):
            temp_array = []
            for number in range(0, len(matrix)): #iterates through each element in the top row
                if number != column_number: #creates a temporary matrix of all numbers not in the chosen number's row or column
                    temp_array.append(matrix[row][number])
            temp_matrix.append(temp_array)
        sum += negative * matrix[0][column_number] * determinant(temp_matrix) #adds the factor to the determinant of the matrix
    return sum


print (determinant([[1, 3, 5, 9],
                    [1, 3, 1, 7],
                    [4, 3, 9, 7],
                    [5, 2, 0, 9]]))
