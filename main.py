import numpy as np


# Define custom colors for printing
class bcolors:
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

def print_matrix(matrix):
    for row in matrix:   # Loop for iteration in each row
        for element in row:   # Loop for iteration in all the elements in each row
            print(element, end=" ")  # Print each element followed by a space
        print()  # Prints an empty to move to the next row
    print()  # Add an extra empty line after printing the entire matrix, providing visual separation.



def swap_rows_elementary_matrix(n, row1, row2):
    elementary_matrix = np.identity(n)
    elementary_matrix[[row1, row2]] = elementary_matrix[[row2, row1]]

    return np.array(elementary_matrix)



def row_addition_elementary_matrix(n, target_row, source_row, scalar=1.0):

    if target_row < 0 or source_row < 0 or target_row >= n or source_row >= n:
        raise ValueError("Invalid row indices.")

    if target_row == source_row:
        raise ValueError("Source and target rows cannot be the same.")

    elementary_matrix = np.identity(n)
    elementary_matrix[target_row, source_row] = scalar

    return np.array(elementary_matrix)



# Partial Pivoting: Find the pivot row with the largest absolute value in the current column
def partial_pivoting(A,i,N):
    pivot_row = i
    v_max = A[pivot_row][i]
    for j in range(i + 1, N):
        if abs(A[j][i]) > abs(v_max):#just in check need to use with abs
            """1 change"""
            v_max = A[j][i]
            pivot_row = j

    # if a principal diagonal element is zero,it denotes that matrix is singular,
    # and will lead to a division-by-zero later.
    if not A[pivot_row][i]:
        return "Singular Matrix"


    # Swap the current row with the pivot row
    if pivot_row != i:
        swap_row(A, i, pivot_row)
    return A

# function for elementary operation of swapping two rows
def swap_row(mat, i, j):
    N = len(mat)
    for k in range(N + 1):
        temp = mat[i][k]
        mat[i][k] = mat[j][k]
        mat[j][k] = temp
def backward_substitution(mat):
    #N = len(mat)
    N, M= np.array(mat).shape #M- number of columns, N- number of lines
    x = np.zeros(N)  # An array to store solution

    # Start calculating from last equation up to the first
    for i in range(N - 1, -1, -1):

        x[i] = mat[i][M-1]

        # Initialize j to i+1 since matrix is upper triangular
        for j in range(i + 1, N):
            x[i] -= mat[i][j] * x[j]

        x[i] = (x[i] / mat[i][i])
    return x
def norm(mat):
    size = len(mat)  # size of matrix
    max_row = 0  # In max_row will be the max sum between rows- יכיל את סכום השורה הגבוה ביותר
    for row in range(size):
        sum_row = 0  # In sum_row will be the sum of each row -סכום כל שורה
        for col in range(size):  # On this for make sum to one row- בלולאה זו סוכמים את השורה
            sum_row += abs(mat[row][col])  # Adding element to sum - מוסיפים את הערך המוחלט של האיבר לסכום השורה בה הוא נמצא
        if sum_row > max_row:  # check if the sum of current row is bigger than sum of previous row- בודקים האם הסכום שחישבנו עכשיו גדול יותר מהסכום השמור
            max_row = sum_row  # enter to max_row the biggest sum- עדכון סכום השורה הגבוה
    return max_row  #return max sum of row-מחזירים את סכום השורה הגבוה ביותר

def gaussianElimination(mat):
    """if np.linalg.det(mat)==0:
        print("The matrix is singular")
        return"""
    N = len(mat)

    singular_flag = forward_substitution(mat)


    if singular_flag != -1:

        if mat[singular_flag][N]:
            return "Singular Matrix (Inconsistent System)"
        else:
            return "Singular Matrix (May have infinitely many solutions)"

    # if matrix is non-singular: get solution to system using backward substitution
    # אם המטריצה אינה סינגולרית: קבל פתרון למערכת באמצעות התאמה לאחור

    return backward_substitution(mat)


def forward_substitution(mat):
    N=len(mat)
    A=mat.copy()
    for k in range (N):
        result1 = partial_pivoting(A, k, N)
        if isinstance(result1, str):  # Check if partial pivoting returns an error message
            print(result1)
            return k  # Matrix is singular
        #print("Matrix after pivoting: ", result1)
        A = result1  # Update A with the result of partial pivoting

        for i in range(k + 1, N):
          m = -A[i][k] / A[k][k] #the multiple
          B=row_addition_elementary_matrix(N, i, k,  m)
          C=np.dot(B, A)
          #C=matrix_multiply(B, A)

          """print('Elementary matrix:')
          print_matrix(B)
          print('*')
          print('Original matrix:')
          print_matrix(A)
          print('=')
          print('Result matrix:')
          print_matrix(C)
          print("------------------------------------------------------------------")"""

          A = C
    mat[:] = A.tolist()  # Update the original matrix with the modified one
    return -1
def is_diagonally_dominant(mat):
    if mat is None:
        return False

    d = np.diag(np.abs(mat))  # Find diagonal coefficients-מחזיר את איברי הציר
    s = np.sum(np.abs(mat), axis=1) - d  # sum() of numpy- made sum of each line (axis=1). Find row sum without diagonal
    return (d > s).all()

def gauss_seidel(A, b, X0, TOL=1e-16, N=200):
    n = len(A)
    k = 1

    if not is_diagonally_dominant(A):
        print("The matrix has not been diagonally Notice to result!!")
        return None
    else:
        print('Matrix is diagonally dominant - preforming gauss seidel algorithm\n')

    #print( "Iteration" + "\t\t\t".join([" {:>12}".format(var) for var in ["x{}".format(i) for i in range(1, len(A) + 1)]]))
    print("-----------------------------------------------------------------------------------------------")
    x = np.zeros(n, dtype=np.double)
    while k <= N:

        for i in range(n):
            sigma = 0
            for j in range(n):
                if j != i:
                    sigma += A[i][j] * x[j]
            x[i] = (b[i] - sigma) / A[i][i]

        #print("{:<15} ".format(k) + "\t\t".join(["{:<15} ".format(val) for val in x]))

        if norm(x - X0, np.inf) < TOL:
            return tuple(x)

        k += 1
        X0 = x.copy()

    print("Maximum number of iterations exceeded")
    return tuple(x)



def linearInterpolation(table_points, point):
    """This function check the behavior of the function at point that is not in the table_points
    :param table_points: table of points which are in the function
    :param point: the point to interpolate- A point that is not in the table, and we would like to evaluate the behavior
     of the function at this point
    :return: None
    """
    if len(table_points) < 2:
        #print("At least two points are required for interpolation.")
        linear_Interpolation_for_two_points(table_points, point)
        return 0

    for i in range(len(table_points) - 1):
        # linear interpolation is just for point that in range of table point
        #if i <= point <= i + 1:
        if table_points[i][0] <= point <= table_points[i + 1][0]:

            x1 = table_points[i][0]
            x2 = table_points[i + 1][0]
            y1 = table_points[i][1]
            y2 = table_points[i + 1][1]
            result = (((y1 - y2) / (x1 - x2)) * point) + ((y2 * x1) - (y1 * x2)) / (x1 - x2)
            print(bcolors.OKGREEN, "\nThe approximation (interpolation) of the point ", point, " is: ", bcolors.ENDC, round(result, 4))

            return 1

    print("The point is not between table_points values, check with another method")
    return 0

def linear_Interpolation_for_two_points(table_points, x):

    x0 = table_points[0][0]
    x1 = table_points[1][0]
    y0 = table_points[0][1]
    y1 = table_points[1][1]

    # Calculating interpolated value
    yp = y0 + ((y1-y0)/(x1-x0)) * (x - x0)

    # Displaying result
    print(bcolors.OKBLUE, 'Interpolated value at %0.4f is %0.4f' %(x,yp), bcolors.ENDC)




def solveMatrix(matrixA,vectorb):
    detA = np.linalg.det(matrixA)
    #print(bcolors.YELLOWBG, "\nDET(A) = ", detA)

    if detA != 0:
        #print("CondA = ", Cond(matrixA, InverseMatrix(matrixA, vectorb)), bcolors.ENDC)
        print(bcolors.OKBLUE, "\nnon-Singular Matrix - Perform GaussJordanElimination",bcolors.ENDC)
        for line in range(len(matrixA)):  # Adding vectorb to matrix
            matrixA[line].append(vectorb[line][0])
        #print(matrixA)
        result = gaussianElimination(matrixA)
        #print("The result is: ", np.array(result))
        return result
    else:
        X0 = np.zeros_like(vectorb)
        result = gauss_seidel(matrixA, vectorb, X0, 1e-16, 200)
        if result is None:
            print("Failed to converge. Adjust the initial guess or matrix A.")
            return None
        vector_result =[]
        for element in result:
            vector_result.append(element)
        print("The result is: ", np.array(vector_result))
        return vector_result
        """print("Singular Matrix - Perform LU Decomposition\n")
        print("Matrix U: \n")
        print(np.array(UMatrix(matrixA, vectorb)))
        print("\nMatrix L: \n")
        print(np.array(LMatrix(matrixA, vectorb)))
        print("\nMatrix A=LU: \n")
        result = matrix_multiply(LMatrix(matrixA, vectorb), UMatrix(matrixA, vectorb)) #change
        print(np.array(result))
        return result"""


def polynomialInterpolation(table_points, x):
    matrix = [[point[0] ** i for i in range(len(table_points))] for point in table_points] # Makes the initial matrix

    b = [[point[1]] for point in table_points]

    print(bcolors.OKBLUE, "The matrix obtained from the points: ", bcolors.ENDC,'\n', np.array(matrix))
    print(bcolors.OKBLUE, "\nb vector: ", bcolors.ENDC,'\n',np.array(b))
    matrixSol = solveMatrix(matrix, b)
    if matrixSol is None:
        print(bcolors.OKBLUE,"\nNo solution !", bcolors.ENDC)
        return

    result = sum([matrixSol[i] * (x ** i) for i in range(len(matrixSol))])
    print(bcolors.OKBLUE, "\nThe polynom:", bcolors.ENDC)
    print('P(X) = '+'+'.join([ '('+str(matrixSol[i])+') * x^' + str(i) + ' ' for i in range(len(matrixSol))])  )
    print(bcolors.OKGREEN, f"\nThe Result of P(X={x}) is:", bcolors.ENDC)
    print(result)
    return result


def lagrange_interpolation(table_points, x):
    """
    Lagrange Interpolation

    Parameters:
    x_data (list): List of x-values for data points.
    y_data (list): List of y-values for data points.
    x (float): The x-value where you want to evaluate the interpolated polynomial.

    Returns:
    float: The interpolated y-value at the given x.
    """
    n = len(table_points)
    result = 0.0

    for i in range(n):

        term = 1
        #term = y_data[i]

        for j in range(n):

            if i != j :
                if table_points[i][0] == table_points[j][0]:
                    print("Impossible to interpolate!")
                    return None

                term *= (x - table_points[j][0]) / (table_points[i][0] - table_points[j][0])
        result += term * table_points[i][1]
        #result += term
    result1= round(result, 1)

    if result is not None:
         print(bcolors.OKBLUE, "\nInterpolated value at x =", x, "is y =", result1, bcolors.ENDC)


    return result1

if __name__ == '__main__':

    table_points = [(-5, -2), (-1, 6)]
    x = -2
    table_points = [(0, 0),(1, 0.8415),(2, 0.9093),(3, 0.1411),(4, -0.7568),(5, -0.9589),(6, -0.2794)]
    x = 1.28
    table_points = [(1.2, -1.2),(1.3, -2.3),(1.4, -0.5),(1.5, 0.89),(1.6, 1.37)]
    xa = 1.35
    xb=1.55



    print(bcolors.OKBLUE, "Interpolation & Extrapolation Methods\n", bcolors.ENDC)
    print(bcolors.OKBLUE, "Table Points", table_points, bcolors.ENDC)
    print(bcolors.OKGREEN, "Finding an approximation to the point ", x, bcolors.ENDC)
    choice = int(input(
        "Which method do you want? \n\t1.Linear Method \n\t2.Polynomial Method\n\t3.Lagrange Method\n"))
    #choice = int(input(
        #"Which method do you want? \n\t1.Linear Method \n\t2.Polynomial Method\n\t3.Lagrange Method\n\t4.cubicSpline Method\n"))

    if choice == 1:
        linearInterpolation(table_points, x)
    elif choice == 2:
        polynomialInterpolation(table_points, x)
    elif choice == 3:
        lagrange_interpolation(table_points, x)
    #elif choice == 4:
        #cubicSpline(table_points,x)

    else:
        print(bcolors.FAIL, "Invalid input", bcolors.ENDC)



