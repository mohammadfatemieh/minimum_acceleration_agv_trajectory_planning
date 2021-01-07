import numpy
import scipy
from numpy import ones, zeros
from numpy import flip, fliplr, flipud
from numpy import poly1d, polyder, polyint
from numpy.polynomial.polynomial import polypow
from scipy.linalg import block_diag
from copy import deepcopy
import cvxopt

#   [-2, 2, 2],
#   [-2, 3, 3],
keyframes = numpy.array([
    [0, 0, 1],
    [0, 4, 4],
])

vel_bound = [1, 2]
acc_bound = [1, 2]

vel_state = [1, 1]
acc_state = [0, 0]

total_time = 10
time_step = 0.1

class sub_cost_matrix(object):
    def __init__(self, poly_order: int, derivative_order: int):
        self.poly_order = poly_order
        self.derivative_order = derivative_order
        self.poly: poly1d = polyder(numpy.ones(poly_order), derivative_order)
        self.matrix: numpy.array = zeros((poly_order, poly_order))
        poly_array: numpy.array = numpy.array(self.poly).reshape((1, len(self.poly)))
        self.matrix[0 : poly_order - derivative_order,
                    0 : poly_order - derivative_order] = \
            poly_array.transpose() * poly_array
        for i in range(0, poly_order - derivative_order + 1):
            self.matrix[i, 0 : poly_order - i] = \
                polyint(self.matrix[i, 0 : poly_order - i])[0 : -1] # throw away the constant term
        print(self.matrix)

    def __str__(self):
        return str(self.matrix)

    def substitute(self, time: tuple[float, float]):
        for i in range(0, self.poly_order):
            for j in range(0, i + 1):
                self.matrix[j, i - j] *= \
                    abs(
                        pow(time[0], (self.poly_order - self.derivative_order - 1) ** 2 - i + 1) - \
                        pow(time[1], (self.poly_order - self.derivative_order - 1) ** 2 - i + 1)
                    )
        return self

class cost_matrix(object):
    def __init__(self, poly_order: int, derivative_order: int, keyframe_list: numpy.array):
        self.poly_order = poly_order
        self.derivative_order = derivative_order
        self.matrix = sub_cost_matrix(poly_order, derivative_order).substitute((keyframe_list[0, 2], keyframe_list[1, 2])).matrix
        for i in range(1, len(keyframe_list) - 1):
            self.matrix = block_diag(
                self.matrix,
                sub_cost_matrix(poly_order, derivative_order).substitute((keyframe_list[i, 2], keyframe_list[i + 1, 2])).matrix
            )

    def __str__(self):
        return str(self.matrix)


def path_constrain(poly_order: int, keyframe_list: numpy.array):
    poly_coef_matrix = 0
    for t in range(0, len(keyframe_list) - 1):
        sub_poly_coef_matrix = [ones(poly_order), ones(poly_order)]
        for j in range(0, 2):
            for i in range(0, poly_order):
                sub_poly_coef_matrix[j][i] = pow(keyframe_list[t + j, 2], poly_order - i - 1)
        if poly_coef_matrix.__class__ == int:
            poly_coef_matrix = deepcopy(sub_poly_coef_matrix)
        else:
            poly_coef_matrix = block_diag(poly_coef_matrix, sub_poly_coef_matrix)
    b = [keyframe_list[0, 0]]
    for i in range(1, len(keyframe_list) - 1):
        b.append(keyframe_list[i, 0])
        b.append(keyframe_list[i, 0])
    b.append(keyframe_list[-1, 0])
    return (poly_coef_matrix, numpy.array(b).reshape(len(b), 1))

def constrain(poly_order: int, keyframe_list: numpy.array):
    def substitute(p: poly1d, time: float):
        for i in range(0, len(p)):
            p[i] *= pow(time, len(p) - 1 - i)
        return p
    A, b = path_constrain(5, keyframe_list)
    A = numpy.concatenate((
        A,
        numpy.array([
            numpy.pad(substitute(polyder(ones(poly_order), 1), keyframe_list[ 0, 2]), (0, (len(keyframe_list) - 2) * poly_order + 1), constant_values=0),
            numpy.pad(substitute(polyder(ones(poly_order), 1), keyframe_list[-1, 2]), ((len(keyframe_list) - 2) * poly_order + 0, 1), constant_values=0),
            numpy.pad(substitute(polyder(ones(poly_order), 2), keyframe_list[ 0, 2]), (0, (len(keyframe_list) - 2) * poly_order + 2), constant_values=0),
            numpy.pad(substitute(polyder(ones(poly_order), 2), keyframe_list[-1, 2]), ((len(keyframe_list) - 2) * poly_order + 0, 2), constant_values=0),
        ])
    ))
    b = numpy.concatenate((
        b,
        numpy.array([
            [1],
            [1],
            [0],
            [0],
        ]),
    ))
    return A, b

def cvxopt_solve_qp(P, q, G=None, h=None, A=None, b=None):
    P = .5 * (P + P.transpose())  # make sure P is symmetric
    args = [cvxopt.matrix(P), cvxopt.matrix(q)]
    if G is not None:
        args.extend([cvxopt.matrix(G), cvxopt.matrix(h)])
        if A is not None:
            args.extend([cvxopt.matrix(A), cvxopt.matrix(b)])
    sol = cvxopt.solvers.qp(*args)
    if 'optimal' not in sol['status']:
        return None
    return numpy.array(sol['x']).reshape((P.shape[1],))

if __name__ == '__main__':
    P = cost_matrix(5, 2, keyframes).matrix
    A, b = constrain(5, keyframes)
    print(numpy.linalg.matrix_rank(A))
    print(numpy.linalg.matrix_rank(numpy.concatenate([P, A])))
    print(numpy.linalg.eigvals(P))
    print(0.5 * (P + P.transpose()))
    print(A)
    print(b)
    q = zeros((5, 1))
    print(q)
    cvxopt_solve_qp(P, zeros((5, 1)), A=A, b=b)
