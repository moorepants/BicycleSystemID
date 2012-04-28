#/usr/bin/env python

from sympy import Matrix, Symbol, eye, symbols

# This computes the one step ahead predictor for the state space with process
# noise.

A = Matrix(2, 2, lambda i, j: Symbol('a' + str(i + 1) + str(j + 1)))
B = Matrix(2, 1, lambda i, j: Symbol('b' + str(i + 1) + str(j + 1)))
C = Matrix(2, 2, lambda i, j: Symbol('c' + str(i + 1) + str(j + 1)))
K = Matrix(2, 2, lambda i, j: Symbol('k' + str(i + 1) + str(j + 1)))

u = Matrix(1, 1, lambda i, j: Symbol('u' + str(i + 1) + str(j + 1)))
y = Matrix(2, 1, lambda i, j: Symbol('y' + str(i + 1) + str(j + 1)))

q = Symbol('q')

yhat = C * (q * eye(2) - A + K * C).inv() * (B * u + K * y)

# At this point figuring out the polynomial coefficients of q would hopefully
# reveal a pattern. These coeffcients can be used when converting back to the
# time domain and setting up the cost function.

##############################################################################

# This is the example in Ljung 1999 that shows the relationship of the state
# space predictor to the ARMAX model from the observer canonical form.

a1, a2, a3 = symbols('a1 a2 a3')
A = Matrix([[-a1, 1, 0],
            [-a2, 0, 1],
            [-a3, 0, 0]])

b1, b2, b3 = symbols('b1 b2 b3')
B = Matrix([[b1],
            [b2],
            [b3]])

C = Matrix([[1, 0, 0]])

k1, k2, k3 = symbols('k1 k2 k3')
K = Matrix([[k1],
            [k2],
            [k3]])

u, y = symbols('u y')

yhat = C * (q * eye(3) - A + K * C).inv() * (B * u + K * y)
