#!/usr/bin/env python

# This script calculates the numerators and denominator polynomials of the
# Whipple bicycle model tranfer functions as a function of the canonical matrix
# entries given in Meijaard2007.

from sympy import Symbol, Matrix, symbols, eye, zeros, roots, Poly, cos, sin

# v = speed, g = acceleration due to gravity, s = Laplace variable
v, g, s = symbols('v g s')

# build the canoncial matrices, setting some entries to zero
# I did not include these realtionships yet: M[0, 1] = M[1, 0] and K0[0, 1] =
# K0[1, 0]. It is possible that they help simplify things.
# I did set C[0, 0] = 0 and K2[0, 0] = K2[1, 0] = 0.
M = Matrix(2, 2, lambda i, j: Symbol('m' + str(i + 1) + str(j + 1)))
C1 = Matrix(2, 2, lambda i, j: 0 if i == 0 and j == 0 else
        Symbol('c1' + str(i + 1) + str(j + 1)))
K0 = Matrix(2, 2, lambda i, j: Symbol('k0' + str(i + 1) + str(j + 1)))
K2 = Matrix(2, 2, lambda i, j: 0 if j == 0 else
        Symbol('k2' + str(i + 1) + str(j + 1)))

# Build the A, B and C matrices such that steer torque is the only input and
# phi and delta are the only outputs.
Minv = M.inv()

Abl = -Minv * (g * K0 + v**2 * K2)
Abr = -Minv * v * C1

A = zeros(2).row_join(eye(2)).col_join(Abl.row_join(Abr))

#B = zeros((2, 1)).col_join(Minv[:, 1])
B = zeros(2).col_join(Minv)

C = eye(4)

apd = A[2, 1]
ITxx, ITxz, ITzz, IAlamx, IAlamlam, IAlamz, mu  = symbols('ITxx ITxz ITzz IAlamx IAlamlam IAlamz mu')
SA, ST, SF, mT, zT, w, lam = symbols('SA ST SF mT zT w lam')
test = apd.subs({M[0, 0]: ITxx,
          M[0, 1]: IAlamx + mu * ITxz,
          M[1, 0]: IAlamx + mu * ITxz,
          M[1, 1]: IAlamlam + 2 * mu * IAlamz + mu**2 * ITzz,
          K0[0, 1]: -SA,
          K0[1, 1]: -SA * sin(lam),
          K2[0, 1]: (ST - mT * zT) / w * cos(lam),
          K2[1, 1]: (SA + SF * sin(lam)) / w * cos(lam),
          })

## Transfer Functions ##

# Calculate the transfer function numerator polynomial. See Hoagg2007 for the
# equation.
N = C * (s * eye(4) - A).adjugate() * B
D = (s * eye(4) - A).det()

# These are the roll and steer numerator polynomials, respectively.
phiTphi = N[0, 0]
deltaTphi = N[1, 0]
phiTdelta = N[0, 1]
deltaTdelta = N[1, 1]

# Find the zeros
phiTphiZeros = roots(phiTphi, s, multiple=True)
deltaTphiZeros = roots(deltaTphi, s, multiple=True)
phiTdeltaZeros = roots(phiTdelta, s, multiple=True)
deltaTdeltaZeros = roots(deltaTdelta, s, multiple=True)

print('The characteristic equation.')
print(Poly(D, s))
print('\n')

inputs = ['Tphi', 'Tdelta']
outputs = ['phi', 'delta', 'phidot', 'deltadot']

for j, i in enumerate(inputs):
    for k, o in enumerate(outputs):
        print(i + '-' + o)
        print(Poly(N[k, j], s))
        print('\n')
