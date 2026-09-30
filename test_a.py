"""
Test A
1D slab, 0<x<8 cm, 5 regions

left boundary: reflective
right boundary: vacuum

S_8 Gauss-Legendre Quadrature
40 uniform spatial cells

epsilon = 1e-6

Region 1: [0, 2] cm
    sigma_t = 50 cm-1
    sigma_s = sigma_f = nu_f = 0
    q       = 50 cm-3 s-1
Region 2: [2, 3] cm
    sigma_t = 5 cm-1
    sigma_s = sigma_f = nu_f = 0
    q       = 0
Region 3: [3, 5] cm
    (vacuum)
Region 4: [5, 6] cm
    sigma_t = 1.1 cm-1
    sigma_s = 0.7 cm-1
    sigma_f = 0.1 cm-1
    nu_f    = 2
    q       = 0
Region 5: [6, 8] cm
    sigma_t = 1.1 cm-1
    sigma_s = 0.5 cm-1
    sigma_f = 0.16 cm-1
    nu      = 2.5 cm-1
    q       = 1 cm-3, s-1
"""

import numpy as np

a = 8
n_x = 40
n_mu = 8
epsilon_phi = 1e-6
is_reflective = True


# spatial mesh
x_boundary = np.linspace(0, a, n_x+1)
dx = x_boundary[1:] - x_boundary[:-1]
x_center = 0.5*(x_boundary[1:] + x_boundary[:-1])

# angular quadrature
[mu, w] = np.polynomial.legendre.leggauss(n_mu)

# boundary condition
# (reflective bc left)
bc_right = np.zeros(n_mu/2)

# xs assignment
sigma_t = np.zeros(n_x)
sigma_s = np.zeros(n_x)
sigma_f = np.zeros(n_x)
nu = np.zeros(n_x)
q = np.zeros(n_x)

# Region 1
range = np.nonzero( (x_center < 2))
sigma_t[range] = 50
q[range] = 50
# Region 2
range = np.nonzero((x_center > 2) * (x_center < 3))
sigma_t[range] = 5
# Region 3
range = np.nonzero((x_center > 2) * (x_center < 3))
# Region 4
range = np.nonzero((x_center > 5) * (x_center < 6))
sigma_t[range] = 1.1
sigma_s[range] = 0.7
sigma_f[range] = 0.1
nu[range] = 2
# Region 5
range = np.nonzero((x_center > 6))
sigma_t[range] = 1.1
sigma_s[range] = 0.5
sigma_f[range] = 0.16
nu[range] = 2.5
q[range] = 1
