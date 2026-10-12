"""
Test B
1D slab, 0<x<11 cm, 2 regions

left boundary: 
    psi(mu) = 1 for mu>0
right boundary: vacuum

S_8 Gauss-Legendre Quadrature
10 0.1-cm cells, 10 1-cm cells.

Region 1: [0, 1] cm
    sigma_t = 2 cm-1
    sigma_s = sigma_f = nu_f = 0
    q       = 0
Region 2: [1, 11] cm
    sigma_t = 100 cm-1
    sigma_s = 100 cm-1
    sigma_f = nu_f = 0
    q       = 0
"""


import numpy as np

from cp1.cornervalues import CornerValues
from cp1.quadrature import Quadrature
from cp1.source_iteration import solve_fixed_source
from cp1.xs import CrossSection

a = 11
n_x = 20
n_mu = 8
epsilon_phi = 1e-6
is_reflective = False

# spatial mesh
r1_bound = np.linspace(0, 1, 11)
r2_bound = np.linspace(2, 11, 10)
x_boundary = np.concatenate((r1_bound, r2_bound))
dx = x_boundary[1:] - x_boundary[:-1]
x_center = 0.5*(x_boundary[1:] + x_boundary[:-1])

# angular quadrature
[mu, w] = np.polynomial.legendre.leggauss(n_mu)
quad = Quadrature(mu, w)

# boundary condition
bc = np.zeros(n_mu)
bc[mu>0] = 1.0

# xs assignment
sigma_t = np.zeros(n_x)
sigma_s = np.zeros(n_x)
sigma_f = np.zeros(n_x)
nu = np.zeros(n_x)
q = np.zeros(n_x)

# Region 1
range = np.nonzero(x_center < 1)
sigma_t[range] = 2
# Region 2
range = np.nonzero(x_center > 1)
sigma_t[range] = 100
sigma_s[range] = 100

xs = CrossSection(sigma_t, sigma_s, 0*sigma_s, 0*sigma_s)


phi_initial = CornerValues(np.zeros(n_x), np.zeros(n_x))


print("solving test A fixed source problem...")
solution, convergence = solve_fixed_source(epsilon_phi, phi_initial, q, xs, quad, dx, bc, is_reflective)

