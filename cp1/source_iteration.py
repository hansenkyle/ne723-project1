import numpy as np
from cp1.cornervalues import CornerValues
from cp1.quadrature import Quadrature
from cp1.sweep import sweep
from cp1.xs import CrossSection


# solve the fixed-source problem with source iteration
def solve_fixed_source(
    epsilon_phi : float,
    phi_initial : CornerValues,
    q : np.array,
    xs : CrossSection,
    quad : Quadrature,
    dx : np.array,
    BC : np.array,
    lhs_reflective : bool
):
    # initialize guess for scalar flux
    phi = phi_initial.copy()

    delta = np.array([1.0])
    rho = np.array([1.0])

    # while not converged:
    l = 0
    while (delta[-1] > epsilon_phi * ((1/rho[-1]) - 1)):
        l += 1
        print(f"Iteration {l}...")
        # solve for angular and scalar flux given known scattering source
        psi, phi_new = sweep(phi, q, xs, quad, dx, BC, lhs_reflective)

        # compute Delta
        delta = np.append(delta, (phi_new - phi).norm_inf())
        rho = np.append(rho, delta[l]/delta[l-1])
        phi = phi_new.copy()

    # compute residuals
    #