import numpy as np

from cp1.cornervalues import CornerValues
from cp1.quadrature import Quadrature
from cp1.sweep import sweep
from cp1.xs import CrossSection


class Solution:
    def __init__(self, phi, psi):
        self.phi = phi
        self.psi = psi

class Convergence:
    def __init__(self, delta, rho):
        self.delta = delta
        self.rho = rho

class Residuals:
    def __init__(self, balance, first_moment):
        self.balance = balance
        self.first_moment = first_moment


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
    while (delta[-1] > (epsilon_phi * ((1/rho[-1]) - 1)) or l<3):
        l += 1
        print(f"Iteration {l}...")
        # solve for angular and scalar flux given known scattering source
        psi, phi_new = sweep(phi, q, xs, quad, dx, BC, lhs_reflective)

        # compute Delta
        delta = np.append(delta, (phi_new - phi).norm_inf())
        rho = np.append(rho, delta[l]/delta[l-1])

        # if converged; compute residuals and exit.
        
        phi_old = phi.copy()
        phi = phi_new.copy()

        print(f"finished iteration {l}; Delta = {delta[-1]:.4e}, rho = {rho[-1]:.4e}")

    # cell-average angular flux
    cell_average_angular = np.array([psi_m.cell_average() for psi_m in psi])
    cell_edge_angular = np.zeros((quad.M, len(dx)+1))

    for m in range(quad.M):
        if quad.mu[m] < 0:
            cell_edge_angular[m, 0:-1] = psi[m].L
            cell_edge_angular[m, -1] = BC[m]
        if quad.mu[m] > 0:
            cell_edge_angular[m, 1:] = psi[m].R
            if lhs_reflective == True:
                cell_edge_angular[m, 0] = psi[quad.M-(m+1)].L[0]
            else:
                cell_edge_angular[m, 0] = BC[m]

    legendre_angular = np.array([0.75*(psi_m.R - psi_m.L)/dx for psi_m in psi])

    cell_average_scalar = np.transpose(cell_average_angular) @ quad.w
    cell_edge_scalar = np.transpose(cell_edge_angular) @ quad.w
    cell_average_scalar_old = phi_old.cell_average()
    legendre_scalar_old = 0.5*(phi_old.R - phi_old.L)
    legendre_scalar = 0.5*(phi.R - phi.L)
    cell_average_current = np.transpose(cell_average_angular) @ (quad.w * quad.mu)
    cell_edge_current = np.transpose(cell_edge_angular) @ (quad.w * quad.mu)
    legendre_current = np.transpose(legendre_angular) @ (quad.w * quad.mu)

    # calculate residuals
    scheme_balance = np.longdouble(cell_edge_current[1:]) - np.longdouble(cell_edge_current[:-1]) \
        + np.longdouble(xs.total)*np.longdouble(cell_average_scalar)*np.longdouble(dx) \
        - (np.longdouble(xs.scatter) + (np.longdouble(xs.fission) * np.longdouble(xs.multiplicity)))*np.longdouble(cell_average_scalar_old)*np.longdouble(dx) \
        - np.longdouble(q*dx)


    scheme_firstmoment = 3*np.longdouble(cell_edge_current[1:]) + 3*np.longdouble(cell_edge_current[:-1]) - 6*np.longdouble(cell_average_current) \
        + np.longdouble(xs.total)*np.longdouble(legendre_scalar)*np.longdouble(dx) \
        - (np.longdouble(xs.scatter) + np.longdouble(xs.fission)*np.longdouble(xs.multiplicity))*np.longdouble(legendre_scalar_old)*np.longdouble(dx)



    iterative_balance = np.longdouble(cell_edge_current[1:]) - np.longdouble(cell_edge_current[:-1]) \
        + (np.longdouble(xs.total)- np.longdouble(xs.scatter + (xs.fission * xs.multiplicity)))*np.longdouble(cell_average_scalar)*np.longdouble(dx) \
        - np.longdouble(q*dx)


    iterative_firstmoment = 3*np.longdouble(cell_edge_current[1:] + cell_edge_current[:-1] - 2*cell_average_current) \
        + np.longdouble(xs.total - xs.scatter - (xs.fission * xs.multiplicity))*np.longdouble(legendre_scalar)*np.longdouble(dx)


    scheme_residuals = Residuals(scheme_balance, scheme_firstmoment)
    iterative_residuals = Residuals(iterative_balance, iterative_firstmoment)


    return Solution(phi, psi), Convergence(delta, rho), scheme_residuals, iterative_residuals