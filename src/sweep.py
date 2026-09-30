import numpy as np

def ldSolve(
        mu : float, 
        dx : float,
        psi_in : float, # upwind angular flux (BC or neighboring cell)
        xs_t : float,
        xs_s : float,
        xs_f : float,
        nu_f : float,
        q    : float,
        phi  : tuple   # expected in (L, R) format. this function re-orders as necessary
):
    s_eff = (xs_s + nu_f*xs_f)
    m = abs(mu)
    if mu<0:
        phi_l = phi[1]
        phi_r = phi[0]
    else:
        phi_l = phi[0]
        phi_r = phi[1]
    psi_l = (
                dx**2 * (q+s_eff*phi_l) 
                + (12*m**2 * psi_in) 
                + (dx*m*(s_eff*(phi_l - phi_r) + 8*xs_t*psi_in))
             )/(
                2*(6*m**2 + 4*dx*m*xs_t + dx**2 * xs_t**2)
             )
    psi_r = (
                0.5*dx*(s_eff*(phi_l/6 + phi_r/3)+0.5*q) + (0.5*m - (dx*xs_t)/6)*psi_l
            )/(
                0.5*m + (xs_t*dx)/3
            )
    if mu>0:
        return psi_l, psi_r
    else:
        return psi_r, psi_l
