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