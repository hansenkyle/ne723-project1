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