from numba.experimental import jitclass
from numba import float64

# 1. Define the specifications for your arrays.
# Format: type[colon_dimensions] -> float64[:] means 1D float array, int32[:,:] means 2D int array
spec = [
    ('total', float64[:]),     # 1D array of floats
    ('scatter', float64[:]),
    ('fission', float64[:]), 
    ('multiplicity', float64[:])
]

# class to hold all xs data

@jitclass(spec)
class CrossSection:

    # Constructor method to initialize instance attributes
    def __init__(self, total, scatter, fission, multiplicity):
        self.total = total 
        self.scatter = scatter
        self.fission = fission
        self.multiplicity = multiplicity

        if (len(total) != len(scatter)):
            raise ValueError(f"'total' and 'scatter' must be same length (got {len(total)} and {len(scatter)})")

        if (len(total) != len(fission)):
            raise ValueError(f"'total' and 'fission' must be same length (got {len(total)} and {len(fission)})")

        if (len(total) != len(multiplicity)):
            raise ValueError(f"'total' and 'multiplicity' must be same length (got {len(total)} and {len(multiplicity)})")