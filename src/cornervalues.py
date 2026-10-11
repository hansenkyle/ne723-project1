class CornerValues:
    def __init__(self, L, R):
        self.L = L
        self.R = R

        if len(L) != len(R):
            raise ValueError(f"left and right CV vectors must have equal length (got {len(L)} and {len(R)})")

    def integrate(self, dx: np.array):
        if (len(dx) != len(self.L)):
            raise ValueError(f"len(dx) = {len(dx)} incompatible with len(L)=len(R)={len(self.L)}")