


    def __str__(self):
        """
        __str__(self: Peak1D) -> str
        
        Return a string representation of the Peak1D object.
        Delegates to __repr__ for consistency.
        """
        return self.__repr__()

    def __repr__(self):
        """
        __repr__(self: Peak1D) -> str
        
        Return a string representation of the Peak1D object.

        Returns key properties in a readable format:
        Peak1D(mz=445.678, intensity=100000.0)
        """
        cdef double mz = self.getMZ()
        cdef float intensity = self.getIntensity()

        return f"Peak1D(mz={mz:.4f}, intensity={intensity:.1f})"
