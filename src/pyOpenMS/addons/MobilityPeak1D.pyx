


    def __str__(self):
        """
        __str__(self: MobilityPeak1D) -> str
        
        __str__(self: MobilityPeak1D) -> str
        
        Return a string representation of the MobilityPeak1D object.
        Delegates to __repr__ for consistency.
        """
        return self.__repr__()

    def __repr__(self):
        """
        __repr__(self: MobilityPeak1D) -> str
        
        __repr__(self: MobilityPeak1D) -> str
        
        Return a string representation of the MobilityPeak1D object.

        Returns key properties in a readable format:
        MobilityPeak1D(mobility=1.234, intensity=100000.0)
        """
        cdef double mobility = self.getMobility()
        cdef float intensity = self.getIntensity()

        return f"MobilityPeak1D(mobility={mobility:.4f}, intensity={intensity:.1f})"
