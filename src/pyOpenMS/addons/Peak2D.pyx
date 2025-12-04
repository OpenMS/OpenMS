


    def __str__(self):
        """
        Return a string representation of the Peak2D object.
        Delegates to __repr__ for consistency.
        """
        return self.__repr__()

    def __repr__(self):
        """
        Return a string representation of the Peak2D object.

        Returns key properties in a readable format:
        Peak2D(rt=1234.5, mz=445.678, intensity=100000.0)
        """
        cdef double rt = self.getRT()
        cdef double mz = self.getMZ()
        cdef float intensity = self.getIntensity()

        return f"Peak2D(rt={rt:.2f}, mz={mz:.4f}, intensity={intensity:.1f})"
