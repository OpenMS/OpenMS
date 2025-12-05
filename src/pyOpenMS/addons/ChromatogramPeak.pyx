


    def __str__(self):
        """
        Return a string representation of the ChromatogramPeak object.
        Delegates to __repr__ for consistency.
        """
        return self.__repr__()

    def __repr__(self):
        """
        Return a string representation of the ChromatogramPeak object.

        Returns key properties in a readable format:
        ChromatogramPeak(rt=1234.5, intensity=100000.0)
        """
        cdef double rt = self.getRT()
        cdef double intensity = self.getIntensity()

        return f"ChromatogramPeak(rt={rt:.2f}, intensity={intensity:.1f})"
