


    def __str__(self):
        """
        __str__(self: ChromatogramPeak) -> str
        
        Return a string representation of the ChromatogramPeak object.
        Delegates to __repr__ for consistency.
        """
        return self.__repr__()

    def __repr__(self):
        """
        __repr__(self: ChromatogramPeak) -> str
        
        Return a string representation of the ChromatogramPeak object.

        Returns key properties in a readable format:
        ChromatogramPeak(rt=1234.5, intensity=100000.0)
        """
        cdef double rt = self.getRT()
        cdef double intensity = self.getIntensity()

        return f"ChromatogramPeak(rt={rt:.2f}, intensity={intensity:.1f})"
