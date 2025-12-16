



    def __repr__(self):
        """
        __repr__(self) -> str

        Return a string representation of the DRange1 object.
        """
        cdef double min_val = self.minX()
        cdef double max_val = self.maxX()
        return f"DRange1([{min_val}, {max_val}))"

    def __str__(self):
        """
        __str__(self) -> str

        Return a string representation of the DRange1 object.
        """
        return self.__repr__()
