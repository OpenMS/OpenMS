


    def __str__(self):
        """
        Return a string representation of the ConsensusFeature object.
        Delegates to __repr__ for consistency.
        """
        return self.__repr__()

    def __repr__(self):
        """
        Return a string representation of the ConsensusFeature object.

        Returns key properties in a readable format:
        ConsensusFeature(rt=1234.5, mz=445.678, intensity=100000.0, charge=2, num_features=5)
        """
        cdef double rt = self.getRT()
        cdef double mz = self.getMZ()
        cdef float intensity = self.getIntensity()
        cdef int charge = self.getCharge()
        cdef size_t num_features = self.size()

        parts = []
        parts.append(f"rt={rt:.2f}")
        parts.append(f"mz={mz:.4f}")
        parts.append(f"intensity={intensity:.1f}")
        if charge != 0:
            parts.append(f"charge={charge}")
        parts.append(f"num_features={num_features}")

        return f"ConsensusFeature({', '.join(parts)})"
