


    def __str__(self):
        """
        Return a string representation of the MRMFeature object.
        Delegates to __repr__ for consistency.
        """
        return self.__repr__()

    def __repr__(self):
        """
        Return a string representation of the MRMFeature object.

        Returns key properties in a readable format:
        MRMFeature(rt=1234.5, mz=445.678, intensity=100000.0, num_transitions=6)
        """
        cdef double rt = self.getRT()
        cdef double mz = self.getMZ()
        cdef float intensity = self.getIntensity()
        cdef int charge = self.getCharge()

        # Get number of transitions using Python method
        features = self.getFeatures()
        num_transitions = len(features)

        parts = []
        parts.append(f"rt={rt:.2f}")
        parts.append(f"mz={mz:.4f}")
        parts.append(f"intensity={intensity:.1f}")
        if charge != 0:
            parts.append(f"charge={charge}")
        parts.append(f"num_transitions={num_transitions}")

        return f"MRMFeature({', '.join(parts)})"
