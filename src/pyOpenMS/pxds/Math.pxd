cdef extern from "<OpenMS/MATH/MathFunctions.h>" namespace "OpenMS::Math":

    cdef cppclass Math:
        # wrap-doc:
        #  Math namespace wrapper providing mathematical utility functions.
        #
        #  This class provides static methods for various mathematical operations
        #  commonly used in mass spectrometry data analysis, including PPM calculations,
        #  statistical functions, and other numerical utilities.
        #
        #  Example usage:
        #    >>> ppm = pyopenms.Math.getPPM(1000.001, 1000.0)  # Returns 1.0 ppm
        #    >>> mass_diff = pyopenms.Math.ppmToMass(5.0, 1000.0)  # Returns 0.005 Da

        # Dummy class to attach Math namespace functions
        # This class should not be instantiated directly
        Math() except + nogil  # wrap-ignore
        Math(Math &) except + nogil  # wrap-ignore

        @staticmethod
        double getPPM(double mz_obs, double mz_ref) except + nogil  # wrap-doc:Compute parts-per-million difference between two m/z values. The returned ppm value can be positive (mz_obs > mz_ref) or negative (mz_obs < mz_ref)

        @staticmethod
        double getPPMAbs(double mz_obs, double mz_ref) except + nogil  # wrap-doc:Compute absolute parts-per-million difference between two m/z values. The returned ppm value is always >= 0

        @staticmethod
        double ppmToMass(double ppm, double mz_ref) except + nogil  # wrap-doc:Compute the mass difference in Dalton [Da] given a ppm value and a reference m/z. The returned mass difference can be positive (ppm > 0) or negative (ppm < 0)

        @staticmethod
        double ppmToMassAbs(double ppm, double mz_ref) except + nogil  # wrap-doc:Compute the absolute mass difference in Dalton [Da] given a ppm value and a reference m/z. The returned mass difference is always positive
