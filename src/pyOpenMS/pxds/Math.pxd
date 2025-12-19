cdef extern from "<OpenMS/MATH/MathFunctions.h>" namespace "OpenMS::Math":

    cdef cppclass Math:
        # wrap-doc:
        #  Math namespace wrapper for mathematical utility functions.
        #  
        #  Provides static methods for PPM (parts-per-million) calculations commonly
        #  used in mass spectrometry for comparing observed and theoretical m/z values.
        #  
        #  Example usage in Python:
        #    >>> ppm = pyopenms.Math.getPPM(1000.001, 1000.0)  # Returns 1.0 ppm
        #    >>> mass_diff = pyopenms.Math.ppmToMass(5.0, 1000.0)  # Returns 0.005 Da
        
        # Dummy class to attach Math namespace functions
        # This class should not be instantiated directly
        Math() except + nogil  # wrap-ignore
        Math(Math &) except + nogil  # wrap-ignore

    double getPPM(double mz_obs, double mz_ref) except + nogil  # wrap-attach:Math wrap-doc:Compute parts-per-million difference between two m/z values. The returned ppm value can be positive (mz_obs > mz_ref) or negative (mz_obs < mz_ref)
    
    double getPPMAbs(double mz_obs, double mz_ref) except + nogil  # wrap-attach:Math wrap-doc:Compute absolute parts-per-million difference between two m/z values. The returned ppm value is always >= 0
    
    double ppmToMass(double ppm, double mz_ref) except + nogil  # wrap-attach:Math wrap-doc:Compute the mass difference in Dalton [Da] given a ppm value and a reference m/z. The returned mass difference can be positive (ppm > 0) or negative (ppm < 0)
    
    double ppmToMassAbs(double ppm, double mz_ref) except + nogil  # wrap-attach:Math wrap-doc:Compute the absolute mass difference in Dalton [Da] given a ppm value and a reference m/z. The returned mass difference is always positive
