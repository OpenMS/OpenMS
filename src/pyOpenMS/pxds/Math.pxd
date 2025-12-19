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

    # Template function getPPM instantiated for double
    double getPPM(double mz_obs, double mz_ref) except + nogil  # wrap-attach:Math
    
    # Template function getPPMAbs instantiated for double
    double getPPMAbs(double mz_obs, double mz_ref) except + nogil  # wrap-attach:Math
    
    # Template function ppmToMass instantiated for double
    double ppmToMass(double ppm, double mz_ref) except + nogil  # wrap-attach:Math
    
    # Template function ppmToMassAbs instantiated for double
    double ppmToMassAbs(double ppm, double mz_ref) except + nogil  # wrap-attach:Math
