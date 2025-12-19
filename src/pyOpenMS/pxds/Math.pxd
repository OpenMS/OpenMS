cdef extern from "<OpenMS/MATH/MathFunctions.h>" namespace "OpenMS::Math":

    cdef cppclass Math:
        # wrap-doc:
        #  Math namespace wrapper for mathematical utility functions
        
        # Dummy class to attach Math namespace functions
        # This class should not be instantiated directly
        Math() except + nogil  # wrap-ignore
        Math(Math &) except + nogil  # wrap-ignore

    # Template function getPPM instantiated for double
    double getPPM(double mz_obs, double mz_ref) except + nogil  # wrap-attach:Math wrap-as:getPPM
    
    # Template function getPPMAbs instantiated for double
    double getPPMAbs(double mz_obs, double mz_ref) except + nogil  # wrap-attach:Math wrap-as:getPPMAbs
    
    # Template function ppmToMass instantiated for double
    double ppmToMass(double ppm, double mz_ref) except + nogil  # wrap-attach:Math wrap-as:ppmToMass
    
    # Template function ppmToMassAbs instantiated for double
    double ppmToMassAbs(double ppm, double mz_ref) except + nogil  # wrap-attach:Math wrap-as:ppmToMassAbs
