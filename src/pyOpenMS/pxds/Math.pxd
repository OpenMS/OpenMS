cdef extern from "<OpenMS/MATH/MathFunctions.h>" namespace "OpenMS":

    cdef cppclass Math:
        # wrap-doc:
        #  Math namespace wrapper for mathematical utility functions
        
        # Dummy class to attach Math namespace functions
        # This class should not be instantiated directly
        Math() except + nogil  # wrap-ignore
        Math(Math &) except + nogil  # wrap-ignore


cdef extern from "<OpenMS/MATH/MathFunctions.h>" namespace "OpenMS::Math":

    double getPPM(double mz_obs, double mz_ref) except + nogil  # wrap-attach:Math
    double getPPMAbs(double mz_obs, double mz_ref) except + nogil  # wrap-attach:Math
    double ppmToMass(double ppm, double mz_ref) except + nogil  # wrap-attach:Math
    double ppmToMassAbs(double ppm, double mz_ref) except + nogil  # wrap-attach:Math
