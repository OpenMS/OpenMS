from Types cimport *

cdef extern from "<OpenMS/PROCESSING/MISC/SplinePackage.h>" namespace "OpenMS":
    
    cdef cppclass SplinePackage "OpenMS::SplinePackage":
        # wrap-doc:
        #  fundamental data structure for SplineInterpolatedPeaks * * In many
        #  cases, data points in MS spectra (or chromatograms) are not
        #  equidistant in m/z (or RT) * but consist of packages of data points
        #  separated by wide m/z (or RT) ranges with zero intensity. *
        #  SplinePackage contains the spline fit of a single set of such data
        #  points. * * @see SplineInterpolatedPeaks

        SplinePackage(libcpp_vector[double] pos, libcpp_vector[double] intensity) except + nogil 
        SplinePackage(SplinePackage &) except + nogil  # compiler

        double getPosMin() except + nogil  # wrap-doc:Returns the minimum position for which the spline fit is valid
        double getPosMax() except + nogil  # wrap-doc:Returns the maximum position for which the spline fit is valid
        double getPosStepWidth() except + nogil  # wrap-doc:Returns a sensible position step width for the package
        bool isInPackage(double pos) except + nogil  # wrap-doc:Returns true if position in [posMin:posMax] interval else false
        double eval(double pos) except + nogil  # wrap-doc:Returns interpolated intensity position `pos`
