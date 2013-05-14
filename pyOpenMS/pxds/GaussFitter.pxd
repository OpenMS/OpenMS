from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool
from DPosition cimport *
from String cimport *

cdef extern from "<OpenMS/MATH/STATISTICS/GaussFitter.h>" namespace "OpenMS::Math":

    cdef cppclass GaussFitter:

        GaussFitter() nogil except +
        GaussFitter(GaussFitter) nogil except +   # wrap-ignore

        # sets the initial parameters used by the fit method as inital guess for the gaussian
        void setInitialParameters(GaussFitResult & result) nogil except +

        #   @brief Fits a gaussian distribution to the given data points
        #   @param points the data points used for the gaussian fitting
        #   @exception Exception::UnableToFit is thrown if fitting cannot be performed
        GaussFitResult fit(libcpp_vector[DPosition2] points) nogil except +

        # return the gnuplot formula of the gaussian
        String getGnuplotFormula() nogil except +

cdef extern from "<OpenMS/MATH/STATISTICS/GaussFitter.h>" namespace "OpenMS::Math::GaussFitter":

    cdef cppclass GaussFitResult:

        GaussFitResult() nogil except +
        GaussFitResult(GaussFitResult) nogil except +   # wrap-ignore

        # parameter A of gaussian distribution (amplitude)
        double A
        # parameter x0 of gaussian distribution (left/right shift)
        double x0
        # parameter sigma of gaussian distribution (width)
        double sigma
