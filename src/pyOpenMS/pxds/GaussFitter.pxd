from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool
from DPosition cimport *
from String cimport *

cdef extern from "<OpenMS/MATH/STATISTICS/GaussFitter.h>" namespace "OpenMS::Math":

    cdef cppclass GaussFitter:

        GaussFitter() nogil except + # wrap-doc:Implements a fitter for Gaussian functions
        # private
        GaussFitter(GaussFitter) nogil except + # wrap-ignore

        # sets the initial parameters used by the fit method as inital guess for the gaussian
        void setInitialParameters(GaussFitResult & result) nogil except + # wrap-doc:Sets the initial parameters used by the fit method as initial guess for the Gaussian

        #   @brief Fits a gaussian distribution to the given data points
        #   @param points the data points used for the gaussian fitting
        #   @exception Exception::UnableToFit is thrown if fitting cannot be performed
        GaussFitResult fit(libcpp_vector[DPosition2] points) nogil except + # wrap-doc:Fits a Gaussian distribution to the given data points

cdef extern from "<OpenMS/MATH/STATISTICS/GaussFitter.h>" namespace "OpenMS::Math::GaussFitter":

    cdef cppclass GaussFitResult:

        GaussFitResult() nogil except +
        GaussFitResult(double, double, double) nogil except +
        GaussFitResult(GaussFitResult &) nogil except + # compiler

        double eval(double) nogil except +

        # parameter A of gaussian distribution (amplitude)
        double A
        # parameter x0 of gaussian distribution (left/right shift)
        double x0
        # parameter sigma of gaussian distribution (width)
        double sigma
