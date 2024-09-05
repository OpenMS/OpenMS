from Types cimport *
from MSSpectrum cimport *
from MSChromatogram cimport *
from SplinePackage cimport *
from Peak1D cimport *

cdef extern from "<OpenMS/PROCESSING/MISC/SplineInterpolatedPeaks.h>" namespace "OpenMS":
    
    cdef cppclass SplineInterpolatedPeaks "OpenMS::SplineInterpolatedPeaks":

        # private
        SplineInterpolatedPeaks() except + nogil  # wrap-ignore

        SplineInterpolatedPeaks(libcpp_vector[double] mz, libcpp_vector[double] intensity) except + nogil 

        SplineInterpolatedPeaks(MSSpectrum raw_spectrum) except + nogil 

        SplineInterpolatedPeaks(MSChromatogram raw_chromatogram) except + nogil 

        SplineInterpolatedPeaks(SplineInterpolatedPeaks &) except + nogil  # compiler

        double getPosMin() except + nogil 

        double getPosMax() except + nogil 

        int size() except + nogil 
        
        SplineSpectrum_Navigator getNavigator(double scaling) except + nogil 


cdef extern from "<OpenMS/PROCESSING/MISC/SplineInterpolatedPeaks.h>" namespace "OpenMS::SplineInterpolatedPeaks":
    
    cdef cppclass SplineSpectrum_Navigator "OpenMS::SplineInterpolatedPeaks::Navigator":
        
            SplineSpectrum_Navigator() except + nogil 
            SplineSpectrum_Navigator(SplineSpectrum_Navigator) except + nogil  # compiler
            
            SplineSpectrum_Navigator(libcpp_vector[SplinePackage]* packages, double posMax, double scaling)  except + nogil 

            double eval(double pos) except + nogil 

            double getNextPos(double pos) except + nogil 
