from Types cimport *
from MSSpectrum cimport *
from MSChromatogram cimport *
from SplinePackage cimport *
from Peak1D cimport *

cdef extern from "<OpenMS/PROCESSING/MISC/SplineInterpolatedPeaks.h>" namespace "OpenMS":
    
    cdef cppclass SplineInterpolatedPeaks "OpenMS::SplineInterpolatedPeaks":
        # wrap-doc:
        #  Data structure for spline interpolation of MS1 spectra and
        #  chromatograms * * The data structure consists of a set of splines,
        #  each interpolating the MS1 spectrum (or chromatogram) in a * certain
        #  m/z (or RT) range. Between these splines no raw data points exist and
        #  the intensity is identical to zero. * * A spline on non-equi-distant
        #  input data is not well supported in regions without data points.
        #  Hence, a spline tends to * swing wildly in these regions and cannot be
        #  used for reliable interpolation. We assume that in m/z (or RT) regions
        #  * without data points, the spectrum (or chromatogram) is identical to
        #  zero. * * @see SplinePackage * @see MSSpectrum * @see MSChromatogram

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
