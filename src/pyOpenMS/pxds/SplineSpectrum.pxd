from Types cimport *
from MSSpectrum cimport *
from SplinePackage cimport *
from Peak1D cimport *

cdef extern from "<OpenMS/FILTERING/DATAREDUCTION/SplineSpectrum.h>" namespace "OpenMS":
    
    cdef cppclass SplineSpectrum "OpenMS::SplineSpectrum":

        SplineSpectrum(libcpp_vector[double] mz, libcpp_vector[double] intensity) nogil except +
        SplineSpectrum(libcpp_vector[double] mz, libcpp_vector[double] intensity, double scaling) nogil except +

        SplineSpectrum(MSSpectrum raw_spectrum) nogil except +
        SplineSpectrum(MSSpectrum raw_spectrum, double scaling) nogil except +

        double getMzMin() nogil except +

        double getMzMax() nogil except +

        int getSplineCount() nogil except +
        
        SplineSpectrum_Navigator getNavigator() nogil except +


cdef extern from "<OpenMS/FILTERING/DATAREDUCTION/SplineSpectrum.h>" namespace "OpenMS::SplineSpectrum":
    
    cdef cppclass SplineSpectrum_Navigator "OpenMS::SplineSpectrum::Navigator":
        
            SplineSpectrum_Navigator() nogil except + #wrap-ignore
            SplineSpectrum_Navigator(SplineSpectrum_Navigator) nogil except + #wrap-ignore
            
            SplineSpectrum_Navigator(libcpp_vector[SplinePackage]* packages, double mzMin, double mzMax)  nogil except +

            double eval(double mz) nogil except +

            double getNextMz(double mz) nogil except +
