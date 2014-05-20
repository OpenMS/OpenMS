from Types cimport *

cdef extern from "<OpenMS/FILTERING/DATAREDUCTION/SplineSpectrum.h>" namespace "OpenMS":
    
    cdef cppclass SplineSpectrum "OpenMS::SplineSpectrum":
        SplineSpectrum() nogil except +
        SplineSpectrum(SplineSpectrum) nogil except + #wrap-ignore
        double mzMin_
        double mzMax_
	vector<SplinePackage> packages_

