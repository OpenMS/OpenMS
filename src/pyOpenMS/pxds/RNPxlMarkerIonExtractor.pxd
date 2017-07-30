from Types cimport *
from MSSpectrum cimport *

# typedef std::map<String, std::vector<std::pair<double, double> > > MarkerIonsType;
cdef extern from "<OpenMS/ANALYSIS/RNPXL/RNPxlMarkerIonExtractor.h>" namespace "OpenMS":
    
    cdef cppclass RNPxlMarkerIonExtractor "OpenMS::RNPxlMarkerIonExtractor":
        RNPxlMarkerIonExtractor() nogil except + 
        RNPxlMarkerIonExtractor(RNPxlMarkerIonExtractor) nogil except + #wrap-ignore
        # MarkerIonsType extractMarkerIons(MSSpectrum[Peak1D] & s, double marker_tolerance) nogil except +

