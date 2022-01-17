from Types cimport *
from MSSpectrum cimport *
from String cimport *

# typedef std::map<String, std::vector<std::pair<double, double> > > MarkerIonsType;

ctypedef libcpp_map[String, libcpp_vector[ libcpp_pair[double, double] ] ] MarkerIonsType

cdef extern from "<OpenMS/ANALYSIS/RNPXL/RNPxlMarkerIonExtractor.h>" namespace "OpenMS":
    
    cdef cppclass RNPxlMarkerIonExtractor "OpenMS::RNPxlMarkerIonExtractor":
        RNPxlMarkerIonExtractor() nogil except + # compiler
        RNPxlMarkerIonExtractor(RNPxlMarkerIonExtractor &) nogil except + # compiler
        # MarkerIonsType extractMarkerIons(MSSpectrum & s, double marker_tolerance) nogil except +
        libcpp_map[String, libcpp_vector[ libcpp_pair[double, double] ] ] extractMarkerIons(MSSpectrum & s,
                                                                                            double marker_tolerance) nogil except +
