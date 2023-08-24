from Types cimport *
from MSSpectrum cimport *
from String cimport *
from libcpp.map cimport map as libcpp_map

# typedef std::map<String, std::vector<std::pair<double, double> > > MarkerIonsType;

ctypedef libcpp_map[String, libcpp_vector[ libcpp_pair[double, double] ] ] MarkerIonsType

cdef extern from "<OpenMS/ANALYSIS/RNPXL/RNPxlMarkerIonExtractor.h>" namespace "OpenMS":
    
    cdef cppclass RNPxlMarkerIonExtractor "OpenMS::RNPxlMarkerIonExtractor":
        RNPxlMarkerIonExtractor() except + nogil  # compiler
        RNPxlMarkerIonExtractor(RNPxlMarkerIonExtractor &) except + nogil  # compiler
        # MarkerIonsType extractMarkerIons(MSSpectrum & s, double marker_tolerance) except + nogil 
        libcpp_map[String, libcpp_vector[ libcpp_pair[double, double] ] ] extractMarkerIons(MSSpectrum & s,
                                                                                            double marker_tolerance) except + nogil 
