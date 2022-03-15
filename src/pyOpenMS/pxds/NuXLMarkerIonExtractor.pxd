from Types cimport *
from MSSpectrum cimport *
from String cimport *
from libcpp.map cimport map as libcpp_map

# typedef std::map<String, std::vector<std::pair<double, double> > > MarkerIonsType;

ctypedef libcpp_map[String, libcpp_vector[ libcpp_pair[double, double] ] ] MarkerIonsType

cdef extern from "<OpenMS/ANALYSIS/NUXL/NuXLMarkerIonExtractor.h>" namespace "OpenMS":
    
    cdef cppclass NuXLMarkerIonExtractor "OpenMS::NuXLMarkerIonExtractor":
        NuXLMarkerIonExtractor() nogil except + 
        NuXLMarkerIonExtractor(NuXLMarkerIonExtractor &) nogil except + #wrap-ignore
        # MarkerIonsType extractMarkerIons(MSSpectrum & s, double marker_tolerance) nogil except +
        libcpp_map[String, libcpp_vector[ libcpp_pair[double, double] ] ] extractMarkerIons(MSSpectrum & s,
                                                                                            double marker_tolerance) nogil except +
