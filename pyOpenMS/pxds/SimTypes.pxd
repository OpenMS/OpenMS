from libcpp.map cimport map as libcpp_map
from libcpp.pair cimport pair as libcpp_pair
from libcpp cimport bool
from Types cimport *
from ProgressLogger cimport *
from DefaultParamHandler cimport *
from FASTAFile cimport *
from MetaInfoInterface cimport *

cdef extern from "<OpenMS/SIMULATION/SimTypes.h>" namespace "OpenMS":

    cdef cppclass SimRandomNumberGenerator:

        SimRandomNumberGenerator() nogil except +
        SimRandomNumberGeneratorSimTypes(SimRandomNumberGeneratorSimTypes) nogil except + # wrap-ignore
        void initialize(bool biological_random, bool technical_random) nogil except +
        

    # /// Container for FASTAEntry & abundance information
    # class OPENMS_DLLAPI SampleProteins : public std::vector<std::pair<FASTAFile::FASTAEntry, MetaInfoInterface> > { }; 
    cdef cppclass SampleProteins:

        SampleProteins() nogil except +
        SampleProteins(SampleProteins) nogil except + # wrap-ignore
        void push_back(libcpp_pair[FASTAEntry, MetaInfoInterface]) nogil except +

    # /// Container for multiple channels of SampleProteins
    # class OPENMS_DLLAPI SampleChannels : public std::vector<SampleProteins> { }; 
    cdef cppclass SampleChannels:

        SampleChannels() nogil except +
        SampleChannels(SampleChannels) nogil except + # wrap-ignore
        void push_back(SampleProteins p) nogil except +

