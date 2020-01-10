from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from libcpp.string cimport string as libcpp_string
from StringList cimport *
from MSSpectrum cimport *


cdef extern from "<OpenMS/CHEMISTRY/Tagger.h>" namespace "OpenMS":

    cdef cppclass Tagger:

        Tagger(Tagger) nogil except +

        Tagger(size_t min_tag_length,
               double ppm,
               size_t max_tag_length,
               size_t min_charge,
               size_t max_charge,
               const StringList& fixed_mods,
               const StringList& var_mods) nogil except +

        void getTag(const libcpp_vector[ double ]& mzs,
                    libcpp_vector[ libcpp_string ]& tags) nogil except +

        void getTag(const MSSpectrum& spec,
                    libcpp_vector[ libcpp_string ]& tags) nogil except +

        void setMaxCharge(size_t max_charge) nogil except +
