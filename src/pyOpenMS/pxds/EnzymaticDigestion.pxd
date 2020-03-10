from Types cimport *
from DigestionEnzyme cimport *
from StringView cimport *
from libcpp cimport bool

cdef extern from "<OpenMS/CHEMISTRY/EnzymaticDigestion.h>" namespace "OpenMS":

    cdef cppclass EnzymaticDigestion "OpenMS::EnzymaticDigestion":
        EnzymaticDigestion(EnzymaticDigestion) nogil except +
        EnzymaticDigestion() nogil except +

        # const String NamesOfSpecificity[SIZE_OF_SPECIFICITY]
        # const String UnspecificCleavage

        Size getMissedCleavages() nogil except +
        void setMissedCleavages(Size missed_cleavages) nogil except +
        String getEnzymeName() nogil except +
        void setEnzyme(DigestionEnzyme* enzyme) nogil except +
        Specificity getSpecificity() nogil except +
        void setSpecificity(Specificity spec) nogil except +
        Specificity getSpecificityByName(const String& name) nogil except +

        Size digestUnmodified(StringView sequence,
                              libcpp_vector[ StringView ]& output,
                              Size min_length,
                              Size max_length) nogil except +

        bool isValidProduct(const String& sequence,
                            int pos, int length,
                            bool ignore_missed_cleavages) nogil except +

        # bool filterByMissedCleavages(const String& sequence,
        #                              std::function<bool(Int)> filter) nogil except +


cdef extern from "<OpenMS/CHEMISTRY/EnzymaticDigestion.h>" namespace "OpenMS::EnzymaticDigestion":
    cdef enum Specificity "OpenMS::EnzymaticDigestion::Specificity":
        #wrap-attach:
        #    EnzymaticDigestion
        CROSS
        MONO
        LOOP
        NUMBER_OF_CROSS_LINK_TYPES
