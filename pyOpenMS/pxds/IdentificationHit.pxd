from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool

from Types cimport *
from MetaInfoInterface cimport *

cdef extern from "<OpenMS/METADATA/IdentificationHit.h>" namespace "OpenMS":

    cdef cppclass IdentificationHit(MetaInfoInterface):
        # wrap-inherits:
        #   MetaInfoInterface

        IdentificationHit()   nogil except +
        IdentificationHit(IdentificationHit) nogil except + # wrap-ignore

        # /// sets the identifier
        void setId(String id) nogil except +

        # /// returns the id
        String getId() nogil except +

        # /// sets the charge state of the peptide
        void setCharge(Int charge) nogil except +

        # /// returns the charge state
        Int getCharge() nogil except +

        # /// sets the calculated mass to charge ratio
        void setCalculatedMassToCharge(DoubleReal mz) nogil except +

        # /// returns the calculated mass to charge ratio
        DoubleReal getCalculatedMassToCharge() nogil except +

        # /// sets the experimental mass to charge ratio
        void setExperimentalMassToCharge(DoubleReal mz) nogil except +

        # /// returns the experimental mass to charge
        DoubleReal getExperimentalMassToCharge() nogil except +

        # /// sets the name
        void setName(String name) nogil except +

        # /// returns the name
        String getName() nogil except +

        # /// sets whether the peptide passed the threshold
        void setPassThreshold(bool) nogil except +

        # /// returns whether the peptide passed the threshold
        bool getPassThreshold() nogil except +

        # /// set the rank of the peptide
        void setRank(Int rank) nogil except +

        # /// returns the rank of the peptide
        Int getRank() nogil except +

        # COPY-PASTE from MetaInfoInterface
        void getKeys(libcpp_vector[String] & keys) nogil except +
        void getKeys(libcpp_vector[unsigned int] & keys) nogil except +
        DataValue getMetaValue(unsigned int) nogil except +
        DataValue getMetaValue(String) nogil except +
        void setMetaValue(unsigned int, DataValue) nogil except +
        void setMetaValue(String, DataValue) nogil except +
        bool metaValueExists(String) nogil except +
        bool metaValueExists(unsigned int) nogil except +
        void removeMetaValue(String) nogil except +
        void removeMetaValue(unsigned int) nogil except +
