from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool

from Types cimport *
from MetaInfoInterface cimport *

cdef extern from "<OpenMS/METADATA/IdentificationHit.h>" namespace "OpenMS":

    cdef cppclass IdentificationHit(MetaInfoInterface):
        # wrap-inherits:
        #   MetaInfoInterface

        IdentificationHit()   nogil except + # wrap-doc:Represents a object which can store the information of an analysisXML instance
        IdentificationHit(IdentificationHit &) nogil except +

        void setId(String id) nogil except + # wrap-doc:Sets the identifier

        String getId() nogil except + # wrap-doc:Returns the id

        void setCharge(Int charge) nogil except + # wrap-doc:Sets the charge state of the peptide

        Int getCharge() nogil except + # wrap-doc:Returns the charge state

        void setCalculatedMassToCharge(double mz) nogil except + # wrap-doc:Sets the calculated mass to charge ratio

        double getCalculatedMassToCharge() nogil except + # wrap-doc:Returns the calculated mass to charge ratio

        void setExperimentalMassToCharge(double mz) nogil except + # wrap-doc:Sets the experimental mass to charge ratio

        double getExperimentalMassToCharge() nogil except + # wrap-doc:Returns the experimental mass to charge

        void setName(String name) nogil except + # wrap-doc:Sets the name

        String getName() nogil except + # wrap-doc:Returns the name

        void setPassThreshold(bool) nogil except + # wrap-doc:Sets whether the peptide passed the threshold

        bool getPassThreshold() nogil except + # wrap-doc:Returns whether the peptide passed the threshold

        void setRank(Int rank) nogil except + # wrap-doc:Sets the rank of the peptide

        Int getRank() nogil except + # wrap-doc:Returns the rank of the peptide
