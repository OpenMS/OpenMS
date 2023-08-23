from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool

from Types cimport *
from MetaInfoInterface cimport *

cdef extern from "<OpenMS/METADATA/IdentificationHit.h>" namespace "OpenMS":

    cdef cppclass IdentificationHit(MetaInfoInterface):
        # wrap-inherits:
        #  MetaInfoInterface

        IdentificationHit()   except + nogil  # wrap-doc:Represents a object which can store the information of an analysisXML instance
        IdentificationHit(IdentificationHit &) except + nogil 

        void setId(String id) except + nogil  # wrap-doc:Sets the identifier

        String getId() except + nogil  # wrap-doc:Returns the id

        void setCharge(Int charge) except + nogil  # wrap-doc:Sets the charge state of the peptide

        Int getCharge() except + nogil  # wrap-doc:Returns the charge state

        void setCalculatedMassToCharge(double mz) except + nogil  # wrap-doc:Sets the calculated mass to charge ratio

        double getCalculatedMassToCharge() except + nogil  # wrap-doc:Returns the calculated mass to charge ratio

        void setExperimentalMassToCharge(double mz) except + nogil  # wrap-doc:Sets the experimental mass to charge ratio

        double getExperimentalMassToCharge() except + nogil  # wrap-doc:Returns the experimental mass to charge

        void setName(String name) except + nogil  # wrap-doc:Sets the name

        String getName() except + nogil  # wrap-doc:Returns the name

        void setPassThreshold(bool) except + nogil  # wrap-doc:Sets whether the peptide passed the threshold

        bool getPassThreshold() except + nogil  # wrap-doc:Returns whether the peptide passed the threshold

        void setRank(Int rank) except + nogil  # wrap-doc:Sets the rank of the peptide

        Int getRank() except + nogil  # wrap-doc:Returns the rank of the peptide
