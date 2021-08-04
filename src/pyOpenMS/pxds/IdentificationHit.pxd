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

        # /// sets the identifier
        void setId(String id) nogil except + # wrap-doc:Sets the identifier

        # /// returns the id
        String getId() nogil except + # wrap-doc:Returns the id

        # /// sets the charge state of the peptide
        void setCharge(Int charge) nogil except + # wrap-doc:Sets the charge state of the peptide

        # /// returns the charge state
        Int getCharge() nogil except + # wrap-doc:Returns the charge state

        # /// sets the calculated mass to charge ratio
        void setCalculatedMassToCharge(double mz) nogil except + # wrap-doc:Sets the calculated mass to charge ratio

        # /// returns the calculated mass to charge ratio
        double getCalculatedMassToCharge() nogil except + # wrap-doc:Returns the calculated mass to charge ratio

        # /// sets the experimental mass to charge ratio
        void setExperimentalMassToCharge(double mz) nogil except + # wrap-doc:Sets the experimental mass to charge ratio

        # /// returns the experimental mass to charge
        double getExperimentalMassToCharge() nogil except + # wrap-doc:Returns the experimental mass to charge

        # /// sets the name
        void setName(String name) nogil except + # wrap-doc:Sets the name

        # /// returns the name
        String getName() nogil except + # wrap-doc:Returns the name

        # /// sets whether the peptide passed the threshold
        void setPassThreshold(bool) nogil except + # wrap-doc:Sets whether the peptide passed the threshold

        # /// returns whether the peptide passed the threshold
        bool getPassThreshold() nogil except + # wrap-doc:Returns whether the peptide passed the threshold

        # /// set the rank of the peptide
        void setRank(Int rank) nogil except + # wrap-doc:Set the rank of the peptide

        # /// returns the rank of the peptide
        Int getRank() nogil except + # wrap-doc:Returns the rank of the peptide

