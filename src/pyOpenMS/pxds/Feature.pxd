from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from DataValue cimport *
from ConvexHull2D cimport *
from PeptideIdentification cimport *

from UniqueIdInterface cimport *

cdef extern from "<OpenMS/KERNEL/Feature.h>" namespace "OpenMS":

    cdef cppclass Feature(UniqueIdInterface):
        #
        # wrap-inherits:
        #    UniqueIdInterface
        #
        # wrap-doc:
        #   An LC-MS feature.
        #   -----
        #   The Feature class is used to describe the two-dimensional signal caused by an
        #   analyte. It can store a charge state and a list of peptide identifications
        #   (for peptides). The area occupied by the Feature in the LC-MS data set is
        #   represented by a list of convex hulls (one for each isotopic peak). There is
        #   also a convex hull for the entire Feature. The model description can store
        #   the parameters of a two-dimensional theoretical model of the underlying
        #   signal in LC-MS. Currently, non-peptide compounds are also represented as
        #   features.

        Feature() nogil except +
        Feature(Feature &) nogil except +

        void setMZ(double) nogil except +
        void setRT(double) nogil except +
        void setIntensity(double) nogil except +

        double getMZ() nogil except +
        double getRT() nogil except +
        double getIntensity() nogil except +

        float getQuality(Size index) nogil except +
        void setQuality(Size index, float q) nogil except +
        float getOverallQuality() nogil except +
        void setOverallQuality(float q) nogil except +

        float getWidth() nogil except +
        void setWidth(float q) nogil except +

        Int getCharge() nogil except +
        void setCharge(Int q) nogil except +

        libcpp_vector[Feature] getSubordinates() nogil except +
        void setSubordinates(libcpp_vector[Feature]) nogil except +

        bool encloses(double rt, double mz) nogil except +
        ConvexHull2D getConvexHull() nogil except +
        libcpp_vector[ConvexHull2D] getConvexHulls() nogil except +
        void setConvexHulls(libcpp_vector[ConvexHull2D]) nogil except +

        # returns a mutable reference to the PeptideIdentification vector
        libcpp_vector[PeptideIdentification] getPeptideIdentifications() nogil except +
        # sets the PeptideIdentification vector
        void setPeptideIdentifications(libcpp_vector[PeptideIdentification] & peptides) nogil except +

        bool operator==(Feature) nogil except +
        bool operator!=(Feature) nogil except +

        void getKeys(libcpp_vector[String] & keys) nogil except +
        void getKeys(libcpp_vector[unsigned int] & keys) nogil except + # wrap-as:getKeysAsIntegers
        DataValue getMetaValue(unsigned int) nogil except +
        DataValue getMetaValue(String) nogil except +
        void setMetaValue(unsigned int, DataValue) nogil except +
        void setMetaValue(String, DataValue) nogil except +
        bool metaValueExists(String) nogil except +
        bool metaValueExists(unsigned int) nogil except +
        void removeMetaValue(String) nogil except +
        void removeMetaValue(unsigned int) nogil except +
        void clearMetaInfo() nogil except +

