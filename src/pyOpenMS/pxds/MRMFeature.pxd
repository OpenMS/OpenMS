from libcpp.vector cimport vector as libcpp_vector
from libcpp.map cimport map as libcpp_map
from String cimport *
from Feature cimport *

cdef extern from "<OpenMS/KERNEL/MRMFeature.h>" namespace "OpenMS":

    cdef cppclass MRMFeature:

        MRMFeature() nogil except +
        MRMFeature(MRMFeature &) nogil except +

        # TODO STL map with wrapped key
        # libcpp_map[String, double] getScores() nogil except +
        double getScore(String name) nogil except +
        void addScore(String name, double score) nogil except +

        Feature getFeature(String key) nogil except +
        void addFeature(Feature & f, String key) nogil except +
        libcpp_vector[Feature] getFeatures() nogil except +
        void getFeatureIDs(libcpp_vector[String] & result) nogil except +

        Feature getPrecursorFeature(String key) nogil except +
        void addPrecursorFeature(Feature & f, String key) nogil except +
        void getPrecursorFeatureIDs(libcpp_vector[String] & result) nogil except +

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

        bool operator==(MRMFeature) nogil except +
        bool operator!=(MRMFeature) nogil except +

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

