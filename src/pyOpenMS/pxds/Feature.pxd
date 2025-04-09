from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from ConvexHull2D cimport *
from PeptideIdentification cimport *
from BaseFeature cimport *

from UniqueIdInterface cimport *

cdef extern from "<OpenMS/KERNEL/Feature.h>" namespace "OpenMS":

    cdef cppclass Feature(UniqueIdInterface, RichPeak2D):
        #
        # wrap-inherits:
        #   UniqueIdInterface
        #   RichPeak2D
        #
        # wrap-doc:
        #  An LC-MS feature
        #  
        #  The Feature class is used to describe the two-dimensional signal caused by an
        #  analyte. It can store a charge state and a list of peptide identifications
        #  (for peptides). The area occupied by the Feature in the LC-MS data set is
        #  represented by a list of convex hulls (one for each isotopic peak). There is
        #  also a convex hull for the entire Feature. The model description can store
        #  the parameters of a two-dimensional theoretical model of the underlying
        #  signal in LC-MS. Currently, non-peptide compounds are also represented as
        #  features

        Feature() except + nogil 
        Feature(Feature &) except + nogil 

        float getQuality(Size index) except + nogil  # wrap-doc:Returns the quality in dimension c
        void setQuality(Size index, float q) except + nogil  # wrap-doc:Sets the quality in dimension c
        float getOverallQuality() except + nogil  # wrap-doc:Model and quality methods
        void setOverallQuality(float q) except + nogil  # wrap-doc:Sets the overall quality

        libcpp_vector[Feature] getSubordinates() except + nogil  # wrap-doc:Returns the subordinate features
        void setSubordinates(libcpp_vector[Feature]) except + nogil  # wrap-doc:Returns the subordinate features

        bool encloses(double rt, double mz) except + nogil  
            # wrap-doc:
            #  Returns if the mass trace convex hulls of the feature enclose the position specified by `rt` and `mz`
            #  
            #  
            #  :param rt: Sequence to digest
            #  :param mz: Digestion products
            
        ConvexHull2D getConvexHull() except + nogil 
        libcpp_vector[ConvexHull2D] getConvexHulls() except + nogil 
        void setConvexHulls(libcpp_vector[ConvexHull2D]) except + nogil 

        bool operator==(Feature) except + nogil 
        bool operator!=(Feature) except + nogil 

        # from BaseFeature

        # float getQuality()  except + nogil 
        # void setQuality(float q) except + nogil 

        float getWidth() except + nogil 
        void setWidth(float q) except + nogil 

        Int getCharge() except + nogil 
        void setCharge(Int q) except + nogil 
        AnnotationState getAnnotationState() except + nogil 

   
        libcpp_vector[PeptideIdentification] getPeptideIdentifications() except + nogil  # wrap-doc:Returns a reference to the PeptideIdentification vector
        
        void setPeptideIdentifications(libcpp_vector[PeptideIdentification] & peptides) except + nogil  # wrap-doc:Sets the PeptideIdentification vector

