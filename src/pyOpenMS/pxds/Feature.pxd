from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from ConvexHull2D cimport *
from PeptideIdentification cimport *
from PeptideIdentificationList cimport *
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
        #  An LC-MS feature representing a detected analyte signal
        #  
        #  The Feature class represents a two-dimensional (RT and m/z) signal from an analyte
        #  in LC-MS data. It is one of the core data structures in OpenMS for representing
        #  detected peaks or compounds.
        #  
        #  A Feature stores:
        #  - Position: retention time (RT) and mass-to-charge ratio (m/z)
        #  - Intensity: the signal strength (typically total ion count)
        #  - Quality metrics: scores indicating detection confidence
        #  - Charge state: the charge of the ion
        #  - Convex hulls: the 2D area occupied by the feature in RT-m/z space
        #  - Peptide identifications: for identified peptides (optional)
        #  - Subordinate features: for isotopic peaks or related signals
        #  
        #  By convention, the feature's position is at the maximum of the elution profile
        #  (RT dimension) and at the monoisotopic peak (m/z dimension).
        #  
        #  Example usage:
        #    >>> feature = oms.Feature()
        #    >>> feature.setRT(1234.5)  # Set retention time in seconds
        #    >>> feature.setMZ(445.678)  # Set m/z value
        #    >>> feature.setIntensity(100000.0)  # Set intensity
        #    >>> feature.setCharge(2)  # Set charge state
        #    >>> feature.setOverallQuality(0.95)  # Set quality score (0-1)
        #    >>> # Access the values
        #    >>> print(f"RT: {feature.getRT()}, m/z: {feature.getMZ()}, charge: {feature.getCharge()}")

        Feature() except + nogil 
        Feature(Feature &) except + nogil 

        float getQuality(Size index) except + nogil 
            # wrap-doc:
            #  Returns the quality score in a specific dimension
            #  
            #  Args:
            #    index (int): The dimension index (0 for RT, 1 for m/z)
            #  
            #  Returns:
            #    float: Quality score for the specified dimension (typically 0-1 range)

        void setQuality(Size index, float q) except + nogil 
            # wrap-doc:
            #  Sets the quality score for a specific dimension
            #  
            #  Args:
            #    index (int): The dimension index (0 for RT, 1 for m/z)
            #    q (float): Quality score to set (typically 0-1 range)

        float getOverallQuality() except + nogil 
            # wrap-doc:
            #  Returns the overall quality score of the feature
            #  
            #  Returns:
            #    float: Overall quality score (typically 0-1, where 1 is highest quality)
            #  
            #  Note:
            #    This score represents the overall confidence in the feature detection

        void setOverallQuality(float q) except + nogil 
            # wrap-doc:
            #  Sets the overall quality score of the feature
            #  
            #  Args:
            #    q (float): Overall quality score (typically 0-1, where 1 is highest quality)

        libcpp_vector[Feature] getSubordinates() except + nogil 
            # wrap-doc:
            #  Returns subordinate features (e.g., isotopic peaks)
            #  
            #  Returns:
            #    list of Feature: List of subordinate features associated with this feature
            #  
            #  Note:
            #    Subordinate features often represent individual isotopic peaks of the same compound

        void setSubordinates(libcpp_vector[Feature]) except + nogil 
            # wrap-doc:
            #  Sets the subordinate features
            #  
            #  Args:
            #    subordinates (list of Feature): List of subordinate features to associate with this feature

        bool encloses(double rt, double mz) except + nogil  
            # wrap-doc:
            #  Checks if the feature's convex hulls enclose a given position
            #  
            #  Args:
            #    rt (float): Retention time in seconds
            #    mz (float): Mass-to-charge ratio
            #  
            #  Returns:
            #    bool: True if the position (rt, mz) is within the feature's convex hulls, False otherwise
            #  
            #  Note:
            #    This uses the feature's convex hull representation to determine spatial containment
            
        ConvexHull2D getConvexHull() except + nogil 
            # wrap-doc:
            #  Returns the overall convex hull of the feature
            #  
            #  Returns:
            #    ConvexHull2D: The overall 2D convex hull encompassing all mass traces
            #  
            #  Note:
            #    This is the union of all individual mass trace convex hulls

        libcpp_vector[ConvexHull2D] getConvexHulls() except + nogil 
            # wrap-doc:
            #  Returns the convex hulls of individual mass traces
            #  
            #  Returns:
            #    list of ConvexHull2D: List of convex hulls, one for each isotopic mass trace
            #  
            #  Note:
            #    Each isotopic peak typically has its own convex hull in RT-m/z space

        void setConvexHulls(libcpp_vector[ConvexHull2D]) except + nogil 
            # wrap-doc:
            #  Sets the convex hulls for individual mass traces
            #  
            #  Args:
            #    hulls (list of ConvexHull2D): List of convex hulls to set for this feature 

        bool operator==(Feature) except + nogil 
        bool operator!=(Feature) except + nogil 

        # from BaseFeature

        # float getQuality()  except + nogil 
        # void setQuality(float q) except + nogil 

        float getWidth() except + nogil 
            # wrap-doc:
            #  Returns the width (FWHM) of the feature in RT dimension
            #  
            #  Returns:
            #    float: Full Width at Half Maximum (FWHM) in seconds
            #  
            #  Note:
            #    Represents the elution peak width

        void setWidth(float q) except + nogil 
            # wrap-doc:
            #  Sets the width (FWHM) of the feature in RT dimension
            #  
            #  Args:
            #    q (float): Full Width at Half Maximum in seconds

        Int getCharge() except + nogil 
            # wrap-doc:
            #  Returns the charge state of the feature
            #  
            #  Returns:
            #    int: Charge state (e.g., 2 for doubly charged ions, 0 if unknown)

        void setCharge(Int q) except + nogil 
            # wrap-doc:
            #  Sets the charge state of the feature
            #  
            #  Args:
            #    q (int): Charge state (e.g., 2 for doubly charged ions)

        AnnotationState getAnnotationState() except + nogil 
            # wrap-doc:
            #  Returns the annotation state of the feature
            #  
            #  Returns:
            #    AnnotationState: Enum indicating the annotation status of this feature 

   
        PeptideIdentificationList getPeptideIdentifications() except + nogil 
            # wrap-doc:
            #  Returns the peptide identifications associated with this feature
            #  
            #  Returns:
            #    PeptideIdentificationList: List of peptide identifications from database search
            #  
            #  Note:
            #    Only relevant for peptide features. Contains results from peptide identification tools

        void setPeptideIdentifications(PeptideIdentificationList & peptides) except + nogil 
            # wrap-doc:
            #  Sets the peptide identifications for this feature
            #  
            #  Args:
            #    peptides (PeptideIdentificationList): List of peptide identifications to associate with this feature

