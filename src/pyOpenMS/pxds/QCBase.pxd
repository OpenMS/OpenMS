from libcpp cimport bool
from libcpp.map cimport map as libcpp_map
from Types cimport *
from String cimport *
from MSExperiment cimport *
from ConsensusMap cimport *

cdef extern from "<OpenMS/QC/QCBase.h>" namespace "OpenMS":
    
    cdef cppclass QCBase:
        # wrap-doc:
        #  Base class for all QC classes
        #  
        #  This class serves as an abstract base class for all QC classes.
        #  It contains the important feature of encoding the input requirements
        #  for a certain QC.

        # QCBase() is abstract, no constructor needed
        
        # Enum for input file requirements
        # Note: This is wrapped as nested class
        
        const String & getName() except + nogil  # wrap-doc:Returns the name of the metric
        
        # Status requirements() except + nogil  # wrap-doc:Returns the input data requirements of the compute(...) function
        
        # bool isRunnable(Status & s) except + nogil  # wrap-doc:Tests if a metric has the required input files
        
        @staticmethod
        bool isLabeledExperiment(ConsensusMap & cm) except + nogil  # wrap-doc:Check if the IsobaricAnalyzer TOPP tool was used to create this ConsensusMap

cdef extern from "<OpenMS/QC/QCBase.h>" namespace "OpenMS::QCBase":
    
    cdef enum Requires "OpenMS::QCBase::Requires":
        #wrap-attach:
        #   QCBase
        NOTHING,
        RAWMZML,
        POSTFDRFEAT,
        PREFDRFEAT,
        CONTAMINANTS,
        TRAFOALIGN,
        ID,
        SIZE_OF_REQUIRES
    
    cdef enum ToleranceUnit "OpenMS::QCBase::ToleranceUnit":
        #wrap-attach:
        #   QCBase
        AUTO,
        PPM,
        DA,
        SIZE_OF_TOLERANCEUNIT

    cdef cppclass SpectraMap "OpenMS::QCBase::SpectraMap":
        # wrap-doc:
        #  Map to find a spectrum via its NativeID
        
        SpectraMap() except + nogil 
        SpectraMap(SpectraMap &) except + nogil 
        SpectraMap(MSExperiment & exp) except + nogil  # wrap-doc:Constructor which allows immediate indexing of an MSExperiment
        
        void calculateMap(MSExperiment & exp) except + nogil  # wrap-doc:Calculate a new map, delete the old one
        
        UInt64 at(String & identifier) except + nogil  # wrap-doc:Get index from identifier
        
        void clear() except + nogil  # wrap-doc:Clear the map
        
        bool empty() except + nogil  # wrap-doc:Check if empty
        
        Size size() except + nogil  # wrap-doc:Get size of map
