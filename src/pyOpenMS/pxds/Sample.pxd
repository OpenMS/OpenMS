from SourceFile cimport *
from DateTime cimport *
from DocumentIdentifier cimport *
from MetaInfoInterface cimport *
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/METADATA/Sample.h>" namespace "OpenMS":

    cdef cppclass Sample(MetaInfoInterface):
        # wrap-inherits:
        #   MetaInfoInterface

        Sample() except + nogil 
        Sample(Sample &) except + nogil 

        #returns the sample name (default: "")
        String getName() except + nogil 
        #sets the sample name
        void setName(String name) except + nogil 

        #returns the sample name (default: "")
        String getOrganism() except + nogil 
        #sets the sample name
        void setOrganism(String organism) except + nogil 

        
        String getNumber() except + nogil  # wrap-doc:Returns the sample number
        
        void setNumber(String number) except + nogil  # wrap-doc:Sets the sample number (e.g. sample ID)

        
        String getComment() except + nogil  # wrap-doc:Returns the comment (default "")
        
        void setComment(String comment) except + nogil  # wrap-doc:Sets the comment (may contain newline characters)

        
        SampleState getState() except + nogil  # wrap-doc:Returns the state of aggregation (default SAMPLENULL)
        
        void setState(SampleState state) except + nogil  # wrap-doc:Sets the state of aggregation

        
        double getMass() except + nogil  # wrap-doc:Returns the mass (in gram) (default 0.0)
        
        void setMass(double mass) except + nogil  # wrap-doc:Sets the mass (in gram)

        
        double getVolume() except + nogil  # wrap-doc:Returns the volume (in ml) (default 0.0)
        
        void setVolume(double volume) except + nogil  # wrap-doc:Sets the volume (in ml)

        
        double getConcentration() except + nogil  # wrap-doc:Returns the concentration (in g/l) (default 0.0)
        
        void setConcentration(double concentration) except + nogil  # wrap-doc:Sets the concentration (in g/l)

        
        libcpp_vector[Sample] getSubsamples() except + nogil  # wrap-doc:Returns a reference to the vector of subsamples that were combined to create this sample
        
        void setSubsamples(libcpp_vector[Sample] subsamples) except + nogil  # wrap-doc:Sets the vector of subsamples that were combined to create this sample

cdef extern from "<OpenMS/METADATA/Sample.h>" namespace "OpenMS::Sample":

    cdef enum SampleState:
        # wrap-attach:
        #    Sample
    
        SAMPLENULL, SOLID, LIQUID, GAS, SOLUTION, EMULSION, SUSPENSION, SIZE_OF_SAMPLESTATE

