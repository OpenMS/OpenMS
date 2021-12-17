from SourceFile cimport *
from DateTime cimport *
from DocumentIdentifier cimport *
from MetaInfoInterface cimport *
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/METADATA/Sample.h>" namespace "OpenMS":

    cdef cppclass Sample(MetaInfoInterface):
        # wrap-inherits:
        #    MetaInfoInterface

        Sample() nogil except +
        Sample(Sample &) nogil except +

        #returns the sample name (default: "")
        String getName() nogil except +
        #sets the sample name
        void setName(String name) nogil except +

        #returns the sample name (default: "")
        String getOrganism() nogil except +
        #sets the sample name
        void setOrganism(String organism) nogil except +

        
        String getNumber() nogil except + # wrap-doc:Returns the sample number
        
        void setNumber(String number) nogil except + # wrap-doc:Sets the sample number (e.g. sample ID)

        
        String getComment() nogil except + # wrap-doc:Returns the comment (default "")
        
        void setComment(String comment) nogil except + # wrap-doc:Sets the comment (may contain newline characters)

        
        SampleState getState() nogil except + # wrap-doc:Returns the state of aggregation (default SAMPLENULL)
        
        void setState(SampleState state) nogil except + # wrap-doc:Sets the state of aggregation

        
        double getMass() nogil except + # wrap-doc:Returns the mass (in gram) (default 0.0)
        
        void setMass(double mass) nogil except + # wrap-doc:Sets the mass (in gram)

        
        double getVolume() nogil except + # wrap-doc:Returns the volume (in ml) (default 0.0)
        
        void setVolume(double volume) nogil except + # wrap-doc:Sets the volume (in ml)

        
        double getConcentration() nogil except + # wrap-doc:Returns the concentration (in g/l) (default 0.0)
        
        void setConcentration(double concentration) nogil except + # wrap-doc:Sets the concentration (in g/l)

        
        libcpp_vector[Sample] getSubsamples() nogil except + # wrap-doc:Returns a reference to the vector of subsamples that were combined to create this sample
        
        void setSubsamples(libcpp_vector[Sample] subsamples) nogil except + # wrap-doc:Sets the vector of subsamples that were combined to create this sample

        # Since SampleTreatment is abstract, we cant wrap it
        ## void addTreatment(SampleTreatment treatment, Int before_position = -1) nogil except +
        ## SampleTreatment getTreatment(UInt position) nogil except +

        
        void removeTreatment(UInt position) nogil except + # wrap-doc:Brief removes the sample treatment at the given position
        
        Int countTreatments() nogil except + # wrap-doc:Returns the number of sample treatments

cdef extern from "<OpenMS/METADATA/Sample.h>" namespace "OpenMS::Sample":

    cdef enum SampleState:
        # wrap-attach:
        #     Sample
    
        SAMPLENULL, SOLID, LIQUID, GAS, SOLUTION, EMULSION, SUSPENSION, SIZE_OF_SAMPLESTATE

