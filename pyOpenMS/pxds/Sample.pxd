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
        Sample(Sample) nogil except + # wrap-ignore

        #returns the sample name (default: "")
        String getName() nogil except +
        #sets the sample name
        void setName(String name) nogil except +

        #returns the sample name (default: "")
        String getOrganism() nogil except +
        #sets the sample name
        void setOrganism(String organism) nogil except +

        # returns the sample number (default: "")
        String getNumber() nogil except +
        # sets the sample number (e.g. sample ID)
        void setNumber(String number) nogil except +

        # returns the comment (default: "")
        String getComment() nogil except +
        # sets the comment (may contain newline characters)
        void setComment(String comment) nogil except +

        # returns the state of aggregation (default: SAMPLENULL)
        SampleState getState() nogil except +
        # sets the state of aggregation
        void setState(SampleState state) nogil except +

        # returns the mass (in gram) (default: 0.0)
        DoubleReal getMass() nogil except +
        # sets the mass (in gram)
        void setMass(DoubleReal mass) nogil except +

        # returns the volume (in ml) (default: 0.0)
        DoubleReal getVolume() nogil except +
        # sets the volume (in ml)
        void setVolume(DoubleReal volume) nogil except +

        # returns the concentration (in g/l) (default: 0.0)
        DoubleReal getConcentration() nogil except +
        # sets the concentration (in g/l)
        void setConcentration(DoubleReal concentration) nogil except +

        # returns a reference to the vector of subsamples that were combined to create this sample
        libcpp_vector[Sample] getSubsamples() nogil except +
        # sets the vector of subsamples that were combined to create this sample
        void setSubsamples(libcpp_vector[Sample] subsamples) nogil except +

        # Since SampleTreatment is abstract, we cant wrap it
        ## void addTreatment(SampleTreatment treatment, Int before_position = -1) nogil except +
        ## SampleTreatment getTreatment(UInt position) nogil except +

        # brief removes the sample treatment at the given position
        void removeTreatment(UInt position) nogil except +
        # returns the number of sample treatments
        Int countTreatments() nogil except +

        # since this class can be derived again, we will not declare the
        # MetaInfoInterface methods here

cdef extern from "<OpenMS/METADATA/Sample.h>" namespace "OpenMS::Sample":

    cdef enum SampleState:
    
            SAMPLENULL, SOLID, LIQUID, GAS, SOLUTION, EMULSION, SUSPENSION, SIZE_OF_SAMPLESTATE

