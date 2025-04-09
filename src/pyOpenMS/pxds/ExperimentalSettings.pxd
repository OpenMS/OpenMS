from SourceFile cimport *
from ContactPerson cimport *
from DateTime cimport *
from Sample cimport *
from DocumentIdentifier cimport *
from Instrument cimport *
from IonDetector cimport *
from HPLC cimport *
from ProteinIdentification cimport *
from MetaInfoInterface cimport *
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/METADATA/ExperimentalSettings.h>" namespace "OpenMS":

    cdef cppclass ExperimentalSettings(MetaInfoInterface, DocumentIdentifier):
        # wrap-inherits:
        #   DocumentIdentifier
        #   MetaInfoInterface
        #
        # wrap-doc:
        #  Description of the experimental settings, provides meta-information
        #  about an LC-MS/MS injection.

        ExperimentalSettings() except + nogil 
        ExperimentalSettings(ExperimentalSettings &) except + nogil 

        
        libcpp_vector[SourceFile] getSourceFiles() except + nogil  # wrap-doc:Returns a reference to the source data file
        
        void setSourceFiles(libcpp_vector[SourceFile] source_files) except + nogil  # wrap-doc:Sets the source data file

        
        DateTime getDateTime() except + nogil  # wrap-doc:Returns the date the experiment was performed
        
        void setDateTime(DateTime date_time) except + nogil  # wrap-doc:Sets the date the experiment was performed

        
        Sample getSample() except + nogil  # wrap-doc:Returns a reference to the sample description
        
        void setSample(Sample sample) except + nogil  # wrap-doc:Sets the sample description

        
        libcpp_vector[ContactPerson] getContacts() except + nogil  # wrap-doc:Returns a reference to the list of contact persons
        
        void setContacts(libcpp_vector[ContactPerson] contacts) except + nogil  # wrap-doc:Sets the list of contact persons

        
        Instrument getInstrument() except + nogil  # wrap-doc:Returns a reference to the MS instrument description
        
        void setInstrument(Instrument instrument) except + nogil  # wrap-doc:Sets the MS instrument description

        
        HPLC getHPLC() except + nogil  # wrap-doc:Returns a reference to the description of the HPLC run
        
        void setHPLC(HPLC hplc) except + nogil  # wrap-doc:Sets the description of the HPLC run

        
        String getComment() except + nogil  # wrap-doc:Returns the free-text comment
        
        void setComment(String comment) except + nogil  # wrap-doc:Sets the free-text comment

        
        libcpp_vector[ProteinIdentification] getProteinIdentifications() except + nogil  # wrap-doc:Returns a reference to the protein ProteinIdentification vector
        
        void setProteinIdentifications(libcpp_vector[ProteinIdentification] protein_identifications) except + nogil  # wrap-doc:Sets the protein ProteinIdentification vector

        
        String getFractionIdentifier() except + nogil  # wrap-doc:Returns fraction identifier
        
        void setFractionIdentifier(String fraction_identifier) except + nogil  # wrap-doc:Sets the fraction identifier

