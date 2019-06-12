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
        #    DocumentIdentifier
        #    MetaInfoInterface
        #
        # wrap-doc:
        #   Description of the experimental settings, provides meta-information
        #   about an LC-MS/MS injection.

        ExperimentalSettings() nogil except +
        ExperimentalSettings(ExperimentalSettings) nogil except + # wrap-ignore

        # returns a mutable reference to the source data file
        libcpp_vector[SourceFile] getSourceFiles() nogil except +
        # sets the source data file
        void setSourceFiles(libcpp_vector[SourceFile] source_files) nogil except +

        # returns the date the experiment was performed
        DateTime getDateTime() nogil except +
        # sets the date the experiment was performed
        void setDateTime(DateTime date_time) nogil except +

        # returns a mutable reference to the sample description
        Sample getSample() nogil except +
        # sets the sample description
        void setSample(Sample sample) nogil except +

        # returns a mutable reference to the list of contact persons
        libcpp_vector[ContactPerson] getContacts() nogil except +
        # sets the list of contact persons
        void setContacts(libcpp_vector[ContactPerson] contacts) nogil except +

        # returns a mutable reference to the MS instrument description
        Instrument getInstrument() nogil except +
        # sets the MS instrument description
        void setInstrument(Instrument instrument) nogil except +

        # returns a mutable reference to the description of the HPLC run
        HPLC getHPLC() nogil except +
        # sets the description of the HPLC run
        void setHPLC(HPLC hplc) nogil except +

        # returns the free-text comment
        String getComment() nogil except +
        # sets the free-text comment
        void setComment(String comment) nogil except +

        # returns a mutable reference to the protein ProteinIdentification vector
        libcpp_vector[ProteinIdentification] getProteinIdentifications() nogil except +
        # sets the protein ProteinIdentification vector
        void setProteinIdentifications(libcpp_vector[ProteinIdentification] protein_identifications) nogil except +

        # returns fraction identifier
        String getFractionIdentifier() nogil except +
        # sets the fraction identifier
        void setFractionIdentifier(String fraction_identifier) nogil except +

