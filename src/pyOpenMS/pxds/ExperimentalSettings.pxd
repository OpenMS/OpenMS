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

        ExperimentalSettings() nogil except + # wrap-doc:Description of the experimental settings
        ExperimentalSettings(ExperimentalSettings &) nogil except +

        # returns a mutable reference to the source data file
        libcpp_vector[SourceFile] getSourceFiles() nogil except + # wrap-doc:Returns a reference to the source data file
        # sets the source data file
        void setSourceFiles(libcpp_vector[SourceFile] source_files) nogil except + # wrap-doc:Sets the source data file

        # returns the date the experiment was performed
        DateTime getDateTime() nogil except + # wrap-doc:Returns the date the experiment was performed
        # sets the date the experiment was performed
        void setDateTime(DateTime date_time) nogil except + # wrap-doc:Sets the date the experiment was performed

        # returns a mutable reference to the sample description
        Sample getSample() nogil except + # wrap-doc:Returns a reference to the sample description
        # sets the sample description
        void setSample(Sample sample) nogil except + # wrap-doc:Sets the sample description

        # returns a mutable reference to the list of contact persons
        libcpp_vector[ContactPerson] getContacts() nogil except + # wrap-doc:Returns a reference to the list of contact persons
        # sets the list of contact persons
        void setContacts(libcpp_vector[ContactPerson] contacts) nogil except + # wrap-doc:Sets the list of contact persons

        # returns a mutable reference to the MS instrument description
        Instrument getInstrument() nogil except + # wrap-doc:Returns a reference to the MS instrument description
        # sets the MS instrument description
        void setInstrument(Instrument instrument) nogil except + # wrap-doc:Sets the MS instrument description

        # returns a mutable reference to the description of the HPLC run
        HPLC getHPLC() nogil except + # wrap-doc:Returns a reference to the description of the HPLC run
        # sets the description of the HPLC run
        void setHPLC(HPLC hplc) nogil except + # wrap-doc:Sets the description of the HPLC run

        # returns the free-text comment
        String getComment() nogil except + # wrap-doc:Returns the free-text comment
        # sets the free-text comment
        void setComment(String comment) nogil except + # wrap-doc:Sets the free-text comment

        # returns a mutable reference to the protein ProteinIdentification vector
        libcpp_vector[ProteinIdentification] getProteinIdentifications() nogil except + # wrap-doc:Returns a reference to the protein ProteinIdentification vector
        # sets the protein ProteinIdentification vector
        void setProteinIdentifications(libcpp_vector[ProteinIdentification] protein_identifications) nogil except + # wrap-doc:Sets the protein ProteinIdentification vector

        # returns fraction identifier
        String getFractionIdentifier() nogil except + # wrap-doc:Returns fraction identifier
        # sets the fraction identifier
        void setFractionIdentifier(String fraction_identifier) nogil except + # wrap-doc:Sets the fraction identifier

