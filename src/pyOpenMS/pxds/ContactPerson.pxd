from SourceFile cimport *
from DateTime cimport *
from DocumentIdentifier cimport *
from libcpp.vector cimport vector as libcpp_vector
from MetaInfoInterface cimport *

cdef extern from "<OpenMS/METADATA/ContactPerson.h>" namespace "OpenMS":

    cdef cppclass ContactPerson(MetaInfoInterface):
        # wrap-inherits:
        #    MetaInfoInterface

        ContactPerson() nogil except +
        ContactPerson(ContactPerson &) nogil except +

        # returns the first name of the person
        String getFirstName() nogil except + # wrap-doc:Returns the first name of the person
        # sets the first name of the person
        void setFirstName(String name) nogil except + # wrap-doc:Sets the first name of the person

        # returns the last name of the person
        String getLastName() nogil except + # wrap-doc:Returns the last name of the person
        # sets the last name of the person
        void setLastName(String name) nogil except + # wrap-doc:Sets the last name of the person

        # sets the full name of the person (gets split into first and last name internally)
        void setName(String name) nogil except + # wrap-doc:Sets the full name of the person (gets split into first and last name internally)

        # returns the affiliation
        String getInstitution() nogil except + # wrap-doc:Returns the affiliation
        # sets the affiliation
        void setInstitution(String institution) nogil except + # wrap-doc:Sets the affiliation

        # returns the email address
        String getEmail() nogil except + # wrap-doc:Returns the email address
        # sets the email address
        void setEmail(String email) nogil except + # wrap-doc:Sets the email address

        # returns the email address
        String getURL() nogil except + # wrap-doc:Returns the URL associated with the contact person (e.g., the institute webpage "https://www.higglesworth.edu/")
        # sets the email address
        void setURL(String email) nogil except + # wrap-doc:Sets the URL associated with the contact person (e.g., the institute webpage "https://www.higglesworth.edu/")

        # returns the address
        String getAddress() nogil except + # wrap-doc:Returns the address
        # sets the address
        void setAddress(String email) nogil except + # wrap-doc:Sets the address

        # returns miscellaneous info about the contact person
        String getContactInfo() nogil except + # wrap-doc:Returns miscellaneous info about the contact person
        # sets miscellaneous info about the contact person
        void setContactInfo(String contact_info) nogil except + # wrap-doc:Sets miscellaneous info about the contact person


