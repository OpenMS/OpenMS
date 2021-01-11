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
        ContactPerson(ContactPerson) nogil except + # wrap-ignore

        # returns the first name of the person
        String getFirstName() nogil except +
        # sets the first name of the person
        void setFirstName(String name) nogil except +

        # returns the last name of the person
        String getLastName() nogil except +
        # sets the last name of the person
        void setLastName(String name) nogil except +

        # sets the full name of the person (gets split into first and last name internally)
        void setName(String name) nogil except +

        # returns the affiliation
        String getInstitution() nogil except +
        # sets the affiliation
        void setInstitution(String institution) nogil except +

        # returns the email address
        String getEmail() nogil except +
        # sets the email address
        void setEmail(String email) nogil except +

        # returns the email address
        String getURL() nogil except +
        # sets the email address
        void setURL(String email) nogil except +

        # returns the address
        String getAddress() nogil except +
        # sets the address
        void setAddress(String email) nogil except +

        # returns miscellaneous info about the contact person
        String getContactInfo() nogil except +
        # sets miscellaneous info about the contact person
        void setContactInfo(String contact_info) nogil except +


