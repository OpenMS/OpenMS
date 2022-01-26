from String cimport *
from MetaData cimport *
from libcpp.set cimport set as libcpp_set


cdef extern from "<OpenMS/METADATA/ID/InputFile.h>" namespace "OpenMS::IdentificationDataInternal":
  cdef cppclass InputFile:

    InputFile() nogil except +

    InputFile(InputFile & other) nogil except +

    String name

    String experimental_design_id

    libcpp_set[String] primary_files

    InputFile(String & name, String & experimental_design_id, libcpp_set[String] primary_files) # Note that default args don't work

    InputFile & merge(InputFile & other)

  cdef cppclass InputFiles:
    InputFiles() nogil except +


    
  cdef cppclass InputFileRef:
      InputFileRef() nogil except +
      InputFileRef(InputFileRef other) nogil except +
      bool operator!=(InputFileRef) nogil except +
      bool operator<(InputFileRef) nogil except +
      InputFile deref()
      #InputFileRef operator*() nogil
      #InputFileRef operator++() nogil