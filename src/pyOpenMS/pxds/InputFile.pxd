from String cimport *
from MetaData cimport *
from libcpp.set cimport set as libcpp_set


cdef extern from "<OpenMS/METADATA/ID/InputFile.h>" namespace "OpenMS::IdentificationDataInternal":
  cdef cppclass InputFile:

    InputFile() nogil except +

    InputFile(InputFile other) nogil except +

    String name

    String experimental_design_id

    libcpp_set[String] primary_files

    #InputFile(String & name, String & experimental_design_id = "", libcpp_set[String] primary_files = libcpp_set[String]())

    InputFile & merge(InputFile & other)

  ctypedef libcpp_set[ InputFile ] InputFiles #FIXME this is a placeholder
  ctypedef libcpp_set[ InputFile ].iterator setIFit # FIXME TOO 