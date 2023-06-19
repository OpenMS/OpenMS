from Types  cimport *
from smart_ptr cimport shared_ptr
from libcpp.string cimport string as libcpp_utf8_output_string
from libcpp.vector cimport vector as libcpp_vector
from OpenSwathDataStructures cimport *

# ctypedef shared_ptr[ISpectrumAccess] SpectrumAccessPtr
cdef extern from "<OpenMS/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>" namespace "OpenSwath":

  # This is an abstract base class in C++ with only pure virtual functions -> we can thus not create an instance of it
  cdef cppclass ISpectrumAccess:
      # wrap-ignore
      # ABSTRACT class
      # no-pxd-import

      ISpectrumAccess() nogil except + # compiler
      ISpectrumAccess(ISpectrumAccess &) nogil except + # compiler

      # virtual boost::shared_ptr<ISpectrumAccess> lightClone() const = 0;

      shared_ptr[OSSpectrum] getSpectrumById(int id_) nogil except + # wrap-doc:Returns a pointer to a spectrum at the given string id
      libcpp_vector[size_t] getSpectraByRT(double RT, double deltaRT) nogil except + # wrap-doc:Returns a vector of ids of spectra that are within RT +/- deltaRT
      size_t getNrSpectra() nogil except + # wrap-doc:Returns the number of spectra available
      # virtual SpectrumMeta getSpectrumMetaById(int id) const = 0;

      shared_ptr[OSChromatogram] getChromatogramById(int id_) nogil except + # wrap-doc:Returns a pointer to a chromatogram at the given id
      size_t getNrChromatograms() nogil except + # wrap-doc:Returns the number of chromatograms available
      libcpp_utf8_output_string getChromatogramNativeID(int id_) nogil except +
