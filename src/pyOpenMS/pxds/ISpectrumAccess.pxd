from Types  cimport *
from smart_ptr cimport shared_ptr
from libcpp.string cimport string as libcpp_string
from libcpp.vector cimport vector as libcpp_vector
from OpenSwathDataStructures cimport *

# ctypedef shared_ptr[ISpectrumAccess] SpectrumAccessPtr
cdef extern from "<OpenMS/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>" namespace "OpenSwath":

  # This is an abstract base class in C++ with only pure virtual functions -> we can thus not create an instance of it
  cdef cppclass ISpectrumAccess:
      # wrap-ignore
      # ABSTRACT class
      # no-pxd-import

      ISpectrumAccess() nogil except +
      ISpectrumAccess(ISpectrumAccess) nogil except +

      # virtual boost::shared_ptr<ISpectrumAccess> lightClone() const = 0;

      shared_ptr[OSSpectrum] getSpectrumById(int id_) nogil except +
      libcpp_vector[size_t] getSpectraByRT(double RT, double deltaRT) nogil except +
      size_t getNrSpectra() nogil except +
      # virtual SpectrumMeta getSpectrumMetaById(int id) const = 0;

      shared_ptr[OSChromatogram] getChromatogramById(int id_) nogil except +
      size_t getNrChromatograms() nogil except +
      libcpp_string getChromatogramNativeID(int id_) nogil except +

