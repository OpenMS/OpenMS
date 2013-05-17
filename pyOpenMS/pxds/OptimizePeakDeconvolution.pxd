from Types cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from PeakShape cimport *
from OptimizePick cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePeakDeconvolution.h>" namespace "OpenMS":
    
    cdef cppclass OptimizePeakDeconvolution(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        OptimizePeakDeconvolution() nogil except +
        OptimizePeakDeconvolution(OptimizePeakDeconvolution) nogil except +
        # NAMESPACE # OptimizationFunctions::PenaltyFactorsIntensity  getPenalties() nogil except +
        # NAMESPACE # void setPenalties(OptimizationFunctions::PenaltyFactorsIntensity & penalties) nogil except +
        Int getCharge() nogil except +
        void setCharge(Int charge) nogil except +
        # TODO doesnt work ...
        # bool optimize(libcpp_vector[ PeakShape ] & peaks, OptimizePeakDeconvolution_Data & data) nogil except +


cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePeakDeconvolution.h>" namespace "OpenMS::OptimizationFunctions":
    
    cdef cppclass PenaltyFactorsIntensity(OptimizationFunctions_PenaltyFactors):
        # wrap-inherits:
        #  OptimizationFunctions_PenaltyFactors
        PenaltyFactorsIntensity() nogil except +
        PenaltyFactorsIntensity(PenaltyFactorsIntensity) nogil except +
        DoubleReal height


cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePeakDeconvolution.h>" namespace "OpenMS::OptimizePeakDeconvolution":

    cdef cppclass OptimizePeakDeconvolution_Data "OpenMS::OptimizePeakDeconvolution::Data":
      # TODO will not work
      # libcpp_vector[PeakShape] peaks
      libcpp_vector[double] positions
      libcpp_vector[double] signal
      # PenaltyFactorsIntensity penalties 
      Int charge 


#     property peaks:
#     
#         def __set__(self, list peaks):
#             cdef libcpp_vector[_PeakShape] * v0 = new libcpp_vector[_PeakShape]()
#             cdef PeakShape item0
#             for item0 in peaks:
#                 v0.push_back(deref(item0.inst.get()))
#             self.inst.get().peaks = deref(v0) ### Cannot convert 'vector[PeakShape]' to Python object
#             del v0
#     
#         def __get__(self):
#             _r = self.inst.get().peaks
#             py_result = []
#             cdef libcpp_vector[_PeakShape].iterator it__r = _r.begin() ### Cannot convert Python object to 'iterator'
#             cdef PeakShape item_py_result
#             while it__r != _r.end(): ### Cannot convert Python object to 'iterator'
#                item_py_result = PeakShape.__new__(PeakShape)
#                item_py_result.inst = shared_ptr[_PeakShape](new _PeakShape(deref(it__r)))
#                py_result.append(item_py_result)
#                inc(it__r)
#             return py_result

