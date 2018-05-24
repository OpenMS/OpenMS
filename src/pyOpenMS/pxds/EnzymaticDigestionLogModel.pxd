from Types cimport *
from String cimport *
from AASequence cimport *

cdef extern from "<OpenMS/CHEMISTRY/EnzymaticDigestionLogModel.h>" namespace "OpenMS":


    cdef cppclass EnzymaticDigestionLogModel:

      EnzymaticDigestionLogModel() nogil except +

      # not wrapped due to name clash with Enzyme.h
      # Enzyme getEnzyme()  nogil except +
      # void setEnzyme(Enzyme enzyme) nogil except +
      # Enzyme getEnzymeByName(String name) nogil except +

      String getEnzymeName() nogil except +
      void setEnzyme(String name) nogil except +

      double getLogThreshold() nogil except +
      void setLogThreshold(double threshold) nogil except +

      void digest(AASequence & protein, libcpp_vector[AASequence] & output) nogil except +
      Size peptideCount(AASequence & protein) nogil except +

# cdef extern from "<OpenMS/CHEMISTRY/EnzymaticDigestionLogModel.h>" namespace "OpenMS::EnzymaticDigestionLogModel":
# 
#     # protected
#     # cdef cppclass BindingSite:
#     #   BindingSite()
#     #   Size position
#     #   String AAname
# 
#     # cdef cppclass CleavageModel:
#     #   CleavageModel()
#     #   double p_cleave
#     #   double p_miss
# 
