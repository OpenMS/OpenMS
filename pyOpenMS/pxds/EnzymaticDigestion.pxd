from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from String cimport *
from AASequence cimport *

cdef extern from "<OpenMS/CHEMISTRY/EnzymaticDigestion.h>" namespace "OpenMS":


    cdef cppclass EnzymaticDigestion:

      EnzymaticDigestion() nogil except +

      SignedSize getMissedCleavages() nogil except +
      void setMissedCleavages(SignedSize missed_cleavages) nogil except +

      Enzyme getEnzyme()  nogil except +
      void setEnzyme(Enzyme enzyme) nogil except +
      Enzyme getEnzymeByName(String & name) nogil except +

      void digest(AASequence & protein, libcpp_vector[AASequence] & output) nogil except +
      Size peptideCount(AASequence & protein) nogil except +

      bool isLogModelEnabled() nogil except +
      void setLogModelEnabled(bool enabled) nogil except +

      DoubleReal getLogThreshold() nogil except +
      void setLogThreshold(DoubleReal threshold) nogil except +

cdef extern from "<OpenMS/CHEMISTRY/EnzymaticDigestion.h>" namespace "OpenMS::EnzymaticDigestion":

    # protected
    # cdef cppclass BindingSite:
    #   BindingSite()
    #   Size position
    #   String AAname

    # cdef cppclass CleavageModel:
    #   CleavageModel()
    #   DoubleReal p_cleave
    #   DoubleReal p_miss

    cdef enum Enzyme:
        # wrap-attach:
        #    EnzymaticDigestion
        TRYPSIN,
        SIZE_OF_TRYPSIN

