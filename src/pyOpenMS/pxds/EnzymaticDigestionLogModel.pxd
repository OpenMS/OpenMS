from Types cimport *
from String cimport *
from AASequence cimport *

cdef extern from "<OpenMS/CHEMISTRY/EnzymaticDigestionLogModel.h>" namespace "OpenMS":


    cdef cppclass EnzymaticDigestionLogModel:
    # wrap-doc:
    # Class for the Log L model of enzymatic digestion of proteins
    
      EnzymaticDigestionLogModel() except + nogil 
      EnzymaticDigestionLogModel(EnzymaticDigestionLogModel &) except + nogil  

      # not wrapped due to name clash with Enzyme.h
      # Enzyme getEnzyme()  except + nogil 
      # void setEnzyme(Enzyme enzyme) except + nogil 
      # Enzyme getEnzymeByName(String name) except + nogil 

      String getEnzymeName() except + nogil  # wrap-doc:Returns the enzyme for the digestion
      void setEnzyme(String name) except + nogil  # wrap-doc:Sets the enzyme for the digestion

      double getLogThreshold() except + nogil  # wrap-doc:Returns the threshold which needs to be exceeded to call a cleavage (only for the trained cleavage model on real data)
      void setLogThreshold(double threshold) except + nogil  # wrap-doc:Sets the threshold which needs to be exceeded to call a cleavage (only for the trained cleavage model on real data). Default is 0.25

      void digest(AASequence & protein, libcpp_vector[AASequence] & output) except + nogil  # wrap-doc:Performs the enzymatic digestion of a protein
      Size peptideCount(AASequence & protein) except + nogil  
          # wrap-doc:
          #  Returns the number of peptides a digestion of `protein` would yield under the current enzyme and missed cleavage settings
          #  
          #  
          #  :param protein: Name of protein

# cdef extern from "<OpenMS/CHEMISTRY/EnzymaticDigestionLogModel.h>" namespace "OpenMS::EnzymaticDigestionLogModel":
# 
#    # protected
#    # cdef cppclass BindingSite:
#    #   BindingSite()
#    #   Size position
#    #   String AAname
# 
#    # cdef cppclass CleavageModel:
#    #   CleavageModel()
#    #   double p_cleave
#    #   double p_miss
# 
