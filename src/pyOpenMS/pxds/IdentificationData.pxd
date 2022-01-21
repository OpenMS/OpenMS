from Types cimport *
from libcpp cimport bool
from libcpp.pair cimport pair as libcpp_pair
from libcpp.vector cimport vector as libcpp_vector
from libcpp.set cimport set as libcpp_set
from MetaInfoInterface cimport *
from MetaData cimport *
from ScoreType cimport *
from ProcessingSoftware cimport *
from ProcessingStep cimport *
from InputFile cimport *



cdef extern from "<OpenMS/METADATA/ID/IdentificationData.h>" namespace "OpenMS":
    
  cdef cppclass IdentificationData(MetaInfoInterface):
      # wrap-inherits:
      #    MetaInfoInterface

      IdentificationData() nogil except + # wrap-doc:Representation of spectrum identification results and associated data

      IdentificationData(IdentificationData &) nogil except + # wrap-doc:Copy constructor

      IteratorWrapper[setSTit, ScoreType] registerScoreType(ScoreType & score) nogil except +

      IteratorWrapper[setPSit, ProcessingSoftware] registerProcessingSoftware(ProcessingSoftware & software) nogil except +

      libcpp_set[ScoreType] getScoreTypes() nogil except +

      libcpp_set[ProcessingSoftware] getProcessingSoftwares() nogil except +