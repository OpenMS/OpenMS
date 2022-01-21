from libcpp.set cimport set as libcpp_set
from CVTerm cimport *
from String cimport *
from libcpp cimport bool
from MetaInfoInterface cimport *

cdef extern from "<OpenMS/METADATA/ID/ScoreType.h>" namespace "OpenMS::IdentificationDataInternal":

  ctypedef libcpp_set[ ScoreType ].iterator setSTit

  cdef cppclass ScoreType(MetaInfoInterface):
    CVTerm cv_term
    
    bool higher_better
    
    ScoreType() nogil except +
    
    ScoreType(CVTerm cv_term, bool higher_better) nogil except +

    ScoreType(String name, bool higher_better) nogil except +

    ScoreType(ScoreType other) nogil except +

    bool operator<(ScoreType other) nogil except +

    bool operator==(ScoreType other) nogil except +

    bool isBetterScore(double first, double second) nogil except +

  ctypedef libcpp_set[ ScoreType ] ScoreTypes

