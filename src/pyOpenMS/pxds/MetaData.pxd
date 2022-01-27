from Types cimport *
from libcpp cimport bool
from libcpp.pair cimport pair as libcpp_pair
from libcpp.vector cimport vector as libcpp_vector
from libcpp.set cimport set as libcpp_set
from cython.operator cimport dereference as deref
from MetaInfoInterface cimport *
cdef extern from "<OpenMS/METADATA/ID/MetaData.h>" namespace "OpenMS::IdentificationDataInternal":

  cdef cppclass IteratorWrapper[I,T](libcpp_set[T].iterator):
    # wrap-instances:
    #  ScoreTypeRef := IteratorWrapper[setSTit, ScoreType]
    #  ProcessingSoftwareRef := IteratorWrapper[setPSit, ProcessingSoftware]
    #  ProcessingStepRef := IteratorWrapper[setPSoftSit, ProcessingStep]
    #  SearchParamRef := IteratorWrapper[setDBSPit, DBSearchParam]
    # wrap-doc:
    #   Class for IteratorWrapper
    IteratorWrapper() nogil except +
    IteratorWrapper(IteratorWrapper[I,T]) nogil except +
    T deref() nogil except +
  
  cdef enum MoleculeType:
    PROTEIN,
    COMPOUND,
    RNA
  
  cdef enum MassType:
    MONOISOTOPIC,
    AVERAGE

cdef extern from "<boost/optional/optional.hpp>" namespace "boost":

  cdef cppclass optional[T]:
    # wrap-instances:
    #   _optional_ProcessingStepRef := optional[ProcessingStepRef]
    optional() nogil except + 
    optional(T val) nogil except +
    optional(optional[T] other) nogil except +
    bool operator==(optional[T] & other) nogil except +
    #bool operator=(optional[T] & other) nogil except +
    #bool operator=(T & other) nogil except +
    T get() nogil except +
    #bool bool() nogil except +