from String cimport *
from MetaData cimport *
from libcpp.set cimport set as libcpp_set
from libcpp.map cimport map as libcpp_map
from ProcessingStep cimport * 
from DigestionEnzyme cimport *
from EnzymaticDigestion cimport *
from MetaData cimport *


cdef extern from "<OpenMS/METADATA/ID/DBSearchParam.h>" namespace "OpenMS::IdentificationDataInternal":
  cdef cppclass DBSearchParam:

    MoleculeType molecule_type
    MassType mass_type
    String database
    String database_version
    String taxonomy

    libcpp_set[int] charges

    libcpp_set[String] fixed_mods
    libcpp_set[String] variable_mods

    double precursor_mass_tolerance
    double fragment_mass_tolerance
    bool precursor_tolerance_ppm
    bool fragment_tolerance_ppm

    const DigestionEnzyme* digestion_enzyme

    Specificity enzyme_term_specificity
    Size missed_cleavages
    Size min_length
    Size max_length

    DBSearchParam() nogil except +

    DBSearchParam(const DBSearchParam & other) nogil except +

    bool operator<(const DBSearchParam & other) nogil except +

    bool operator==(const DBSearchParam & other) nogil except +

  ctypedef libcpp_set[ DBSearchParam ] DBSearchParams
  ctypedef libcpp_set[ DBSearchParam ].iterator setDBSPit
  ctypedef libcpp_map[IteratorWrapper[setPSoftSit, ProcessingStep], IteratorWrapper[setDBSPit, DBSearchParam]] DBSearchSteps