from libcpp.vector cimport vector as libcpp_vector
from libcpp.map cimport map as libcpp_map
from libcpp.pair cimport pair as libcpp_pair
from libcpp cimport bool

from ConsensusMap cimport *
from DefaultParamHandler cimport *
from Feature cimport *
from FeatureMap cimport *
from String cimport *

from ProteinIdentification cimport *
from PeptideIdentification cimport *

from MSExperiment cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/IDRipper.h>" namespace "OpenMS":

    cdef cppclass IDRipper(DefaultParamHandler):
        # wrap-inherits:
        #    DefaultParamHandler

        IDRipper() nogil except +
        IDRipper(IDRipper) nogil except +   # wrap-ignore

        # see additional pyx file in ./addons
        void rip( 
          libcpp_map[String, libcpp_pair[ libcpp_vector[ProteinIdentification],
          libcpp_vector[PeptideIdentification]]] & ripped,
          libcpp_vector[ProteinIdentification] & proteins,
          libcpp_vector[PeptideIdentification] & peptides) # wrap-ignore

