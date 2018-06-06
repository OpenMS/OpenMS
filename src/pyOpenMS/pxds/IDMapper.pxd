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

cdef extern from "<OpenMS/ANALYSIS/ID/IDMapper.h>" namespace "OpenMS":

    cdef cppclass IDMapper(DefaultParamHandler):
        # wrap-inherits:
        #    DefaultParamHandler

        IDMapper() nogil except +
        IDMapper(IDMapper) nogil except +

        void annotate(MSExperiment & map_,
                      libcpp_vector[PeptideIdentification] & ids,
                      libcpp_vector[ProteinIdentification] & protein_ids,
                      bool clear_ids,
                      bool mapMS1) nogil except +

        void annotate(MSExperiment & map_,
                      FeatureMap & fmap,
                      bool clear_ids,
                      bool mapMS1) nogil except +

        void annotate(FeatureMap & map_,
                      libcpp_vector[PeptideIdentification] & ids,
                      libcpp_vector[ProteinIdentification] & protein_ids,
                      bool use_centroid_rt,
                      bool use_centroid_mz,
                      MSExperiment & spectra) nogil except +

        void annotate(ConsensusMap & map_,
                      libcpp_vector[PeptideIdentification] & ids,
                      libcpp_vector[ProteinIdentification] & protein_ids,
                      bool measure_from_subelements,
                      bool annotate_ids_with_subelements, 
                      MSExperiment & spectra) nogil except +

        IDMapper_SpectraIdentificationState mapPrecursorsToIdentifications(MSExperiment spectra,
                                                                           libcpp_vector[ PeptideIdentification ] & ids, 
                                                                           double mz_tol, double rt_tol) nogil except +

cdef extern from "<OpenMS/ANALYSIS/ID/IDMapper.h>" namespace "OpenMS::IDMapper":

    cdef enum Measure:
      MEASURE_PPM = 0,
      MEASURE_DA

cdef extern from "<OpenMS/ANALYSIS/ID/IDMapper.h>" namespace "OpenMS::IDMapper":

    cdef cppclass IDMapper_SpectraIdentificationState "OpenMS::IDMapper::SpectraIdentificationState":
        IDMapper_SpectraIdentificationState()  nogil except +
        IDMapper_SpectraIdentificationState(IDMapper_SpectraIdentificationState) nogil except + #wrap-ignore

        libcpp_vector[size_t] no_precursors
        libcpp_vector[size_t] identified
        libcpp_vector[size_t] unidentified

