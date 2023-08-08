from libcpp.vector cimport vector as libcpp_vector
from libcpp.pair cimport pair as libcpp_pair
from libcpp cimport bool

from ConsensusMap cimport *
from DefaultParamHandler cimport *
from Feature cimport *
from FeatureMap cimport *
from String cimport *
from Types cimport *

from ProteinIdentification cimport *
from PeptideIdentification cimport *

from MSExperiment cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *


cdef extern from "<OpenMS/ANALYSIS/ID/IDRipper.h>" namespace "OpenMS::IDRipper":

    cdef enum OriginAnnotationFormat:
         FILE_ORIGIN, MAP_INDEX, ID_MERGE_INDEX, UNKNOWN_OAF, SIZE_OF_ORIGIN_ANNOTATION_FORMAT

    cdef cppclass IdentificationRuns:

        IdentificationRuns(
                libcpp_vector[ProteinIdentification] & prot_ids
                ) except + nogil 

        IdentificationRuns(IdentificationRuns) except + nogil    # wrap-ignore

    cdef cppclass RipFileIdentifier:

        RipFileIdentifier(
                IdentificationRuns & id_runs,
                PeptideIdentification & pep_id,
                libcpp_map[String, unsigned int] & file_origin_map,
                OriginAnnotationFormat origin_annotation_fmt,
                bool split_ident_runs) except + nogil 

        RipFileIdentifier(RipFileIdentifier) except + nogil    # wrap-ignore

        UInt getIdentRunIdx() except + nogil 

        UInt getFileOriginIdx() except + nogil 

        String getOriginFullname() except + nogil 

        String getOutputBasename() except + nogil 

    cdef cppclass RipFileContent:

        RipFileContent(
                libcpp_vector[ProteinIdentification] & prot_idents,
                libcpp_vector[PeptideIdentification] & pep_idents
                ) except + nogil 

        RipFileContent(RipFileContent) except + nogil    # wrap-ignore

        libcpp_vector[ProteinIdentification] getProteinIdentifications() except + nogil 

        libcpp_vector[PeptideIdentification] getPeptideIdentifications() except + nogil 


cdef extern from "<OpenMS/ANALYSIS/ID/IDRipper.h>" namespace "OpenMS":

    cdef cppclass IDRipper(DefaultParamHandler):
        # wrap-inherits:
        #   DefaultParamHandler
        
        IDRipper() except + nogil  # wrap-doc:Ripping protein/peptide identification according their file origin
        # private

        IDRipper(IDRipper) except + nogil    # wrap-ignore

        void rip(
                libcpp_vector[RipFileIdentifier] & rfis,
                libcpp_vector[RipFileContent] & rfcs,
                libcpp_vector[ProteinIdentification] & proteins,
                libcpp_vector[PeptideIdentification] & peptides,
                bool full_split,
                bool split_ident_runs
                ) except + nogil 

