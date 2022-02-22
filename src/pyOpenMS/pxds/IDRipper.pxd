from libcpp.vector cimport vector as libcpp_vector
from libcpp.map cimport map as libcpp_map
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
                ) nogil except +

        IdentificationRuns(IdentificationRuns) nogil except +   # wrap-ignore

    cdef cppclass RipFileIdentifier:

        RipFileIdentifier(
                IdentificationRuns & id_runs,
                PeptideIdentification & pep_id,
                libcpp_map[String, unsigned int] & file_origin_map,
                OriginAnnotationFormat origin_annotation_fmt,
                bool split_ident_runs) nogil except +

        RipFileIdentifier(RipFileIdentifier) nogil except +   # wrap-ignore

        UInt getIdentRunIdx() nogil except +

        UInt getFileOriginIdx() nogil except +

        String getOriginFullname() nogil except +

        String getOutputBasename() nogil except +

    cdef cppclass RipFileContent:

        RipFileContent(
                libcpp_vector[ProteinIdentification] & prot_idents,
                libcpp_vector[PeptideIdentification] & pep_idents
                ) nogil except +

        RipFileContent(RipFileContent) nogil except +   # wrap-ignore

        libcpp_vector[ProteinIdentification] getProteinIdentifications() nogil except +

        libcpp_vector[PeptideIdentification] getPeptideIdentifications() nogil except +


cdef extern from "<OpenMS/ANALYSIS/ID/IDRipper.h>" namespace "OpenMS":

    cdef cppclass IDRipper(DefaultParamHandler):
        # wrap-inherits:
        #    DefaultParamHandler
        
        IDRipper() nogil except + # wrap-doc:Ripping protein/peptide identification according their file origin
        # private

        IDRipper(IDRipper) nogil except +   # wrap-ignore

        void rip(
                libcpp_vector[RipFileIdentifier] & rfis,
                libcpp_vector[RipFileContent] & rfcs,
                libcpp_vector[ProteinIdentification] & proteins,
                libcpp_vector[PeptideIdentification] & peptides,
                bool full_split,
                bool split_ident_runs
                ) nogil except +

