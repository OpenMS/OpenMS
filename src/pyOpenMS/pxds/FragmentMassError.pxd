from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from String cimport *
from QCBase cimport *
from MSExperiment cimport *
from FeatureMap cimport *
from ProteinIdentification cimport *
from PeptideIdentificationList cimport *

cdef extern from "<OpenMS/QC/FragmentMassError.h>" namespace "OpenMS":

    cdef cppclass FragmentMassError(QCBase):
        # wrap-inherits:
        #    QCBase
        #
        # wrap-doc:
        #  QC metric for fragment mass error (FME)
        #
        #  Calculates the fragment mass error in ppm or Da for each match of a
        #  theoretical spectrum to an experimental spectrum.

        FragmentMassError() except + nogil
        FragmentMassError(FragmentMassError &) except + nogil

        void compute(FeatureMap & fmap,
                     MSExperiment & exp,
                     SpectraMap & map_to_spectrum,
                     ToleranceUnit tolerance_unit,
                     double tolerance) except + nogil
            # wrap-doc:
            #  Computes FragmentMassError (FME) in ppm and Dalton using FeatureMap
            #
            #  Stores average FME over all spectra and its variance in ppm.
            #  Each FME is stored at the first PeptideHit of the corresponding
            #  PeptideIdentification as metavalue.
            #
            #  :param fmap: Input FeatureMap for annotation and data for theoretical spectra
            #  :param exp: Input MSExperiment for MS2 spectra
            #  :param map_to_spectrum: Map to find index of spectrum given by meta value at PepID
            #  :param tolerance_unit: Tolerance in ppm or Dalton (AUTO uses FeatureMap metadata)
            #  :param tolerance: Search window for matching peaks

        void compute(PeptideIdentificationList & pep_ids,
                     SearchParameters & search_params,
                     MSExperiment & exp,
                     SpectraMap & map_to_spectrum,
                     ToleranceUnit tolerance_unit,
                     double tolerance) except + nogil
            # wrap-doc:
            #  Computes FragmentMassError (FME) in ppm and Dalton using PeptideIdentifications
            #
            #  Stores average FME over all spectra and its variance in ppm.
            #  Each FME is stored at the first PeptideHit of the corresponding
            #  PeptideIdentification as metavalue.
            #
            #  :param pep_ids: Input PeptideIdentifications for annotation
            #  :param search_params: Search parameters for finding fragment mass tolerance
            #  :param exp: Input MSExperiment for MS2 spectra
            #  :param map_to_spectrum: Map to find index of spectrum given by meta value at PepID
            #  :param tolerance_unit: Tolerance in ppm or Dalton (AUTO uses search_params)
            #  :param tolerance: Search window for matching peaks

        const String & getName() except + nogil  # wrap-doc:Returns the name of the metric

        libcpp_vector[Statistics] getResults() except + nogil  # wrap-doc:Returns results

cdef extern from "<OpenMS/QC/FragmentMassError.h>" namespace "OpenMS::FragmentMassError":

    cdef cppclass Statistics "OpenMS::FragmentMassError::Statistics":
        # wrap-doc:
        #  Structure for storing fragment mass error statistics

        Statistics() except + nogil
        Statistics(Statistics &) except + nogil

        double average_ppm  # wrap-doc:Average fragment mass error in ppm
        double variance_ppm  # wrap-doc:Variance of fragment mass error in ppm
