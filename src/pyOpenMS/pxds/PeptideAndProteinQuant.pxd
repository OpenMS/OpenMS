from MSSpectrum cimport *
from ExperimentalDesign cimport *
from FeatureMap cimport *
from ConsensusMap cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from Param cimport *
from DefaultParamHandler cimport *
from ProgressLogger cimport *
from ProteinIdentification cimport *
# ctypedef libcpp_map<UInt64, double> SampleAbundances;

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/PeptideAndProteinQuant.h>" namespace "OpenMS":

    cdef cppclass PeptideAndProteinQuant(DefaultParamHandler):
        # wrap-inherits:
        #    DefaultParamHandler

        PeptideAndProteinQuant() nogil except +
        PeptideAndProteinQuant(PeptideAndProteinQuant) nogil except + #wrap-ignore

        void readQuantData(FeatureMap & map_in, ExperimentalDesign & ed) nogil except +
        void readQuantData(ConsensusMap & map_in, ExperimentalDesign & ed) nogil except +
        void readQuantData(libcpp_vector[ProteinIdentification] & proteins,
                           libcpp_vector[PeptideIdentification] & peptides,
                           ExperimentalDesign & ed) nogil except +

        void quantifyPeptides(libcpp_vector[PeptideIdentification] & peptides) nogil except +
        void quantifyProteins(ProteinIdentification & proteins) nogil except +

        PeptideAndProteinQuant_Statistics getStatistics() nogil except +

        # ctypedef libcpp_map<String, ProteinData] ProteinQuant
        # ctypedef libcpp_map<AASequence, PeptideData] PeptideQuant
        # PeptideQuant getPeptideResults() nogil except +
        # ProteinQuant getProteinResults() nogil except +

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/PeptideAndProteinQuant.h>" namespace "OpenMS::PeptideAndProteinQuant":

    cdef cppclass PeptideAndProteinQuant_Statistics "OpenMS::PeptideAndProteinQuant::Statistics":
      # number of samples
      Size n_samples

      # protein statistics
      Size quant_proteins
      Size too_few_peptides

      # peptide statistics
      Size quant_peptides
      Size total_peptides

      # feature statistics
      Size quant_features
      Size total_features
      Size blank_features
      Size ambig_features

      # constructor
      PeptideAndProteinQuant_Statistics()  nogil except +
      PeptideAndProteinQuant_Statistics(PeptideAndProteinQuant_Statistics) # wrap-ignore

    cdef cppclass PeptideAndProteinQuant_PeptideData "OpenMS::PeptideAndProteinQuant::PeptideData":

      # libcpp_map[Int, SampleAbundances] abundances
      # libcpp_map[unsigned long, double] total_abundances

      # protein accessions for this peptide
      libcpp_set[String] accessions

      # number of identifications
      Size id_count

      # constructor
      PeptideAndProteinQuant_PeptideData()  nogil except +
      PeptideAndProteinQuant_PeptideData(PeptideAndProteinQuant_PeptideData) # wrap-ignore

    # Quantitative and associated data for a protein
    cdef cppclass PeptideAndProteinQuant_ProteinData "OpenMS::PeptideAndProteinQuant::ProteinData":

      # libcpp_map[String, SampleAbundances] abundances
 
      # compile issue 
      # libcpp_map[unsigned long, double] total_abundances

      # total number of identifications (of peptides mapping to this protein)
      Size id_count

      # constructor
      PeptideAndProteinQuant_ProteinData()  nogil except +
      PeptideAndProteinQuant_ProteinData(PeptideAndProteinQuant_ProteinData) # wrap-ignore

