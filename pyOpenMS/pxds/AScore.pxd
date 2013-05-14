from ConsensusMap cimport *
from DefaultParamHandler cimport *
from Feature cimport *
from FeatureMap cimport *
from String cimport *

from ProteinIdentification cimport *
from PeptideIdentification cimport *

from PeptideHit cimport *

from MSExperiment cimport *
from MSSpectrum cimport *
from RichPeak1D cimport *
from ChromatogramPeak cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/AScore.h>" namespace "OpenMS":

    cdef cppclass AScore:

        AScore() nogil except +
        AScore(AScore) nogil except +   # wrap-ignore

        PeptideHit compute(PeptideHit & hit, MSSpectrum[RichPeak1D] & real_spectrum,
                     DoubleReal fmt, int n_sites) nogil except +

        # Computes the cumulative binomial probabilities.
        DoubleReal computeCumulativeScore(UInt N, UInt n, DoubleReal p) nogil except +

        # Finds the peptides with the highest PeptideScores and outputs all informations for computing the AScore
        # This function assumes that there are more permutations than the assumed number of phosphorylations!
        # TODO nested std::vector
        #  void computeHighestPeptides(libcpp_vector[libcpp_vector[double] ] & peptide_site_scores, libcpp_vector[ProbablePhosphoSites] & sites, libcpp_vector[libcpp_vector[Size] ] & permutations) nogil except +

        #Computes the site determing_ions for the given AS and sequences in candidates
        void compute_site_determining_ions(libcpp_vector[MSSpectrum[RichPeak1D]] & th_spectra, ProbablePhosphoSites & candidates, Int charge, libcpp_vector[MSSpectrum[RichPeak1D]] & site_determining_ions) nogil except +


    cdef cppclass ProbablePhosphoSites:

        ProbablePhosphoSites() nogil except +
        ProbablePhosphoSites(ProbablePhosphoSites) nogil except +   # wrap-ignore

        Size first
        Size second
        Size seq_1
        Size seq_2
        Size peak_depth
        Size AScore

