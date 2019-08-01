from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool
from libcpp.pair cimport pair as libcpp_pair

from Types cimport *
from DefaultParamHandler cimport *
from String cimport *

from ConsensusMap cimport *
from Feature cimport *
from FeatureMap cimport *
from MSExperiment cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *

from ProteinIdentification cimport *
from PeptideIdentification cimport *
from GaussFitter cimport *
from TextFile cimport *

cdef extern from "<OpenMS/MATH/STATISTICS/PosteriorErrorProbabilityModel.h>" namespace "OpenMS::Math":

    cdef cppclass PosteriorErrorProbabilityModel(DefaultParamHandler):
        # wrap-inherits:
        #    DefaultParamHandler

        PosteriorErrorProbabilityModel() nogil except +
        PosteriorErrorProbabilityModel(PosteriorErrorProbabilityModel) nogil except +   # wrap-ignore

        bool fit(libcpp_vector[double] & search_engine_scores) nogil except +
        bool fit(libcpp_vector[double] & search_engine_scores, libcpp_vector[double] & probabilities) nogil except +

        #Writes the distributions densities into the two vectors for a set of scores. Incorrect_densities represent the incorrectly assigned sequences.
        void fillDensities(libcpp_vector[double] & x_scores, libcpp_vector[double] & incorrect_density, libcpp_vector[double] & correct_density) nogil except +
        void fillLogDensities(libcpp_vector[double] & x_scores, libcpp_vector[double] & incorrect_density, libcpp_vector[double] & correct_density) nogil except +

        #computes the Maximum Likelihood with a log-likelihood function.
        double computeLogLikelihood(libcpp_vector[double] & incorrect_density, libcpp_vector[double] & correct_density) nogil except +

        libcpp_pair[ double, double ] pos_neg_mean_weighted_posteriors(libcpp_vector[double] &x_scores,
                                                                 libcpp_vector[double] &incorrect_posteriors);
        #returns estimated parameters for correctly assigned sequences. Fit should be used before.
        GaussFitResult getCorrectlyAssignedFitResult() nogil except +

        #returns estimated parameters for correctly assigned sequences. Fit should be used before.
        GaussFitResult getIncorrectlyAssignedFitResult() nogil except +

        #returns the estimated negative prior probability.
        double getNegativePrior() nogil except +

        #   Returns the computed posterior error probability for a given score.
        #   @note: fit has to be used before using this function. Otherwise this function will compute nonsense.
        double computeProbability(double score) nogil except +

        #initializes the plots
        TextFile initPlots(libcpp_vector[ double ] & x_scores) nogil except +

        # returns the gnuplot formula of the fitted gumbel distribution. Only
        # x0 and sigma are used as local parameter alpha and scale parameter
        # beta, respectively.
        String getGumbelGnuplotFormula(GaussFitResult & params) nogil except +

        # returns the gnuplot formula of the fitted gauss distribution.
        String getGaussGnuplotFormula(GaussFitResult & params) nogil except +

        # returns the gnuplot formula of the fitted mixture distribution.
        String getBothGnuplotFormula(GaussFitResult & incorrect, GaussFitResult & correct) nogil except +

        #plots the estimated distribution against target and decoy hits
        void plotTargetDecoyEstimation(libcpp_vector[double] & target, libcpp_vector[double] & decoy) nogil except +

        # returns the smallest score used in the last fit
        double getSmallestScore() nogil except +

        void tryGnuplot(const String & gp_file) nogil except +

        # Cannot handle bool by reference ... 
        # void updateScores(PosteriorErrorProbabilityModel PEP_model,
        #                   const String & search_engine,
        #                   const Int charge,
        #                   const bool prob_correct,
        #                   const bool split_charge,
        #                   libcpp_vector[ ProteinIdentification ] & protein_ids,
        #                   libcpp_vector[ PeptideIdentification ] & peptide_ids,
        #                   bool unable_to_fit_data,
        #                   bool data_might_not_be_well_fit) nogil except +

