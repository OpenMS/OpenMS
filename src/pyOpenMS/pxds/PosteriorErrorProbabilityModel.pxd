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
        # private
        PosteriorErrorProbabilityModel(PosteriorErrorProbabilityModel) nogil except +   # wrap-ignore

        bool fit(libcpp_vector[double] & search_engine_scores, String outlier_handling) nogil except +
            # wrap-doc:
                #   Fits the distributions to the data points(search_engine_scores). Estimated parameters for the distributions are saved in member variables
                #   computeProbability can be used afterwards
                #   Uses two Gaussians to fit. And Gauss+Gauss or Gumbel+Gauss to plot and calculate final probabilities
                #   -----
                #   :param search_engine_scores: A vector which holds the data points
                #   :returns: `true` if algorithm has run through. Else false will be returned. In that case no plot and no probabilities are calculated

        bool fit(libcpp_vector[double] & search_engine_scores, libcpp_vector[double] & probabilities, String outlier_handling) nogil except +
            # wrap-doc:
                #   Fits the distributions to the data points(search_engine_scores). Estimated parameters for the distributions are saved in member variables
                #   computeProbability can be used afterwards
                #   Uses two Gaussians to fit. And Gauss+Gauss or Gumbel+Gauss to plot and calculate final probabilities
                #   -----
                #   :param search_engine_scores: A vector which holds the data points
                #   :param probabilities a vector which holds the probability for each data point after running this function. If it has some content it will be overwritten
                #   :returns: `true` if algorithm has run through. Else false will be returned. In that case no plot and no probabilities are calculated

        void fillDensities(libcpp_vector[double] & x_scores, libcpp_vector[double] & incorrect_density, libcpp_vector[double] & correct_density) nogil except + # wrap-doc:Writes the distributions densities into the two vectors for a set of scores. Incorrect_densities represent the incorrectly assigned sequences
        void fillLogDensities(libcpp_vector[double] & x_scores, libcpp_vector[double] & incorrect_density, libcpp_vector[double] & correct_density) nogil except + # wrap-doc:Writes the log distributions densities into the two vectors for a set of scores. Incorrect_densities represent the incorrectly assigned sequences

        double computeLogLikelihood(libcpp_vector[double] & incorrect_density, libcpp_vector[double] & correct_density) nogil except + # wrap-doc:Computes the Maximum Likelihood with a log-likelihood function

        libcpp_pair[ double, double ] pos_neg_mean_weighted_posteriors(libcpp_vector[double] &x_scores,
                                                                 libcpp_vector[double] &incorrect_posteriors);
        GaussFitResult getCorrectlyAssignedFitResult() nogil except + # wrap-doc:Returns estimated parameters for correctly assigned sequences. Fit should be used before

        GaussFitResult getIncorrectlyAssignedFitResult() nogil except + # wrap-doc:Returns estimated parameters for correctly assigned sequences. Fit should be used before

        double getNegativePrior() nogil except + # wrap-doc:Returns the estimated negative prior probability

        #   @note: fit has to be used before using this function. Otherwise this function will compute nonsense.
        double computeProbability(double score) nogil except + # wrap-doc:Returns the computed posterior error probability for a given score

        TextFile initPlots(libcpp_vector[ double ] & x_scores) nogil except + # wrap-doc:Initializes the plots

        # returns the gnuplot formula of the fitted gumbel distribution. Only
        # x0 and sigma are used as local parameter alpha and scale parameter
        # beta, respectively.
        String getGumbelGnuplotFormula(GaussFitResult & params) nogil except + # wrap-doc:Returns the gnuplot formula of the fitted gumbel distribution

        String getGaussGnuplotFormula(GaussFitResult & params) nogil except + # wrap-doc:Returns the gnuplot formula of the fitted gauss distribution

        String getBothGnuplotFormula(GaussFitResult & incorrect, GaussFitResult & correct) nogil except + # wrap-doc:Returns the gnuplot formula of the fitted mixture distribution

        void plotTargetDecoyEstimation(libcpp_vector[double] & target, libcpp_vector[double] & decoy) nogil except + # wrap-doc:Plots the estimated distribution against target and decoy hits

        double getSmallestScore() nogil except + # wrap-doc:Returns the smallest score used in the last fit

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

