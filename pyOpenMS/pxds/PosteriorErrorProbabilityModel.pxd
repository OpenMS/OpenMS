from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool

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

# TODO vector<DoubleReal] doesnt get resolved 

cdef extern from "<OpenMS/MATH/STATISTICS/PosteriorErrorProbabilityModel.h>" namespace "OpenMS::Math":

    cdef cppclass PosteriorErrorProbabilityModel(DefaultParamHandler):
        # wrap-inherits:
        #    DefaultParamHandler

        PosteriorErrorProbabilityModel() nogil except +
        PosteriorErrorProbabilityModel(PosteriorErrorProbabilityModel) nogil except +   # wrap-ignore

        bool fit(libcpp_vector[double] & search_engine_scores) nogil except +
        bool fit(libcpp_vector[double] & search_engine_scores, libcpp_vector[double] & probabilities) nogil except +

        #Writes the distributions densities into the two vectors for a set of scores. Incorrect_densities represent the incorreclty assigned seqeuences.
        void fillDensities(libcpp_vector[double] & x_scores, libcpp_vector[double] & incorrect_density, libcpp_vector[double] & correct_density) nogil except +
        #computes the Maximum Likelihood with a log-likelihood funciotn.
        DoubleReal computeMaxLikelihood(libcpp_vector[double] & incorrect_density, libcpp_vector[double] & correct_density) nogil except +

        #sums (1 - posterior porbabilities)
        DoubleReal one_minus_sum_post(libcpp_vector[double] & incorrect_density, libcpp_vector[double] & correct_density) nogil except +
        #sums  posterior porbabilities
        DoubleReal sum_post(libcpp_vector[double] & incorrect_density, libcpp_vector[double] & correct_density) nogil except +
        #helper function for the EM algorithm (for fitting)
        DoubleReal sum_pos_x0(libcpp_vector[double] & x_scores, libcpp_vector[double] & incorrect_density, libcpp_vector[double] & correct_density) nogil except +
        #helper function for the EM algorithm (for fitting)
        DoubleReal sum_neg_x0(libcpp_vector[double] & x_scores, libcpp_vector[double] & incorrect_density, libcpp_vector[double] & correct_density) nogil except +
        #helper function for the EM algorithm (for fitting)
        DoubleReal sum_pos_sigma(libcpp_vector[double] & x_scores, libcpp_vector[double] & incorrect_density, libcpp_vector[double] & correct_density, DoubleReal positive_mean) nogil except +
        #helper function for the EM algorithm (for fitting)
        DoubleReal sum_neg_sigma(libcpp_vector[double] & x_scores, libcpp_vector[double] & incorrect_density, libcpp_vector[double] & correct_density, DoubleReal positive_mean) nogil except +

        #returns estimated parameters for correctly assigned sequences. Fit should be used before.
        GaussFitResult getCorrectlyAssignedFitResult()

        #returns estimated parameters for correctly assigned sequences. Fit should be used before.
        GaussFitResult getIncorrectlyAssignedFitResult()

        #returns the estimated negative prior probability.
        DoubleReal getNegativePrior()

        #computes the gaussian density at position x with parameters params.
        DoubleReal getGauss(DoubleReal x, GaussFitResult & params)

        #computes the gumbel density at position x with parameters params.
        DoubleReal getGumbel(DoubleReal x, GaussFitResult & params)

        #   Returns the computed posterior error probability for a given score.
        #   @note: fit has to be used before using this function. Otherwise this function will compute nonsense.
        DoubleReal computeProbability(DoubleReal score) nogil except +

        #initializes the plots
        # TODO raw ptr
        ## TextFile * InitPlots(libcpp_vector[double] & x_scores) nogil except +

        # returns the gnuplot formula of the fitted gumbel distribution. Only x0 and sigma are used as local parameter alpha and scale parameter beta, respectively.
        String getGumbelGnuplotFormula(GaussFitResult & params) nogil except +

        # returns the gnuplot formula of the fitted gauss distribution.
        String getGaussGnuplotFormula(GaussFitResult & params) nogil except +

        # returns the gnuplot formula of the fitted mixture distribution.
        String getBothGnuplotFormula(GaussFitResult & incorrect, GaussFitResult & correct) nogil except +

        #plots the estimated distribution against target and decoy hits
        void plotTargetDecoyEstimation(libcpp_vector[double] & target, libcpp_vector[double] & decoy) nogil except +

        # returns the smallest score used in the last fit
        DoubleReal getSmallestScore()

