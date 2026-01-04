from TransitionTSVFile cimport *
from TargetedExperiment cimport *
from DataValue cimport *
from DefaultParamHandler cimport *
from MSExperiment cimport *
from MSSpectrum cimport *
from FeatureMap cimport *
from String cimport *
from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from libcpp.pair cimport pair as libcpp_pair
from libcpp.map cimport map as libcpp_map

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/TargetedSpectraExtractor.h>" namespace "OpenMS":

    cdef cppclass TargetedSpectraExtractor(DefaultParamHandler):
        # wrap-inherits:
        #  DefaultParamHandler
        # wrap-doc:
        #   Filter, annotate, pick, and score spectra based on a target list
        #
        #   This class processes spectra from DDA experiments against a target transition list.
        #   It provides a complete pipeline from raw spectra to scored, selected spectra.
        #
        #   Workflow:
        #   1. annotateSpectra() - Match spectra to targets by precursor MZ and RT
        #   2. pickSpectrum() - Smooth and pick peaks
        #   3. scoreSpectra() - Score by TIC, FWHM, SNR with configurable weights
        #   4. selectSpectra() - Choose the best spectrum per transition
        #
        #   Or use extractSpectra() for the full pipeline in one call.
        #
        #   Example usage:
        #
        #   .. code-block:: python
        #
        #       from pyopenms import *
        #
        #       # Load experiment and targets
        #       exp = MSExperiment()
        #       MzMLFile().load("experiment.mzML", exp)
        #       targets = TargetedExperiment()
        #       TraMLFile().load("targets.traML", targets)
        #
        #       # Extract matching spectra
        #       extractor = TargetedSpectraExtractor()
        #       extracted = []
        #       features = FeatureMap()
        #       extractor.extractSpectra(exp, targets, extracted, features, True)
        #
        #       print(f"Found {len(extracted)} matching spectra")

        TargetedSpectraExtractor() except + nogil
        TargetedSpectraExtractor(TargetedSpectraExtractor &) except + nogil

        void getDefaultParameters(Param &) except + nogil
            # wrap-doc:
            #   Get default parameters for the extractor

        void annotateSpectra(
            libcpp_vector[MSSpectrum] & spectra,
            TargetedExperiment & targeted_exp,
            libcpp_vector[MSSpectrum] & annotated_spectra,
            FeatureMap & features,
            bool compute_features
        ) except + nogil
            # wrap-doc:
            #   Annotate spectra that match target transitions
            #
            #   Filters spectra that fall within the precursor RT window and MZ tolerance.
            #
            #   :param spectra: Input spectra to filter
            #   :param targeted_exp: Target transition list
            #   :param annotated_spectra: Output annotated spectra
            #   :param features: Output FeatureMap with transition info
            #   :param compute_features: If True, populate features output

        void annotateSpectra(
            libcpp_vector[MSSpectrum] & spectra,
            TargetedExperiment & targeted_exp,
            libcpp_vector[MSSpectrum] & annotated_spectra
        ) except + nogil
            # wrap-doc:
            #   Annotate spectra that match target transitions (without features)
            #
            #   :param spectra: Input spectra to filter
            #   :param targeted_exp: Target transition list
            #   :param annotated_spectra: Output annotated spectra

        void annotateSpectra(
            libcpp_vector[MSSpectrum] & spectra,
            FeatureMap & ms1_features,
            FeatureMap & ms2_features,
            libcpp_vector[MSSpectrum] & annotated_spectra
        ) except + nogil
            # wrap-doc:
            #   Annotate spectra using MS1 and MS2 feature maps
            #
            #   :param spectra: Input spectra
            #   :param ms1_features: MS1 feature map
            #   :param ms2_features: MS2 feature map
            #   :param annotated_spectra: Output annotated spectra

        void searchSpectrum(
            FeatureMap & ms1_features,
            FeatureMap & ms2_features,
            bool add_unidentified_features
        ) except + nogil
            # wrap-doc:
            #   Search for matching features between MS1 and MS2
            #
            #   :param ms1_features: MS1 feature map
            #   :param ms2_features: MS2 feature map (modified in place)
            #   :param add_unidentified_features: Add features without matches

        void pickSpectrum(MSSpectrum & spectrum, MSSpectrum & picked_spectrum) except + nogil
            # wrap-doc:
            #   Smooth and pick peaks from a spectrum
            #
            #   Applies Gaussian or Savitzky-Golay smoothing followed by peak picking.
            #
            #   :param spectrum: Input spectrum
            #   :param picked_spectrum: Output picked spectrum

        void scoreSpectra(
            libcpp_vector[MSSpectrum] & annotated_spectra,
            libcpp_vector[MSSpectrum] & picked_spectra,
            FeatureMap & features,
            libcpp_vector[MSSpectrum] & scored_spectra
        ) except + nogil
            # wrap-doc:
            #   Score spectra based on TIC, FWHM, and SNR
            #
            #   Scores are computed using configurable weights for each metric.
            #
            #   :param annotated_spectra: Annotated input spectra
            #   :param picked_spectra: Peak-picked spectra
            #   :param features: Feature map with transition info
            #   :param scored_spectra: Output scored spectra

        void scoreSpectra(
            libcpp_vector[MSSpectrum] & annotated_spectra,
            libcpp_vector[MSSpectrum] & picked_spectra,
            libcpp_vector[MSSpectrum] & scored_spectra
        ) except + nogil
            # wrap-doc:
            #   Score spectra (without feature map)
            #
            #   :param annotated_spectra: Annotated input spectra
            #   :param picked_spectra: Peak-picked spectra
            #   :param scored_spectra: Output scored spectra

        void selectSpectra(
            libcpp_vector[MSSpectrum] & scored_spectra,
            FeatureMap & features,
            libcpp_vector[MSSpectrum] & selected_spectra,
            FeatureMap & selected_features
        ) except + nogil
            # wrap-doc:
            #   Select the best spectrum for each transition
            #
            #   :param scored_spectra: Scored input spectra
            #   :param features: Feature map with transition info
            #   :param selected_spectra: Output selected spectra
            #   :param selected_features: Output selected features

        void selectSpectra(
            libcpp_vector[MSSpectrum] & scored_spectra,
            libcpp_vector[MSSpectrum] & selected_spectra
        ) except + nogil
            # wrap-doc:
            #   Select the best spectrum for each transition (without features)
            #
            #   :param scored_spectra: Scored input spectra
            #   :param selected_spectra: Output selected spectra

        void extractSpectra(
            MSExperiment & experiment,
            TargetedExperiment & targeted_exp,
            libcpp_vector[MSSpectrum] & extracted_spectra,
            FeatureMap & features,
            bool compute_features
        ) except + nogil
            # wrap-doc:
            #   Complete pipeline: annotate, pick, score, and select spectra
            #
            #   Runs the full extraction workflow in a single call.
            #
            #   :param experiment: Input MS experiment
            #   :param targeted_exp: Target transition list
            #   :param extracted_spectra: Output extracted spectra
            #   :param features: Output feature map
            #   :param compute_features: If True, populate features

        void extractSpectra(
            MSExperiment & experiment,
            TargetedExperiment & targeted_exp,
            libcpp_vector[MSSpectrum] & extracted_spectra
        ) except + nogil
            # wrap-doc:
            #   Complete extraction pipeline (without features)
            #
            #   :param experiment: Input MS experiment
            #   :param targeted_exp: Target transition list
            #   :param extracted_spectra: Output extracted spectra

        void extractSpectra(
            MSExperiment & experiment,
            FeatureMap & ms2_features,
            libcpp_vector[MSSpectrum] & extracted_spectra
        ) except + nogil
            # wrap-doc:
            #   Extract spectra using MS2 features
            #
            #   :param experiment: Input MS experiment
            #   :param ms2_features: MS2 feature map
            #   :param extracted_spectra: Output extracted spectra

        void constructTransitionsList(
            FeatureMap & ms1_features,
            FeatureMap & ms2_features,
            TargetedExperiment & t_exp
        ) except + nogil
            # wrap-doc:
            #   Construct a transition list from MS1 and MS2 features
            #
            #   :param ms1_features: MS1 feature map
            #   :param ms2_features: MS2 feature map
            #   :param t_exp: Output TargetedExperiment

        void storeSpectraMSP(const String & filename, MSExperiment & experiment) except + nogil
            # wrap-doc:
            #   Store spectra in MSP format
            #
            #   :param filename: Output filename
            #   :param experiment: Experiment to store

        void mergeFeatures(FeatureMap & fmap_input, FeatureMap & fmap_output) except + nogil
            # wrap-doc:
            #   Merge features from input to output map
            #
            #   :param fmap_input: Input feature map
            #   :param fmap_output: Output feature map (modified in place)

        # Note: matchSpectrum, targetedMatching, and untargetedMatching methods
        # that use TSE_Comparator are not wrapped due to type conversion issues.
        # Use the extractSpectra pipeline methods instead.


cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/TargetedSpectraExtractor.h>" namespace "OpenMS::TargetedSpectraExtractor":

    cdef cppclass TSE_Match "OpenMS::TargetedSpectraExtractor::Match":
        # wrap-doc:
        #   Structure holding a matched spectrum and its score

        TSE_Match() except + nogil
        TSE_Match(TSE_Match &) except + nogil
        TSE_Match(MSSpectrum & spectrum, double score) except + nogil

        MSSpectrum spectrum
        double score

    # Note: TSE_Comparator and TSE_BinnedSpectrumComparator are not wrapped
    # due to type conversion issues with vector<pair<Size, double>> parameters.
    # The matchSpectrum and targetedMatching/untargetedMatching methods that
    # use these comparators are also not wrapped for this reason.
