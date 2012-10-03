// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hannes Roest, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <cmath> // for "exp"
#include <ctime> // for "time" (random number seed)
#include <iostream> // for "cout"
#include <limits> // for "infinity"
#include <numeric> // for "accumulate"
#include <boost/bimap.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/ALGO/Scoring.h"
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>

using namespace OpenMS;
using namespace std;


//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page TOPP_OpenSwathConfidenceScoring OpenSwathConfidenceScoring

    @brief Computes confidence scores for OpenSwath results.

    <CENTER>
        <table>
            <tr>
                <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
                <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ OpenSwathConfidenceScoring \f$ \longrightarrow \f$</td>
                <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
            </tr>
            <tr>
                <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_OpenSwathAnalyzer </td>
                <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_OpenSwathFeatureXMLToTSV </td>
            </tr>
        </table>
    </CENTER>

		This is an implementation of the SRM scoring algorithm described in:

		Malmstroem, L.; Malmstroem, J.; Selevsek, N.; Rosenberger, G. & Aebersold, R.:\n
		<a href="http://dx.doi.org/10.1021/pr200844d">Automated workflow for large-scale selected reaction monitoring experiments.</a>\n
		J. Proteome Res., 2012, 11, 1644-1653

		It has been adapted for the scoring of OpenSwath results.

		The algorithm compares SRM/MRM features (peak groups) to assays and computes scores for the agreements. Every feature is compared not only to the "true" assay that was used to acquire the corresponding ion chromatograms, but also to a number (parameter @p decoys) of unrelated - but real - assays selected at random from the assay library (parameter @p lib). This serves to establish a background distribution of scores, against which the significance of the "true" score can be evaluated. The final confidence value of a feature is the local false discovery rate (FDR), calculated as the fraction of decoy assays that score higher than the "true" assay against the feature. In the output feature map, every feature is annotated with its local FDR in the meta value "local_FDR" (a "userParam" element in the featureXML), and its overall quality is set to "1 - local_FDR".

		The agreement of a feature and an assay is assessed based on the difference in retention time (RT) and on the deviation of relative transition intensities. The score @e S is computed using a binomial generalized linear model (GLM) of the form:

		@f[
		S = \frac{1}{1 + \exp(-(a + b \cdot \Delta_{RT}^2 + c \cdot d_{int}))}
		@f]

		The meanings of the model terms are as follows:

		@f$ \Delta_{RT} @f$: Observed retention times are first mapped to the scale of the assays (parameter @p trafo), then all RTs are scaled to the range 0 to 100 (based on the lowest/highest RT in the assay library). @f$ \Delta_{RT} @f$ is the absolute difference of the scaled RTs; note that this is squared in the scoring model.

		@f$ d_{int} @f$: To compute the intensity distance, the @e n (advanced parameter @p transitions) most intensive transitions of the feature are selected. For comparing against the "true" assay, the same transitions are considered; otherwise, the same number of most intensive transitions from the decoy assay. Transition intensities are scaled to a total of 1 per feature/assay and are ordered by the product (Q3) m/z value. Then the Manhattan distance of the intensity vectors is calculated (Malmstroem et. al used the RMSD instead, which has been replaced here to be independent of the number of transitions).

		@f$ a, b, c @f$: Model coefficients, stored in the advanced parameters @p GLM:intercept, @p GLM:delta_rt, and @p GLM:dist_int. The default values were estimated based on the training dataset used in the Malmstroem et al. study, reprocessed with the OpenSwath pipeline.

		In addition to the local FDRs, the scores of features against their "true" assays are recorded in the output - in the meta value "GLM_score" of the respective feature.


    <B>The command line parameters of this tool are:</B>
    @verbinclude TOPP_OpenSwathConfidenceScoring.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_OpenSwathConfidenceScoring.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPOpenSwathConfidenceScoring :
  public TOPPBase
{
public:

	/// Constructor
  TOPPOpenSwathConfidenceScoring() :
    TOPPBase("OpenSwathConfidenceScoring", 
						 "Compute confidence scores for OpenSwath results"), generator_(),
		rand_gen_(generator_, boost::uniform_int<>())
	{
		if (!test_mode_) rand_gen_.engine().seed(time(0)); // seed with current time
	}

protected:

	/// Mapping: Q3 m/z <-> transition intensity (maybe not unique!)
	typedef boost::bimap<DoubleReal, boost::bimaps::multiset_of<DoubleReal> > 
	BimapType;

	/// Binomial GLM
	struct
	{
		DoubleReal intercept;
		DoubleReal rt_coef;
		DoubleReal int_coef;

		DoubleReal operator()(DoubleReal diff_rt, DoubleReal dist_int)
		{
			DoubleReal lm = intercept + rt_coef * diff_rt * diff_rt + 
				int_coef * dist_int;
			return 1.0 / (1.0 + exp(-lm));
		}
	} glm_;

	/// Helper for RT normalization (range 0-100)
	struct
	{
		DoubleReal min_rt;
		DoubleReal max_rt;
		
		DoubleReal operator()(DoubleReal rt)
		{
			return (rt - min_rt) / (max_rt - min_rt) * 100;
		}
	} rt_norm_;

	TargetedExperiment library_; // assay library

	IntList decoy_index_; // indexes of assays to use as decoys

	Size n_decoys_; // number of decoys to use (per feature/true assay)

	Map<String, IntList> transition_map_; // assay (ID) -> transitions (indexes)

	Size n_transitions_; // number of transitions to consider

	/// RT transformation to map measured RTs to assay RTs
	TransformationDescription rt_trafo_;

	boost::mt19937 generator_; // random number generation engine

	/// Random number generator (must be initialized in init. list of c'tor!)
	boost::variate_generator<boost::mt19937&, boost::uniform_int<> > rand_gen_;

	/// Docu in base class
	void registerOptionsAndFlags_()
	{
		registerInputFile_("in", "<file>", "", "Input file (OpenSwath results)");
		setValidFormats_("in", StringList::create("featureXML"));
		registerInputFile_("lib", "<file>", "", "Assay library");
		setValidFormats_("lib", StringList::create("traML"));
		registerInputFile_("trafo", "<file>", "", "Retention time transformation");
		setValidFormats_("trafo", StringList::create("trafoXML"));
		registerOutputFile_("out", "<file>", "", 
												"Output file (results with confidence scores)");
		setValidFormats_("out", StringList::create("featureXML"));
		registerIntOption_("decoys", "<number>", 1000, "Number of decoy assays to "
											 "select from the library for every true assay (0 for "
											 "\"all\")", false);
		setMinInt_("decoys", 0);
		registerIntOption_("transitions", "<number>", 6, "Number of transitions "
											 "per feature to consider (highest intensities first; "
											 "0 for \"all\")", false);
		setMinInt_("transitions", 0);

		registerTOPPSubsection_("GLM",
														"Parameters of the binomial GLM");
		registerDoubleOption_("GLM:intercept", "<value>", 3.87333466, 
													"Intercept term", false, true);
		registerDoubleOption_("GLM:delta_rt", "<value>", -0.02898629, "Coefficient "
													"of retention time difference", false, true);
		registerDoubleOption_("GLM:dist_int", "<value>", -7.75880768,
													"Coefficient of intensity distance", false, true);
	}

	/// Randomize the list of decoy indexes
	void chooseDecoys_()
	{
		if (n_decoys_ == 0) return; // list is already initialized
		// somewhat inefficient to shuffle the whole list when we only need a random
		// sample, but easy to do...
		random_shuffle(decoy_index_.begin(), decoy_index_.end(), rand_gen_);
	}

	// DoubleReal rmsd_(DoubleList x, DoubleList y)
  // {
	// 	DoubleReal sum_of_squares = 0;
	// 	for (Size i = 0; i < x.size(); i++)
	// 	{
	// 		DoubleReal diff = x[i] - y[i];
	// 		sum_of_squares += diff * diff;
	// 	}
	// 	return sqrt(sum_of_squares / x.size());
	// }

	/// Manhattan distance
	DoubleReal manhattanDist_(DoubleList x, DoubleList y)
	{
		DoubleReal sum = 0;
		for (Size i = 0; i < x.size(); ++i)
		{
			sum += fabs(x[i] - y[i]);
		}
		return sum;
	}

	/// Get the retention time of an assay
	DoubleReal getAssayRT_(const TargetedExperiment::Peptide& assay,
												 const String& cv_accession = "MS:1000896")
	{
		String value = assay.rts[0].getCVTerms()[cv_accession][0].getValue();
		return value.toDouble();
	}

	/// Extract the @p n_transitions highest intensities from @p intensity_map,
	/// store them in @p intensities
	void extractIntensities_(BimapType& intensity_map, Size n_transitions,
													 DoubleList& intensities)
	{
		// keep only as many transitions as needed, remove those with lowest
		// intensities:
		if (n_transitions > 0)
		{
			// use "Int" instead of "Size" to prevent overflows:
			Int diff = intensity_map.size() - n_transitions;
			for (Size i = 0; Int(i) < diff; ++i)
			{
				intensity_map.right.erase(intensity_map.right.begin());
			}
		}
		// fill output list ordered by m/z:
		intensities.clear();
		for (BimapType::left_map::iterator int_it = intensity_map.left.begin();
				 int_it != intensity_map.left.end(); ++int_it)
		{
			intensities << max(0.0, int_it->second); // missing values might be "-1"
		}
	}

	/// Score the assay @p assay against feature data (@p feature_rt,
	/// @p feature_intensities), optionally using only the specified transitions
	/// (@p transition_ids)
	DoubleReal scoreAssay_(const TargetedExperiment::Peptide& assay, 
												 DoubleReal feature_rt, DoubleList& feature_intensities,
												 const set<String>& transition_ids = set<String>())
	{
		// compute RT difference:
		DoubleReal assay_rt = rt_norm_(getAssayRT_(assay));
		DoubleReal diff_rt = assay_rt - feature_rt;

		// collect transition intensities:
		BimapType intensity_map;
		for (IntList::iterator trans_it = transition_map_[assay.id].begin();
				 trans_it != transition_map_[assay.id].end(); ++trans_it)
		{
			const ReactionMonitoringTransition& transition = 
				library_.getTransitions()[*trans_it];
			// for the "true" assay, we need to choose the same transitions as for the
			// feature:
			if (!transition_ids.empty() && 
					(transition_ids.count(transition.getNativeID()) == 0)) continue;
			// seems like Boost's Bimap doesn't support "operator[]"...
			intensity_map.left.insert(make_pair(transition.getProductMZ(), 
																					transition.getLibraryIntensity()));
		}
		DoubleList assay_intensities;
		extractIntensities_(intensity_map, feature_intensities.size(), 
												assay_intensities);

		// compute intensity distance:
		OpenSwath::Scoring::normalize_sum(&feature_intensities[0],
																			feature_intensities.size());
		OpenSwath::Scoring::normalize_sum(&assay_intensities[0],
																			assay_intensities.size());
		DoubleReal dist_int = manhattanDist_(feature_intensities, 
																				 assay_intensities);

		DoubleReal score = glm_(diff_rt, dist_int);

		LOG_DEBUG << "\ndelta_RT:  " << fabs(diff_rt)
							<< "\ndist_int:  " << dist_int
							<< "\nGLM_score: " << score << endl;

		return score;
	}

	/// Score a feature
	void scoreFeature_(Feature& feature)
	{
		// extract predictors from feature:
		DoubleReal feature_rt = rt_norm_(rt_trafo_.apply(feature.getRT()));
		BimapType intensity_map;
		// for the "true" assay, we need to make sure we compare based on the same
		// transitions, so keep track of them:
		Map<DoubleReal, String> trans_id_map; // Q3 m/z -> transition ID
		for (vector<Feature>::iterator sub_it = feature.getSubordinates().begin();
				 sub_it != feature.getSubordinates().end(); ++sub_it)
		{
			// seems like Boost's Bimap doesn't support "operator[]"...
			intensity_map.left.insert(make_pair(sub_it->getMZ(), 
																					sub_it->getIntensity()));
			trans_id_map[sub_it->getMZ()] = sub_it->getMetaValue("native_id");
		}
		DoubleList feature_intensities;
		extractIntensities_(intensity_map, n_transitions_, feature_intensities);
		if ((n_transitions_ > 0) && (feature_intensities.size() < n_transitions_))
		{
			LOG_WARN << "Warning: Feature '" << feature.getUniqueId() 
							 << "' contains only " << feature_intensities.size()
							 << " transitions." << endl;
		}
		// "intensity_map" now only contains the transitions we need later:
		set<String> transition_ids;
		for (BimapType::left_map::iterator int_it = intensity_map.left.begin();
				 int_it != intensity_map.left.end(); ++int_it)
		{
			transition_ids.insert(trans_id_map[int_it->first]);
		}

		DoubleList scores; // "true" score is in "scores[0]", decoy scores follow

		// compare to "true" assay:
		String true_id = feature.getMetaValue("PeptideRef");
		LOG_DEBUG << "True assay (ID '" << true_id << "')" << endl;
		scores << scoreAssay_(library_.getPeptideByRef(true_id), feature_rt, 
													feature_intensities, transition_ids);

		// compare to decoy assays:
		chooseDecoys_();
		Size counter = 0;
		for (IntList::iterator decoy_it = decoy_index_.begin(); 
				 decoy_it != decoy_index_.end(); ++decoy_it)
		{
			const TargetedExperiment::Peptide& decoy_assay = 
				library_.getPeptides()[*decoy_it];

			// skip the "true" assay and assays with too few transitions:
			// TODO: maybe add an option to include assays with too few transitions?
			if ((decoy_assay.id == true_id) || 
					(transition_map_[decoy_assay.id].size() < feature_intensities.size()))
			{
				continue;
			}
			LOG_DEBUG << "Decoy assay " << scores.size() << " (ID '" << decoy_assay.id
								<< "')" << endl;

			scores << scoreAssay_(decoy_assay, feature_rt, feature_intensities);

			if ((n_decoys_ > 0) && (++counter >= n_decoys_)) break; // enough decoys
		}
		
		Size n_scores = scores.size();
		if (n_scores - 1 < n_decoys_)
		{
			LOG_WARN << "Warning: Feature '" << feature.getUniqueId() 
							 << "': Couldn't find enough decoy assays with at least "
							 << feature_intensities.size() << " transitions. "
							 << "Scoring based on " << n_scores - 1 << " decoys." << endl;
		}
		// TODO: this warning may trigger for every feature and get annoying
		if ((n_decoys_ == 0) && (n_scores < library_.getPeptides().size()))
		{
			LOG_WARN << "Warning: Feature '" << feature.getUniqueId() 
							 << "': Skipped some decoy assays with fewer than " 
							 << feature_intensities.size() << " transitions. "
							 << "Scoring based on " << n_scores - 1 << " decoys." << endl;
		}

		// count decoy scores that are greater than the "true" score:
		counter = 0;
		for (DoubleList::iterator it = ++scores.begin(); it != scores.end(); ++it)
		{
			if (*it > scores[0]) counter++;
		}

		// annotate feature:
		feature.setMetaValue("GLM_score", scores[0]);
		DoubleReal local_fdr = counter / (n_scores - 1.0);
		feature.setMetaValue("local_FDR", local_fdr);
		feature.setOverallQuality(1.0 - local_fdr);
	}

	/// Docu in base class
	ExitCodes main_(int, const char**)
	{
		if (debug_level_ > 0) Log_debug.insert(cout);

		LOG_DEBUG << "Reading parameters..." << endl;
    String in = getStringOption_("in");
		String lib = getStringOption_("lib");
		String trafo = getStringOption_("trafo");
    String out = getStringOption_("out");
		n_decoys_ = getIntOption_("decoys");
		n_transitions_ = getIntOption_("transitions");

		glm_.intercept = getDoubleOption_("GLM:intercept");
		glm_.rt_coef = getDoubleOption_("GLM:delta_rt");
		glm_.int_coef = getDoubleOption_("GLM:dist_int");
		
		LOG_DEBUG << "Loading input files..." << endl;
		FeatureMap<> features;
		FeatureXMLFile().load(in, features);
		TraMLFile().load(lib, library_);
		TransformationXMLFile().load(trafo, rt_trafo_);
		if (rt_trafo_.getModelType() == "none") // fit a linear model now
		{
			rt_trafo_.fitModel("linear");
		}

		// are there enough assays in the library?
		Size n_assays = library_.getPeptides().size();
		if (n_assays < 2)
		{
			LOG_FATAL_ERROR << "Error: Not enough assays in the library!" << endl;
			return INCOMPATIBLE_INPUT_DATA; // right exit code for this situation?
		}
		if (n_assays - 1 < n_decoys_)
		{
			LOG_WARN << "Warning: Parameter 'decoys' (" << n_decoys_ 
							 << ") is higher than the number of unrelated assays in the "
							 << "library (" << n_assays - 1 << "). "
							 << "Using all unrelated assays as decoys." << endl;
		}
		if (n_assays - 1 <= n_decoys_) n_decoys_ = 0; // use all available assays

		decoy_index_.resize(n_assays);
		for (Size i = 0; i < n_assays; ++i) decoy_index_[i] = i;

		// build mapping between assays and transitions:
		LOG_DEBUG << "Building transition map..." << endl;
		for (Size i = 0; i < library_.getTransitions().size(); ++i)
		{
			const String& ref = library_.getTransitions()[i].getPeptideRef();
			transition_map_[ref].push_back(i);
		}
		// find min./max. RT in the library:
		LOG_DEBUG << "Determining retention time range..." << endl;
		rt_norm_.min_rt = numeric_limits<double>::infinity();
		rt_norm_.max_rt = -numeric_limits<double>::infinity();
		for (vector<TargetedExperiment::Peptide>::const_iterator it = 
					 library_.getPeptides().begin(); it != library_.getPeptides().end();
				 ++it)
		{
			DoubleReal current_rt = getAssayRT_(*it);
			rt_norm_.min_rt = min(rt_norm_.min_rt, current_rt);
			rt_norm_.max_rt = max(rt_norm_.max_rt, current_rt);
		}

		// log scoring progress:
		ProgressLogger progress;
		progress.setLogType(log_type_);
		LOG_DEBUG << "Scoring features..." << endl;
		progress.startProgress(0, features.size(), "scoring features");

		for (FeatureMap<>::Iterator feat_it = features.begin(); 
				 feat_it != features.end(); ++feat_it)
		{
			LOG_DEBUG << "Feature " << feat_it - features.begin() + 1 
								<< " (ID '" << feat_it->getUniqueId() << "')"<< endl;
			scoreFeature_(*feat_it);
			progress.setProgress(feat_it - features.begin());
		}
		
		progress.endProgress();
		LOG_DEBUG << "Storing results..." << endl;
		addDataProcessing_(features, 
											 getProcessingInfo_(DataProcessing::DATA_PROCESSING));
		FeatureXMLFile().store(out, features);

		return EXECUTION_OK;
	}

};


int main(int argc, const char** argv)
{
  TOPPOpenSwathConfidenceScoring t;
  return t.main(argc, argv);
}

/// @endcond
