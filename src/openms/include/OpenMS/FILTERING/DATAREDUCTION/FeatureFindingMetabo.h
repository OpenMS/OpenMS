// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
// $Maintainer: Erhan Kenar $
// $Authors: Erhan Kenar, Holger Franken $
// --------------------------------------------------------------------------

#ifndef OPENMS_FILTERING_DATAREDUCTION_FEATUREFINDINGMETABO_H
#define OPENMS_FILTERING_DATAREDUCTION_FEATUREFINDINGMETABO_H

#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <vector>
#include <svm.h>

namespace OpenMS
{

  /**
    @brief Internal structure used in @ref FeatureFindingMetabo that keeps
    track of a feature hypothesis (isotope group hypothesis).

    @htmlinclude OpenMS_FeatureFindingMetabo.parameters

    @ingroup Quantitation
  */
  class OPENMS_DLLAPI FeatureHypothesis
  {
public:
    /// default constructor
    FeatureHypothesis();

    /// default destructor
    ~FeatureHypothesis();

    /// copy constructor
    FeatureHypothesis(const FeatureHypothesis&);

    /// assignment operator
    FeatureHypothesis& operator=(const FeatureHypothesis& rhs);

    // getter & setter
    Size getSize() const
    {
      return iso_pattern_.size();
    }

    String getLabel() const
    {
      String label;

      if (iso_pattern_.size() > 0)
      {
        label = iso_pattern_[0]->getLabel();
      }

      for (Size i = 1; i < iso_pattern_.size(); ++i)
      {
        String tmp_str = "_" + iso_pattern_[i]->getLabel();
        label += tmp_str;
      }

      return label;
    }

    std::vector<String> getLabels() const
    {
      std::vector<String> tmp_labels;

      for (Size i = 0; i < iso_pattern_.size(); ++i)
      {
        tmp_labels.push_back(iso_pattern_[i]->getLabel());
      }

      return tmp_labels;
    }

    double getScore() const
    {
      return feat_score_;
    }

    void setScore(const double& score)
    {
      feat_score_ = score;
    }

    SignedSize getCharge() const
    {
      return charge_;
    }

    void setCharge(const SignedSize& ch)
    {
      charge_ = ch;
    }

    std::vector<double> getAllIntensities(bool smoothed = false) const
    {
      std::vector<double> tmp;
      for (Size i = 0; i < iso_pattern_.size(); ++i)
      {
        tmp.push_back(iso_pattern_[i]->getIntensity(smoothed));
      }
      return tmp;
    }

    double getCentroidMZ() const
    {
      if (iso_pattern_.empty())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,
            "FeatureHypothesis is empty, no centroid MZ!", String(iso_pattern_.size()));
      }
      return iso_pattern_[0]->getCentroidMZ();
    }

    double getCentroidRT() const
    {
      if (iso_pattern_.empty())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,
            "FeatureHypothesis is empty, no centroid RT!", String(iso_pattern_.size()));
      }
      return iso_pattern_[0]->getCentroidRT();
    }

    double getFWHM() const
    {
      if (iso_pattern_.empty())
      {
        return 0.0;
      }
      return iso_pattern_[0]->getFWHM();
    }

    /// addMassTrace
    void addMassTrace(const MassTrace&);
    double getMonoisotopicFeatureIntensity(bool) const;
    double getSummedFeatureIntensity(bool) const;

    Size getNumFeatPoints() const;
    std::vector<ConvexHull2D> getConvexHulls() const;

private:

    // pointers of MassTraces contained in isotopic pattern
    std::vector<const MassTrace*> iso_pattern_;

    double feat_score_;

    SignedSize charge_;
  };

  class OPENMS_DLLAPI CmpMassTraceByMZ
  {
public:

    bool operator()(MassTrace x, MassTrace y) const
    {
      return x.getCentroidMZ() < y.getCentroidMZ();
    }

  };

  class OPENMS_DLLAPI CmpHypothesesByScore
  {
public:

    bool operator()(FeatureHypothesis x, FeatureHypothesis y) const
    {
      return x.getScore() > y.getScore();
    }

  };

  /**
    @brief Method for the assembly of mass traces belonging to the same isotope
    pattern, i.e., that are compatible in retention times, mass-to-charge ratios,
    and isotope abundances.

    In @ref FeatureFindingMetabo, mass traces detected by the @ref
    MassTraceDetection method and afterwards split into individual
    chromatographic peaks by the @ref ElutionPeakDetection method are assembled
    to composite features if they are compatible with respect to RTs, m/z ratios,
    and isotopic intensities. To this end, feature hypotheses are formulated
    exhaustively based on the set of mass traces detected within a local RT and
    m/z region. These feature hypotheses are scored by their similarity to real
    metabolite isotope patterns. The score is derived from independent models for
    retention time shifts and m/z differences between isotopic mass traces.
    Hypotheses with correct or false isotopic abundances are distinguished by a
    SVM model. Mass traces that could not be assembled or low-intensity
    metabolites with only a monoisotopic mass trace to observe are left in the
    resulting @ref FeatureMap as singletons with the undefined charge state of 0.

    Reference: Kenar et al., doi: 10.1074/mcp.M113.031278

    @htmlinclude OpenMS_FeatureFindingMetabo.parameters

    @ingroup Quantitation
  */
  class OPENMS_DLLAPI FeatureFindingMetabo :
    public DefaultParamHandler,
    public ProgressLogger
  {
public:
    /// Default constructor
    FeatureFindingMetabo();

    /// Default destructor
    virtual ~FeatureFindingMetabo();

    /// main method of FeatureFindingMetabo
    void run(std::vector<MassTrace>&, FeatureMap&);

protected:
    virtual void updateMembers_();

private:
    /// private member functions

    /** @brief Computes the cosine similarity between two vectors
     *
     * The cosine similarity (or cosine distance) is the cosine of the angle
     * between two vectors or the normalized dot product of two vectors.
     *
     * See also https://en.wikipedia.org/wiki/Cosine_similarity
     *
    */
    double computeCosineSim_(const std::vector<double>&, const std::vector<double>&) const;

    /// unused function ???
    /// TODO: remove
    double computeOLSCoeff_(const std::vector<double>&, const std::vector<double>&) const;

    /// Only used in debugging output ???  -> seems old model with 7 traces
    /// TODO: remove
    bool isLegalIsotopePattern_(const FeatureHypothesis &) const;

    /** @brief Compare intensities of feature hypothesis with model 
     *
     * Use a pre-trained SVM model to evaluate the intensity distribution of a
     * given feature hypothesis. The model is trained on the monoisotopic and
     * the first tree isotopic traces of each feature and uses the scaled
     * ratios between the traces as input.
     *
     * Reference: Kenar et al., doi: 10.1074/mcp.M113.031278
     *
    */
    bool isLegalIsotopePattern2_(const FeatureHypothesis &) const;

    void loadIsotopeModel_(const String&);

    /** @brief Perform mass to charge scoring of two multiple mass traces
     *
     * Scores two mass traces based on the m/z and the hypothesis that one
     * trace is an isotopic trace of the other one. The isotopic position
     * (which trace it is) and the charge for the hypothesis are given as
     * additional parameters. 
     * The scoring is described in Kenar et al., and is based on a random
     * sample of 115 000 compounds drawn from a comprehensive set of 24 million
     * putative sum formulas, of which the isotopic distribution was accurately
     * calculated. Thus, a theoretical mu and sigma are calculated as:
     *
     * mu = 1.000857 * j + 0.001091 u
     * sigma = 0.0016633 j * 0.0004751
     *
     * where j is the isotopic peak considered. A similarity score based on
     * agreement with the model is then computed.
     *
     * Reference: Kenar et al., doi: 10.1074/mcp.M113.031278
     *
    */
    double scoreMZ_(const MassTrace &, const MassTrace &, Size isotopic_position, Size charge) const;

    /// Not used any more ???  -> seems to be old model 
    /// TODO: remove
    double scoreMZ2_(const MassTrace &, const MassTrace &, Size isotopic_position, Size charge) const;

    /** @brief Perform retention time scoring of two multiple mass traces
     *
     * Computes the similarity of the two peak shapes using cosine similarity
     * (see computeCosineSim_) if some conditions are fulfilled. Mainly the
     * overlap between the two peaks at FHWM needs to exceed a certain
     * threshold. The threshold is set at 0.7 (i.e. 70 % overlap) as also
     * described in Kenar et al.
     *
     * @note this only works for equally sampled mass traces, e.g. they need to
     * come from the same map (not for SRM measurements for example).
    */
    double scoreRT_(const MassTrace&, const MassTrace&) const;

    /** @brief Perform intensity scoring using the averagine model (for peptides only)
     *
     * Compare the isotopic intensity distribution with the theoretical one
     * expected for peptides, using the averagine model. Compute the cosine
     * similarity between the two values.
    */
    double computeAveragineSimScore_(const std::vector<double>& intensities, const double& molecular_weight);

    /** @brief Identify groupings of mass traces based on a set of reasonable candidates
     *
     * Takes a set of reasonable candidates for mass trace grouping and checks
     * all combinations of charge and isotopic positions on the candidates. It
     * is assumed that candidates[0] is the monoisotopic trace.
     *
     * The resulting possible groupings are appended to output_hypotheses.
    */
    void findLocalFeatures_(const std::vector<const MassTrace*>& candidates, std::vector<FeatureHypothesis>& output_hypotheses);

    /// SVM parameters
    svm_model* isotope_filt_svm_;
    std::vector<double> svm_feat_centers_;
    std::vector<double> svm_feat_scales_;

    double total_intensity_;

    /// parameter stuff
    double local_rt_range_;
    double local_mz_range_;
    Size charge_lower_bound_;
    Size charge_upper_bound_;
    //double mass_error_ppm_;
    double chrom_fwhm_;

    bool report_summed_ints_;
    bool enable_RT_filtering_;
    bool disable_isotope_filtering_;
    String isotope_model_;
    String metabo_iso_noisemodel_;
    bool use_smoothed_intensities_;

  };

}

#endif // OPENMS_FILTERING_DATAREDUCTION_FEATUREFINDINGMETABO_H
