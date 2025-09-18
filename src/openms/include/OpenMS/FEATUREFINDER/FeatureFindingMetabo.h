// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Erhan Kenar, Holger Franken $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/CHEMISTRY/Element.h>

#include <vector>

struct svm_model;

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
    FeatureHypothesis() = default;

    /// default destructor
    ~FeatureHypothesis() = default;

    /// copy constructor
    FeatureHypothesis(const FeatureHypothesis&) = default;

    /// assignment operator
    FeatureHypothesis& operator=(const FeatureHypothesis& rhs) = default;

    // getter & setter
    Size getSize() const;

    String getLabel() const;

    std::vector<String> getLabels() const;

    double getScore() const;

    void setScore(const double& score);

    SignedSize getCharge() const;

    void setCharge(const SignedSize& ch);

    std::vector<double> getAllIntensities(bool smoothed = false) const;

    std::vector<double> getAllCentroidMZ() const;

    std::vector<double> getAllCentroidRT() const;

    std::vector<double> getIsotopeDistances() const;

    double getCentroidMZ() const;

    double getCentroidRT() const;

    double getFWHM() const;

    /// addMassTrace
    void addMassTrace(const MassTrace&);
    double getMonoisotopicFeatureIntensity(bool) const;
    double getSummedFeatureIntensity(bool) const;

    /// return highest apex of all isotope traces
    double getMaxIntensity(bool smoothed = false) const;

    Size getNumFeatPoints() const;
    std::vector<ConvexHull2D> getConvexHulls() const;
    std::vector< OpenMS::MSChromatogram > getChromatograms(UInt64 feature_id) const;

private:

    // pointers of MassTraces contained in isotopic pattern
    std::vector<const MassTrace*> iso_pattern_;

    double feat_score_{};

    SignedSize charge_{};
  };

  class OPENMS_DLLAPI CmpMassTraceByMZ
  {
public:

    bool operator()(const MassTrace& x, const MassTrace& y) const
    {
      return x.getCentroidMZ() < y.getCentroidMZ();
    }

  };

  class OPENMS_DLLAPI CmpHypothesesByScore
  {
public:

    bool operator()(const FeatureHypothesis& x, const FeatureHypothesis& y) const
    {
      return x.getScore() > y.getScore();
    }

  };

  /**
   * @brief Internal structure to store a lower and upper bound of an m/z range
   */
  struct OPENMS_DLLAPI  Range
{
  double left_boundary;
  double right_boundary;
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
    ~FeatureFindingMetabo() override;

    /// main method of FeatureFindingMetabo
    void run(std::vector<MassTrace>& input_mtraces, FeatureMap& output_featmap, std::vector<std::vector< OpenMS::MSChromatogram > >& output_chromatograms);

protected:
    void updateMembers_() override;

private:
    /**
     * @brief parses a string of element symbols into a vector of Elements
     * @param elements_string string of element symbols without whitespaces or commas. e.g. CHNOPSCl
     * @return vector of Elements
     */
    std::vector<const Element*> elementsFromString_(const std::string& elements_string) const;
    /**
     * Calculate the maximal and minimal mass defects of isotopes for a given set of elements.
     *
     * @param alphabet   chemical alphabet (elements which are expected to be present)
     * @param peakOffset integer distance between isotope peak and monoisotopic peak (minimum: 1)
     * @return an interval which should contain the isotopic peak. This interval is relative to the monoisotopic peak.
     */
    Range getTheoreticIsotopicMassWindow_(const std::vector<Element const *>& alphabet, int peakOffset) const;

    /** @brief Computes the cosine similarity between two vectors
     *
     * The cosine similarity (or cosine distance) is the cosine of the angle
     * between two vectors or the normalized dot product of two vectors.
     *
     * See also https://en.wikipedia.org/wiki/Cosine_similarity
     *
    */
    double computeCosineSim_(const std::vector<double>&, const std::vector<double>&) const;

    /** @brief Compare intensities of feature hypothesis with model 
     *
     * Use a pre-trained SVM model to evaluate the intensity distribution of a
     * given feature hypothesis. The model is trained on the monoisotopic and
     * the first tree isotopic traces of each feature and uses the scaled
     * ratios between the traces as input.
     *
     * Reference: Kenar et al., doi: 10.1074/mcp.M113.031278
     *
     * @param feat_hypo A feature hypotheses containing mass traces
     * @return 0 for 'no'; 1 for 'yes'; -1 if only a single mass trace exists
    */
    int isLegalIsotopePattern_(const FeatureHypothesis& feat_hypo) const;

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
     * An alternative scoring was added which test if isotope m/z distances lie in an expected m/z window.
     * This window is computed from a given set of elements.
     *
    */
    double scoreMZ_(const MassTrace &, const MassTrace &, Size isotopic_position, Size charge, Range isotope_window) const;

    /**
     * @brief score isotope m/z distance based on the expected m/z distances using C13-C12 or Kenar method
     * @param iso_pos
     * @param charge
     * @param diff_mz
     * @param mt_variances
     * @return
     */
    double scoreMZByExpectedMean_(Size iso_pos, Size charge, const double diff_mz, double mt_variances) const;

    /**
     * @brief score isotope m/z distance based on an expected isotope window which was calculated from a set of expected elements
     * @param charge
     * @param diff_mz
     * @param mt_variances m/z variance between the two mass traces which are compared
     * @param isotope_window
     * @return
     */
    double scoreMZByExpectedRange_(Size charge, const double diff_mz, double mt_variances, Range isotope_window) const;

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
    double computeAveragineSimScore_(const std::vector<double>& intensities, const double& molecular_weight) const;

    /** @brief Identify groupings of mass traces based on a set of reasonable candidates
     *
     * Takes a set of reasonable candidates for mass trace grouping and checks
     * all combinations of charge and isotopic positions on the candidates. It
     * is assumed that candidates[0] is the monoisotopic trace.
     *
     * The resulting possible groupings are appended to output_hypotheses.
    */
    void findLocalFeatures_(const std::vector<const MassTrace*>& candidates, double total_intensity, std::vector<FeatureHypothesis>& output_hypotheses) const;

    /// SVM parameters
    svm_model* isotope_filt_svm_ = nullptr;
    std::vector<double> svm_feat_centers_;
    std::vector<double> svm_feat_scales_;

    //unused
    //double total_intensity_;

    /// parameter stuff
    double local_rt_range_;
    double local_mz_range_;
    Size charge_lower_bound_;
    Size charge_upper_bound_;
    double chrom_fwhm_;

    bool report_summed_ints_;
    bool enable_RT_filtering_;
    String isotope_filtering_model_;
    bool use_smoothed_intensities_;
    bool report_smoothed_intensities_;

    bool use_mz_scoring_C13_;
    bool use_mz_scoring_by_element_range_;
    bool report_convex_hulls_;
    bool report_chromatograms_;

    bool remove_single_traces_;
    std::vector<const Element*> elements_;
  };

}

