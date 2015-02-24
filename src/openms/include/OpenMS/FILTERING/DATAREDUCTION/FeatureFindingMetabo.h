// --------------------------------------------------------------------------
//           OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//  notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//  notice, this list of conditions and the following disclaimer in the
//  documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//  may be used to endorse or promote products derived from this software
//  without specific prior written permission.
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
  @brief Method for the assembly of mass traces belonging to the same isotope pattern, i.e., that are compatible in retention times, mass-to-charge ratios, and isotope abundances.

  In @ref FeatureFindingMetabo, mass traces detected by the @ref MassTraceDetection method and afterwards split into individual chromatographic peaks by the
  @ref ElutionPeakDetection method are assembled to composite features if they are compatible with respect to RTs, m/z ratios, and isotopic intensities. To this end,
  feature hypotheses are formulated exhaustively based on the set of mass traces detected within a local RT and m/z region. These feature hypotheses are scored by their similarity to
  real metabolite isotope patterns. The score is derived from independent models for retention time shifts and m/z differences between isotopic mass traces.
  Hypotheses with correct or false isotopic abundances are distinguished by a SVM model. Mass traces that could not be assembled or low-intensity metabolites with only a
  monoisotopic mass trace to observe are left in the resulting @ref FeatureMap as singletons with the undefined charge state of 0.

  @htmlinclude OpenMS_FeatureFindingMetabo.parameters

  @ingroup Quantitation
  */

class OPENMS_DLLAPI CmpMassTraceByMZ
{
public:

  bool operator()(MassTrace x, MassTrace y) const
  {
    return x.getCentroidMZ() < y.getCentroidMZ();
  }

};

class OPENMS_DLLAPI FeatureHypothesis
{
public:
  /// default constructor
  FeatureHypothesis();

  /// default destructor
  ~FeatureHypothesis();

  /// copy constructor
  FeatureHypothesis(const FeatureHypothesis &);

  /// assignment operator
  FeatureHypothesis & operator=(const FeatureHypothesis & rhs);

  // getter & setter
  Size getSize() const;

  String getLabel() const;

  /// Collect all labels of the isotope patterns
  std::vector<String> getLabels() const;

  double getScore() const;

  void setScore(const double & score);

  SignedSize getCharge() const;

  void setCharge(const SignedSize & ch);

  std::vector<double> getAllIntensities(bool smoothed = false);

  double getCentroidMZ() const;

  double getCentroidRT() const;

  double getFWHM(bool use_smoothed_ints = false) const;

  /// addMassTrace
  void addMassTrace(MassTrace &);
  double getMonoisotopicFeatureIntensity(bool);
  double getSummedFeatureIntensity(bool);

  Size getNumFeatPoints() const;
  std::vector<ConvexHull2D> getConvexHulls() const;

private:
  // pointers of MassTraces contained in isotopic pattern
  std::vector<MassTrace *> iso_pattern_;
  double feat_score_;

  SignedSize charge_;

};

class OPENMS_DLLAPI CmpHypothesesByScore
{
public:

  bool operator()(FeatureHypothesis x, FeatureHypothesis y) const
  {
    return x.getScore() > y.getScore();
  }

};

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
  void run(std::vector<MassTrace> &, FeatureMap &);

protected:

  virtual void updateMembers_();

private:

  /// private member functions
  double computeOLSCoeff_(const std::vector<double> &, const std::vector<double> &);
  double computeCosineSim_(const std::vector<double> &, const std::vector<double> &);

  svm_model * isotope_filt_svm_;
  std::vector<double> svm_feat_centers_;
  std::vector<double> svm_feat_scales_;
  // bool isLegalIsotopePattern_(FeatureHypothesis &);
  bool isLegalIsotopePattern2_(FeatureHypothesis &);

  // bool isLegalAveraginePattern(FeatureHypothesis&);
  void loadIsotopeModel_(const String&);

  double total_intensity_;

  double scoreMZ_(const MassTrace &, const MassTrace &, Size, Size);
  // double scoreMZ2_(const MassTrace &, const MassTrace &, Size, Size);
  double scoreRT_(const MassTrace &, const MassTrace &);

  double computeAveragineSimScore_(const std::vector<double> &, const double &);

  // double scoreTraceSim_(MassTrace, MassTrace);
  // double scoreIntRatio_(double, double, Size);
  void findLocalFeatures_(std::vector<MassTrace *> &, std::vector<FeatureHypothesis> &);

  /// parameter stuff
  double local_rt_range_;
  double local_mz_range_;
  Size charge_lower_bound_;
  Size charge_upper_bound_;
  //double mass_error_ppm_;
  double chrom_fwhm_;

  bool report_summed_ints_;
  bool disable_isotope_filtering_;
  String isotope_model_;
  String metabo_iso_noisemodel_;
  bool use_smoothed_intensities_;

};


}




#endif // OPENMS_FILTERING_DATAREDUCTION_FEATUREFINDINGMETABO_H
