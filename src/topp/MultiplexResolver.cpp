// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DBoundingBox.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/FeatureHandle.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/MATH/STATISTICS/LinearRegression.h>
#include <OpenMS/KERNEL/RangeUtils.h>
#include <OpenMS/KERNEL/ChromatogramTools.h>
#include <OpenMS/FORMAT/MzQuantMLFile.h>
#include <OpenMS/FORMAT/PeakTypeEstimator.h>

#include <OpenMS/METADATA/MSQuantifications.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/MATH/STATISTICS/LinearRegression.h>

#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/MzQuantMLFile.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexDeltaMasses.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexDeltaMassesGenerator.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexIsotopicPeakPattern.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilteringCentroided.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilteringProfile.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexClustering.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <OpenMS/DATASTRUCTURES/DBoundingBox.h>

//Contrib includes
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/classification.hpp>

#include <QDir>

//std includes
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <limits>
#include <locale>
#include <iomanip>

using namespace std;
using namespace OpenMS;
using namespace boost::math;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_MultiplexResolver MultiplexResolver

  @brief Completes peptide multiplets and resolves conflicts within them.

<CENTER>
  <table>
    <tr>
      <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
      <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ MultiplexResolver \f$ \longrightarrow \f$</td>
      <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDMapper </td>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_ProteinQuantifier </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDConflictResolver </td>
    </tr>
  </table>
</CENTER>

  Tools such as FeatureFinderMultiplex can detect peptide feature multiplets in labeled experimental data. The multiplets can then be annotated with peptide sequences
  using the IDMapper tool. The MultiplexResolver tool is consolidating these results in two steps.
  - Any multiplets with conflicting quantitative and sequence information are filtered out. As example, let us consider a triple SILAC analyis. Let us assume a sequence
  "LDNLVAIFDINR(Label:13C(6)15N(4))" with a single Arg10 label is mapped to the light feature in a SILAC triplet. Either peptide feature detection or sequence information
  must be incorrect und the triplet is removed.
  - In a second step, any incomplete peptide feature groups are completed with dummy features of zero intensity. As example, let us stay with the triple SILAC analysis.
  But let us now assume the sequence "LDNLVAIFDINR(Label:13C(6)15N(4))" is mapped to the heavy partner of a peptide feature pair. This is no conflict. Medium and heavy
  peptides have been correctly detected. The MultiplexResolver adds a dummy peptide feature of zero intensity at the light position and thereby completes the triplet.

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_MultiplexResolver.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_MultiplexResolver.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPMultiplexResolver :
  public TOPPBase
{
private:

  // input and output files
  String in_;
  String out_;
  String out_conflicts_;

  // section "algorithm"
  String labels_;
  std::vector<std::vector<String> > samples_labels_;
  unsigned charge_min_;
  unsigned charge_max_;
  unsigned missed_cleavages_;
  unsigned isotopes_per_peptide_min_;
  unsigned isotopes_per_peptide_max_;
  double mass_tolerance_;

  // section "labels"
  map<String, double> label_mass_shift_;
  
public:
  TOPPMultiplexResolver() :
    TOPPBase("MultiplexResolver", "Completes peptide multiplets and resolves conflicts within them.", true),
    charge_min_(1), charge_max_(1), missed_cleavages_(0), isotopes_per_peptide_min_(1), isotopes_per_peptide_max_(1), mass_tolerance_(0.1)
  {
  }
  
  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "Peptide multiplets with assigned sequence information");
    setValidFormats_("in", ListUtils::create<String>("consensusXML"));
    registerOutputFile_("out", "<file>", "", "Complete peptide multiplets.", false);
    setValidFormats_("out", ListUtils::create<String>("consensusXML"));
    registerOutputFile_("out_conflicts", "<file>", "", "Optional output containing peptide multiplets with conflicting quant/ID information.", false, true);
    setValidFormats_("out_conflicts", ListUtils::create<String>("consensusXML"));

    registerSubsection_("algorithm", "Parameters for the algorithm.");
    registerSubsection_("labels", "Isotopic labels that can be specified in section \'algorithm:labels\'.");

  }

  // create parameters for sections (set default values and restrictions)
  Param getSubsectionDefaults_(const String& section) const
  {
    Param defaults;

    if (section == "algorithm")
    {
      defaults.setValue("labels", "[][Lys8,Arg10]", "Labels used for labelling the samples. [...] specifies the labels for a single sample. For example\n\n[][Lys8,Arg10]        ... SILAC\n[][Lys4,Arg6][Lys8,Arg10]        ... triple-SILAC\n[Dimethyl0][Dimethyl6]        ... Dimethyl\n[Dimethyl0][Dimethyl4][Dimethyl8]        ... triple Dimethyl\n[ICPL0][ICPL4][ICPL6][ICPL10]        ... ICPL");
      defaults.setValue("charge", "1:4", "Range of charge states in the sample, i.e. min charge : max charge.");
      defaults.setValue("isotopes_per_peptide", "3:6", "Range of isotopes per peptide in the sample. For example 3:6, if isotopic peptide patterns in the sample consist of either three, four, five or six isotopic peaks. ", ListUtils::create<String>("advanced"));
      defaults.setValue("missed_cleavages", 0, "Maximum number of missed cleavages due to incomplete digestion.");
      defaults.setMinInt("missed_cleavages", 0);
      defaults.setValue("mass_tolerance", 0.1, "Mass tolerance for matching the detected to the theoretical mass shifts.", ListUtils::create<String>("advanced"));
    }

    if (section == "labels")
    {
      // create labels that can be chosen in section "algorithm/labels"
      defaults.setValue("Arg6", 6.0201290268, "Label:13C(6)  |  C(-6) 13C(6)  |  unimod #188", ListUtils::create<String>("advanced"));
      defaults.setMinFloat("Arg6", 0.0);
      defaults.setValue("Arg10", 10.008268600, "Label:13C(6)15N(4)  |  C(-6) 13C(6) N(-4) 15N(4)  |  unimod #267", ListUtils::create<String>("advanced"));
      defaults.setMinFloat("Arg10", 0.0);
      defaults.setValue("Lys4", 4.0251069836, "Label:2H(4)  |  H(-4) 2H(4)  |  unimod #481", ListUtils::create<String>("advanced"));
      defaults.setMinFloat("Lys4", 0.0);
      defaults.setValue("Lys6", 6.0201290268, "Label:13C(6)  |  C(-6) 13C(6)  |  unimod #188", ListUtils::create<String>("advanced"));
      defaults.setMinFloat("Lys6", 0.0);
      defaults.setValue("Lys8", 8.0141988132, "Label:13C(6)15N(2)  |  C(-6) 13C(6) N(-2) 15N(2)  |  unimod #259", ListUtils::create<String>("advanced"));
      defaults.setMinFloat("Lys8", 0.0);
      defaults.setValue("Dimethyl0", 28.031300, "Dimethyl  |  H(4) C(2)  |  unimod #36", ListUtils::create<String>("advanced"));
      defaults.setMinFloat("Dimethyl0", 0.0);
      defaults.setValue("Dimethyl4", 32.056407, "Dimethyl:2H(4)  |  2H(4) C(2)  |  unimod #199", ListUtils::create<String>("advanced"));
      defaults.setMinFloat("Dimethyl4", 0.0);
      defaults.setValue("Dimethyl6", 34.063117, "Dimethyl:2H(4)13C(2)  |  2H(4) 13C(2)  |  unimod #510", ListUtils::create<String>("advanced"));
      defaults.setMinFloat("Dimethyl6", 0.0);
      defaults.setValue("Dimethyl8", 36.075670, "Dimethyl:2H(6)13C(2)  |  H(-2) 2H(6) 13C(2)  |  unimod #330", ListUtils::create<String>("advanced"));
      defaults.setMinFloat("Dimethyl8", 0.0);
      defaults.setValue("ICPL0", 105.021464, "ICPL  |  H(3) C(6) N O  |  unimod #365", ListUtils::create<String>("advanced"));
      defaults.setMinFloat("ICPL0", 0.0);
      defaults.setValue("ICPL4", 109.046571, "ICPL:2H(4)  |  H(-1) 2H(4) C(6) N O  |  unimod #687", ListUtils::create<String>("advanced"));
      defaults.setMinFloat("ICPL4", 0.0);
      defaults.setValue("ICPL6", 111.041593, "ICPL:13C(6)  |  H(3) 13C(6) N O  |  unimod #364", ListUtils::create<String>("advanced"));
      defaults.setMinFloat("ICPL6", 0.0);
      defaults.setValue("ICPL10", 115.066700, "ICPL:13C(6)2H(4)  |  H(-1) 2H(4) 13C(6) N O  |  unimod #866", ListUtils::create<String>("advanced"));
      defaults.setMinFloat("ICPL10", 0.0);
      //defaults.setValue("18O", 2.004246, "Label:18O(1)  |  O(-1) 18O  |  unimod #258", ListUtils::create<String>("advanced"));
      //defaults.setMinFloat("18O", 0.0);
    }

    return defaults;
  }

  /**
   * @brief process parameters of 'input/output' section
   */
  void getParameters_in_out_()
  {
    in_ = getStringOption_("in");
    out_ = getStringOption_("out");
    out_conflicts_ = getStringOption_("out_conflicts");
  }

  /**
   * @brief process parameters of 'algorithm' section
   */
  void getParameters_algorithm_()
  {
    // get selected labels
    labels_ = getParam_().getValue("algorithm:labels");
    samples_labels_ = splitLabelString_();

    // get selected charge range
    String charge_string = getParam_().getValue("algorithm:charge");
    double charge_min_temp, charge_max_temp;
    parseRange_(charge_string, charge_min_temp, charge_max_temp);
    charge_min_ = charge_min_temp;
    charge_max_ = charge_max_temp;
    if (charge_min_ > charge_max_)
    {
      swap(charge_min_, charge_max_);
    }

    // get isotopes per peptide range
    String isotopes_per_peptide_string = getParam_().getValue("algorithm:isotopes_per_peptide");
    double isotopes_per_peptide_min_temp, isotopes_per_peptide_max_temp;
    parseRange_(isotopes_per_peptide_string, isotopes_per_peptide_min_temp, isotopes_per_peptide_max_temp);
    isotopes_per_peptide_min_ = isotopes_per_peptide_min_temp;
    isotopes_per_peptide_max_ = isotopes_per_peptide_max_temp;
    if (isotopes_per_peptide_min_ > isotopes_per_peptide_max_)
    {
      swap(isotopes_per_peptide_min_, isotopes_per_peptide_max_);
    }

    missed_cleavages_ = getParam_().getValue("algorithm:missed_cleavages");
    
    mass_tolerance_ = getParam_().getValue("algorithm:mass_tolerance");
  }

  /**
   * @brief process parameters of 'labels' section
   */
  void getParameters_labels_()
  {
    // create map of pairs (label as string, mass shift as double)
    label_mass_shift_.insert(make_pair("Arg6", getParam_().getValue("labels:Arg6")));
    label_mass_shift_.insert(make_pair("Arg10", getParam_().getValue("labels:Arg10")));
    label_mass_shift_.insert(make_pair("Lys4", getParam_().getValue("labels:Lys4")));
    label_mass_shift_.insert(make_pair("Lys6", getParam_().getValue("labels:Lys6")));
    label_mass_shift_.insert(make_pair("Lys8", getParam_().getValue("labels:Lys8")));
    label_mass_shift_.insert(make_pair("Dimethyl0", getParam_().getValue("labels:Dimethyl0")));
    label_mass_shift_.insert(make_pair("Dimethyl4", getParam_().getValue("labels:Dimethyl4")));
    label_mass_shift_.insert(make_pair("Dimethyl6", getParam_().getValue("labels:Dimethyl6")));
    label_mass_shift_.insert(make_pair("Dimethyl8", getParam_().getValue("labels:Dimethyl8")));
    label_mass_shift_.insert(make_pair("ICPL0", getParam_().getValue("labels:ICPL0")));
    label_mass_shift_.insert(make_pair("ICPL4", getParam_().getValue("labels:ICPL4")));
    label_mass_shift_.insert(make_pair("ICPL6", getParam_().getValue("labels:ICPL6")));
    label_mass_shift_.insert(make_pair("ICPL10", getParam_().getValue("labels:ICPL10")));
  }

  /**
   * @brief split labels string
   *
   * @return list of samples containing lists of corresponding labels
   */
  std::vector<std::vector<String> > splitLabelString_()
  {
    std::vector<std::vector<String> > samples_labels;
    std::vector<String> temp_samples;
    
    String labels(labels_);
    boost::replace_all(labels, "[]", "no_label");
    boost::replace_all(labels, "()", "no_label");
    boost::replace_all(labels, "{}", "no_label");
    boost::split(temp_samples, labels, boost::is_any_of("[](){}")); // any bracket allowed to separate samples
    
    for (unsigned i = 0; i < temp_samples.size(); ++i)
    {
      if (!temp_samples[i].empty())
      {
        if (temp_samples[i]=="no_label")
        {
          vector<String> temp_labels;
          temp_labels.push_back("no_label");
          samples_labels.push_back(temp_labels);
        }
        else
        {
          vector<String> temp_labels;
          boost::split(temp_labels, temp_samples[i], boost::is_any_of(",;: ")); // various separators allowed to separate labels
          samples_labels.push_back(temp_labels);
        }
      }
    }
    
    if (samples_labels.empty())
    {
      vector<String> temp_labels;
      temp_labels.push_back("no_label");
      samples_labels.push_back(temp_labels);
    }

    return samples_labels;
  }

  /**
   * @brief comparator of peak patterns
   *
   * @param pattern1    first peak pattern
   * @param pattern2    second peak pattern
   *
   * @return true if pattern1 should be searched before pattern2
   */
  static bool less_pattern(const MultiplexIsotopicPeakPattern& pattern1, const MultiplexIsotopicPeakPattern& pattern2)
  {
    if (pattern1.getMassShiftCount() == pattern2.getMassShiftCount())
    {
      if (pattern1.getCharge() == pattern2.getCharge())
      {
        // The first mass shift is by definition always zero.
        if ((pattern1.getMassShiftCount() > 1) && (pattern2.getMassShiftCount() > 1))
        {
          // 4Da before 8Da etc. (larger miss cleavages last)
          return pattern1.getMassShiftAt(1) < pattern2.getMassShiftAt(1);
        }
        else
        {
          // Should never happen.
          return true;
        }
      }
      else
      {
        // 5+ before 4+ before 3+ etc.
        return pattern1.getCharge() > pattern2.getCharge();
      }
    }
    else
    {
      // triplets before doublets before singlets
      return pattern1.getMassShiftCount() > pattern2.getMassShiftCount();
    }
  }

  /**
   * @brief generate list of m/z shifts
   *
   * @param charge_min    minimum charge
   * @param charge_max    maximum charge
   * @param peaks_per_peptide_max    maximum number of isotopes in peptide
   * @param mass_pattern_list    mass shifts due to labelling
   *
   * @return list of m/z shifts
   */
  std::vector<MultiplexIsotopicPeakPattern> generatePeakPatterns_(int charge_min, int charge_max, int peaks_per_peptide_max, std::vector<MultiplexDeltaMasses> mass_pattern_list)
  {
    std::vector<MultiplexIsotopicPeakPattern> list;

    // iterate over all charge states
    for (int c = charge_max; c >= charge_min; --c)
    {
      // iterate over all mass shifts
      for (unsigned i = 0; i < mass_pattern_list.size(); ++i)
      {
        MultiplexIsotopicPeakPattern pattern(c, peaks_per_peptide_max, mass_pattern_list[i], i);
        list.push_back(pattern);
      }
    }
    
    sort(list.begin(),list.end(),less_pattern);
    
    return list;
  }
  
  /**
   * @brief returns the relative delta mass between the first feature
   * and the feature with the map index idx
   *
   * @param feature_handles    feature handles of a consensus feature
   * @param idx    map index of interest
   */
  double DeltaMassFromMapIndex(const ConsensusFeature::HandleSetType& feature_handles, unsigned idx)
  {
    double first_mass = feature_handles.begin()->getMZ() * feature_handles.begin()->getCharge();
    
    for (ConsensusFeature::HandleSetType::const_iterator it_feat = feature_handles.begin(); it_feat != feature_handles.end(); ++it_feat)
    {
      if (it_feat->getMapIndex() == idx)
      {
        return it_feat->getMZ() * it_feat->getCharge() - first_mass;
      }
    }
    
    // return NaN if no matching index was found
    return numeric_limits<double>::quiet_NaN();
  }
  
  /**
   * @brief check whether the theoretical delta mass pattern
   * contains the label set of the detected pattern
   *
   * @param pattern    theoretical pattern
   * @param label_set    label set of the detected pettern
   * @param index_label_set    index within the pattern at which the label sets were matched
   * 
   * @return mass shift in the theoretical pattern where both label sets match
   */
  double matchLabelSet(const std::vector<MultiplexDeltaMasses::DeltaMass>& pattern, const MultiplexDeltaMasses::LabelSet& label_set, int& index_label_set)
  {
    for (std::vector<MultiplexDeltaMasses::DeltaMass>::const_iterator it_mass_shift = pattern.begin(); it_mass_shift != pattern.end(); ++it_mass_shift)
    {
      if (it_mass_shift->label_set == label_set)
      {
        index_label_set = it_mass_shift - pattern.begin();
        return it_mass_shift->delta_mass;
      }
    }
    
    // return NaN if no matching label set was found
    return numeric_limits<double>::quiet_NaN();
  }

  /**
   * @brief check wether all delta masses in the detected patter
   * match up with a delta mass in the theoretical pattern
   *
   * @param consensus    detected pattern
   * @param pattern    theoretical pattern
   * @param delta_mass_at_label_set    delta mass in the theoretical pattern at which the matching label set was found
   * @param delta_mass_matched    Was this delta mass in the theoretical pattern matched?
   * 
   * @return All delta masses matching?
   */
  bool matchDeltaMasses(const ConsensusMap::ConstIterator consensus, const std::vector<MultiplexDeltaMasses::DeltaMass>& pattern, double theoretical_delta_mass_at_label_set, std::vector<bool>& delta_mass_matched)
  {
    double first_mass = consensus->getFeatures().begin()->getMZ() * consensus->getFeatures().begin()->getCharge();
    double detected_delta_mass_at_label_set = DeltaMassFromMapIndex(consensus->getFeatures(), consensus->getPeptideIdentifications()[0].getMetaValue("map index"));
    if (isnan(detected_delta_mass_at_label_set))
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "No delta mass with this map index could be found.", "");
    }
    
    // loop over features in consensus    
    for (ConsensusFeature::HandleSetType::const_iterator it_feat = consensus->getFeatures().begin(); it_feat != consensus->getFeatures().end(); ++it_feat)
    {
      // delta mass in the detected pattern relative to the feature with the matched label set
      double mass_shift_detected = (it_feat->getMZ() * it_feat->getCharge() - first_mass) - detected_delta_mass_at_label_set;
      bool matched = false;
      
      // loop over delta masses in theoretical pattern
      for (std::vector<MultiplexDeltaMasses::DeltaMass>::const_iterator it_mass_shift = pattern.begin(); it_mass_shift != pattern.end(); ++it_mass_shift)
      {
        // delta mass in the theoretical pattern relative to the feature with the matched label set
        double mass_shift_theoretical = it_mass_shift->delta_mass - theoretical_delta_mass_at_label_set;
        
        if (abs(mass_shift_detected - mass_shift_theoretical) < mass_tolerance_)
        {
          delta_mass_matched[it_mass_shift - pattern.begin()] = true;
          matched = true;
          break;
        }
      }
      
      if (!matched)
      {
        return false;
      }      
    }
    
    return true;
  }
  
  /**
   * @brief find a theoretical delta mass pattern that matches the detected pattern
   *
   * @param consensus    detected pattern
   * @param label set    label set extracted from the detected pattern
   * @param theoretical_patterns    list of theoretical delta mass patterns
   * @param delta_mass_matched    Was this delta mass in the theoretical pattern matched?
   * @param index_label_set    index within the pattern at which the label sets were matched
   * 
   * @return index of matching pattern
   */
  int findMatchingPattern(const ConsensusMap::ConstIterator consensus, const MultiplexDeltaMasses::LabelSet& label_set, const std::vector<MultiplexDeltaMasses>& theoretical_patterns, std::vector<bool>& delta_mass_matched, int& index_label_set)
  {
    // loop over theoretical patterns
    for (std::vector<MultiplexDeltaMasses>::const_iterator it_pattern = theoretical_patterns.begin(); it_pattern != theoretical_patterns.end(); ++it_pattern)
    {
      std::vector<MultiplexDeltaMasses::DeltaMass> pattern = it_pattern->getDeltaMasses();
      
      double shift = matchLabelSet(pattern, label_set, index_label_set);
      if (!isnan(shift))
      {        
        // reset boolean vector
        unsigned i = delta_mass_matched.size();
        delta_mass_matched.assign(i,false);
        
        bool match = matchDeltaMasses(consensus, pattern, shift, delta_mass_matched);
        if (match)
        {
          return (it_pattern - theoretical_patterns.begin());
        }
      }
    }
    
    return -1;
  }
  
  /**
   * @brief find the m/z for the complete consensus
   *
   * @param mz    m/z of the incomplete consensus
   * @param charge    charge of the incomplete consensus
   * @param pattern    matching theoretical delta mass pattern
   * @param delta_mass_matched    Was this delta mass in the theoretical pattern matched?
   * 
   * @return m/z for the complete consensus
   */
  double findNewMZ(double mz, int charge, const std::vector<MultiplexDeltaMasses::DeltaMass>& pattern, const std::vector<bool>& delta_mass_matched)
  {
    // loop over delta masses in theoretical pattern
    std::vector<MultiplexDeltaMasses::DeltaMass>::const_iterator it_mass_shift;
    std::vector<bool>::const_iterator it_delta_mass_matched;
    for (it_mass_shift = pattern.begin(), it_delta_mass_matched = delta_mass_matched.begin();
         it_mass_shift != pattern.end(), it_delta_mass_matched != delta_mass_matched.end();
         ++it_mass_shift, ++it_delta_mass_matched)
    {
      // find the first match
      if (*it_delta_mass_matched)
      {
        return (mz * charge - it_mass_shift->delta_mass)/charge;
      }
    }
    
    // Should never happen.
    return mz;
  }
  
  /**
   * @brief complete consensus
   *
   * @param consensus    (possibly) incomplete consensus
   * @param pattern    matching theoretical delta mass pattern
   * @param delta_mass_matched    Was this delta mass in the theoretical pattern matched?
   * @param index_label_set    index within the pattern at which the label sets were matched
   * 
   * @return completed consensus
   */
  ConsensusFeature completeConsensus(const ConsensusFeature& consensus, const std::vector<MultiplexDeltaMasses::DeltaMass>& pattern, const std::vector<bool>& delta_mass_matched, int index_label_set)
  {
    // Nothing to do. Detected consensus is already complete.
    if (consensus.size() == pattern.size())
    {
      return ConsensusFeature(consensus);
    }
    
    if (pattern.size() != delta_mass_matched.size())
    {
       throw Exception::InvalidSize(__FILE__, __LINE__, __PRETTY_FUNCTION__, delta_mass_matched.size());
    }
    
    // new complete consensus feature
    ConsensusFeature consensus_complete;
    
    int charge = consensus.getCharge();
    double RT = consensus.getRT();
    double mz = consensus.getMZ();
    
    // find m/z of the new complete consensus
    double mz_complete = findNewMZ(mz, charge, pattern, delta_mass_matched);
    
    consensus_complete.setMZ(mz_complete);
    consensus_complete.setRT(consensus.getRT());
    consensus_complete.setCharge(consensus.getCharge());
    consensus_complete.setIntensity(consensus.getIntensity());    // Alternatively, reduce intensity due to new zero-intensity dummy features.
    consensus_complete.setQuality(consensus.getQuality());
    consensus_complete.setPeptideIdentifications(consensus.getPeptideIdentifications());
    consensus_complete.getPeptideIdentifications()[0].getHits()[0].setMetaValue("map index", index_label_set);
    
    // loop over delta masses in theoretical pattern
    std::vector<MultiplexDeltaMasses::DeltaMass>::const_iterator it_mass_shift;
    std::vector<bool>::const_iterator it_delta_mass_matched;
    ConsensusFeature::HandleSetType::const_iterator it_feature;
    for (it_mass_shift = pattern.begin(), it_delta_mass_matched = delta_mass_matched.begin(), it_feature = consensus.getFeatures().begin();
         it_mass_shift != pattern.end(), it_delta_mass_matched != delta_mass_matched.end();
         ++it_mass_shift, ++it_delta_mass_matched)
    {
      
      //std::cout << "    index = " << (it_mass_shift - pattern.begin()) << "    shift = " << it_mass_shift->delta_mass;
      if (*it_delta_mass_matched)
      {
        // copy feature from incomplete consensus
        FeatureHandle feature_handle(*it_feature);
        feature_handle.setMapIndex(it_mass_shift - pattern.begin());
        consensus_complete.insert(feature_handle);
        
        if (it_feature != consensus.getFeatures().end())
        {
          ++it_feature;
        }
      }
      else
      {
        // construct dummy feature
        FeatureHandle feature_handle;
        feature_handle.setMZ(mz_complete + it_mass_shift->delta_mass/charge);
        feature_handle.setRT(RT);
        feature_handle.setIntensity(0.0);
        feature_handle.setCharge(charge);
        feature_handle.setMapIndex(it_mass_shift - pattern.begin());
        consensus_complete.insert(feature_handle);
      }
      
    }
    
    return consensus_complete;
  }
  
  /**
   * @brief construct the new consensus map
   * (1) remove quant/ID conflicts
   * (2) fill in dummy features in order to complete multiplets
   *
   * @param map_in    input consensus map
   * @param map_out    consensus map without conflicts and complete multiplets
   * @param map_conflicts    consensus map with conflicts
   * @param generator    generator for the list of theoretical patterns
   */
  void constructNewConsensusMap_(const ConsensusMap& map_in, ConsensusMap& map_out, ConsensusMap& map_conflicts, MultiplexDeltaMassesGenerator generator)
  {
    unsigned found_pattern_count = 0;
    std::vector<MultiplexDeltaMasses> theoretical_masses = generator.getDeltaMassesList();
    unsigned multiplicity = theoretical_masses[0].getDeltaMasses().size();
    
    for (ConsensusMap::ConstIterator cit = map_in.begin(); cit != map_in.end(); ++cit)
    {
      // extract the label set from the attached peptide sequence (There should be only one, since IDConflictResolver was run first.)
      AASequence sequence = cit->getPeptideIdentifications()[0].getHits()[0].getSequence();      
      MultiplexDeltaMasses::LabelSet label_set = generator.extractLabelSet(sequence);
      std::vector<bool> delta_mass_matched(multiplicity, false);
      int index_label_set = -1;
      
      int index = findMatchingPattern(cit, label_set, theoretical_masses, delta_mass_matched, index_label_set);
      
      /*std::cout << "consensus = " << (cit - map_in.begin());
      std::cout << "    RT = " << cit->getRT();
      std::cout << "    sequence = " << sequence;*/
      
      if (index >= 0)
      {
        //std::cout << "  (Ok)\n\n";
        
        ++found_pattern_count;
        
        ConsensusFeature consensus = completeConsensus(*cit, theoretical_masses[index].getDeltaMasses(), delta_mass_matched, index_label_set);
        map_out.push_back(consensus);
      }
      else
      {
        //std::cout << "  (Conflict)\n\n";
        
        ConsensusFeature consensus(*cit);
        map_conflicts.push_back(consensus);
      }
    }
    
    // update map sizes
    for (unsigned map_index = 0; map_index < multiplicity; ++map_index)
    {
      map_out.getFileDescriptions()[map_index].size = map_out.size();
    }
    
    map_out.applyMemberFunction(&UniqueIdInterface::setUniqueId);
    map_conflicts.applyMemberFunction(&UniqueIdInterface::setUniqueId);
    
    /*std::cout << "\n";
    std::cout << "number of consensuses                   = " << map_in.size() << "\n";
    std::cout << "number of consensuses without conflicts = " << found_pattern_count << "\n";
    std::cout << "\n";*/
  }
  
  ExitCodes main_(int, const char**)
  {
    /**
     * handle parameters
     */
    getParameters_in_out_();
    getParameters_labels_();
    getParameters_algorithm_();

    /**
     * load consensus map
     */
    ConsensusXMLFile file;
    ConsensusMap map_in;
    file.load(in_, map_in);

    /**
     * generate patterns
     */
    MultiplexDeltaMassesGenerator generator = MultiplexDeltaMassesGenerator(labels_, missed_cleavages_, label_mass_shift_);
    generator.printSamplesLabelsList();
    generator.printDeltaMassesList();
        
    /**
     * construct the new consensus map
     */
    ConsensusMap map_out = map_in;
    ConsensusMap map_conflicts = map_in;
    map_out.resize(0);
    map_conflicts.resize(0);
    constructNewConsensusMap_(map_in, map_out, map_conflicts, generator);
         
    /**
     * store consensus maps
     */
    ConsensusXMLFile file_out;
    ConsensusXMLFile file_out_conflicts;
    file_out.store(out_, map_out);
    file_out_conflicts.store(out_conflicts_, map_conflicts);
   
    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPMultiplexResolver tool;
  return tool.main(argc, argv);
}

//@endcond
