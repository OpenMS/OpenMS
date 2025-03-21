// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FEATUREFINDER/MultiplexDeltaMasses.h>
#include <OpenMS/FEATUREFINDER/MultiplexDeltaMassesGenerator.h>
#include <OpenMS/FEATUREFINDER/MultiplexIsotopicPeakPattern.h>

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/replace.hpp>

using namespace std;
using namespace OpenMS;

//#define DEBUG

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_MultiplexResolver MultiplexResolver

@brief Completes peptide multiplets and resolves conflicts within them.

<CENTER>
  <table>
    <tr>
      <th ALIGN = "center"> pot. predecessor tools </td>
      <td VALIGN="middle" ROWSPAN=3> &rarr; MultiplexResolver &rarr;</td>
      <th ALIGN = "center"> pot. successor tools </td>
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
using the IDMapper tool (*). The MultiplexResolver tool is consolidating these results in two steps. 
- Any multiplets with conflicting quantitative and sequence information are filtered out. As example, let us consider a triple SILAC analyis. Let us assume a sequence
"LDNLVAIFDINR(Label:13C(6)15N(4))" with a single Arg10 label is mapped to the light feature in a SILAC triplet. Either peptide feature detection or sequence information
must be incorrect und the triplet is removed.
- In a second step, any incomplete peptide feature groups are completed with dummy features of zero intensity. As example, let us stay with the triple SILAC analysis.
But let us now assume the sequence "LDNLVAIFDINR(Label:13C(6)15N(4))" is mapped to the heavy partner of a peptide feature pair. This is no conflict. Medium and heavy
peptides have been correctly detected. The MultiplexResolver adds a dummy peptide feature of zero intensity at the light position and thereby completes the triplet.

(*) Note that the MultiplexResolver tool takes only a single (the first) peptide sequence annotation into account. By running IDConflictResolver first, it is assured that
each multiplet has only one peptide sequence annotation, the best one. Multiplets without sequence annotation are passed to the optional out_conflicts output.

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
  String in_blacklist_;
  String out_;
  String out_conflicts_;

  // section "algorithm"
  String labels_;
  unsigned missed_cleavages_;
  double mass_tolerance_;
  double mz_tolerance_;
  double rt_tolerance_;

  // section "labels"
  map<String, double> label_mass_shift_;
  
  // blacklist
  MSExperiment exp_blacklist_;
  
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Peptide multiplets with assigned sequence information");
    setValidFormats_("in", ListUtils::create<String>("consensusXML"));
    registerInputFile_("in_blacklist", "<file>", "", "Optional input containing spectral peaks blacklisted during feature detection. Needed for generation of dummy features.", false);
    setValidFormats_("in_blacklist", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "Complete peptide multiplets.");
    setValidFormats_("out", ListUtils::create<String>("consensusXML"));
    registerOutputFile_("out_conflicts", "<file>", "", "Optional output containing peptide multiplets without ID annotation or with conflicting quant/ID information.", false);
    setValidFormats_("out_conflicts", ListUtils::create<String>("consensusXML"));

    registerSubsection_("algorithm", "Parameters for the algorithm.");
    registerSubsection_("labels", "Isotopic labels that can be specified in section \'algorithm:labels\'.");
    
  }
  
  // create parameters for sections (set default values and restrictions)
  Param getSubsectionDefaults_(const String& section) const override
  {
    Param defaults;

    if (section == "algorithm")
    {
      defaults.setValue("labels", "[][Lys8,Arg10]", "Labels used for labelling the samples. [...] specifies the labels for a single sample. For example\n\n[][Lys8,Arg10]        ... SILAC\n[][Lys4,Arg6][Lys8,Arg10]        ... triple-SILAC\n[Dimethyl0][Dimethyl6]        ... Dimethyl\n[Dimethyl0][Dimethyl4][Dimethyl8]        ... triple Dimethyl\n[ICPL0][ICPL4][ICPL6][ICPL10]        ... ICPL");
      defaults.setValue("missed_cleavages", 0, "Maximum number of missed cleavages due to incomplete digestion. (Only relevant if enzymatic cutting site coincides with labelling site. For example, Arg/Lys in the case of trypsin digestion and SILAC labelling.)");
      defaults.setMinInt("missed_cleavages", 0);
      defaults.setValue("mass_tolerance", 0.1, "Mass tolerance in Da for matching the mass shifts in the detected peptide multiplet to the theoretical mass shift pattern.", {"advanced"});
      defaults.setValue("mz_tolerance", 10, "m/z tolerance in ppm for checking if dummy feature vicinity was blacklisted.", {"advanced"});
      defaults.setValue("rt_tolerance", 5, "Retention time tolerance in seconds for checking if dummy feature vicinity was blacklisted.", {"advanced"});
    }

    if (section == "labels")
    {     
      MultiplexDeltaMassesGenerator generator;
      Param p = generator.getParameters();
      
      for (Param::ParamIterator it = p.begin(); it != p.end(); ++it)
      {
        defaults.setValue(it->name, it->value, it->description, {"advanced"});
        defaults.setMinFloat(it->name, 0.0);
      }
    }

    return defaults;
  }

  /**
   * @brief process parameters of 'input/output' section
   */
  void getParameters_in_out_()
  {
    in_ = getStringOption_("in");
    in_blacklist_ = getStringOption_("in_blacklist");
    out_ = getStringOption_("out");
    out_conflicts_ = getStringOption_("out_conflicts");
  }

  /**
   * @brief process parameters of 'algorithm' section
   */
  void getParameters_algorithm_()
  {
    labels_ = getParam_().getValue("algorithm:labels").toString();
    missed_cleavages_ = getParam_().getValue("algorithm:missed_cleavages");
    mass_tolerance_ = getParam_().getValue("algorithm:mass_tolerance");
    mz_tolerance_ = getParam_().getValue("algorithm:mz_tolerance");
    rt_tolerance_ = getParam_().getValue("algorithm:rt_tolerance");
  }

  /**
   * @brief process parameters of 'labels' section
   */
  void getParameters_labels_()
  {
    Param p = getParam_();
    
    // create map of pairs (label as string, mass shift as double)
    for (Param::ParamIterator it = p.begin(); it != p.end(); ++it)
    {
      label_mass_shift_.insert(make_pair(it->name, it->value));
    }
  }

  /**
   * @brief returns the relative delta mass between the first feature
   * and the feature with the map index idx
   *
   * @param feature_handles    feature handles of a consensus feature
   * @param idx    map index of interest
   */
  double deltaMassFromMapIndex_(const ConsensusFeature::HandleSetType& feature_handles, unsigned idx)
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
  double matchLabelSet_(const std::vector<MultiplexDeltaMasses::DeltaMass>& pattern, const MultiplexDeltaMasses::LabelSet& label_set, int& index_label_set)
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
   * @brief check wether all delta masses in the detected pattern
   * match up with a delta mass in the theoretical pattern
   *
   * @param consensus    detected pattern
   * @param pattern    theoretical pattern
   * @param delta_mass_at_label_set    delta mass in the theoretical pattern at which the matching label set was found
   * @param delta_mass_matched    Was this delta mass in the theoretical pattern matched?
   * 
   * @return All delta masses matching?
   */
  bool matchDeltaMasses_(const ConsensusMap::ConstIterator consensus, const std::vector<MultiplexDeltaMasses::DeltaMass>& pattern, double theoretical_delta_mass_at_label_set, std::vector<bool>& delta_mass_matched)
  {
    double first_mass = consensus->getFeatures().begin()->getMZ() * consensus->getFeatures().begin()->getCharge();
    if (!consensus->getPeptideIdentifications()[0].metaValueExists("map_index"))
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The meta value 'map_index' is missing in the input data. In the IDMapper tool, please set the advanced parameter consensus:annotate_ids_with_subelements = true.");
    }
    double detected_delta_mass_at_label_set = deltaMassFromMapIndex_(consensus->getFeatures(), consensus->getPeptideIdentifications()[0].getMetaValue("map_index"));
    if (std::isnan(detected_delta_mass_at_label_set))
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No delta mass with this map_index could be found.", "");
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
  int findMatchingPattern_(const ConsensusMap::ConstIterator consensus, const MultiplexDeltaMasses::LabelSet& label_set, const std::vector<MultiplexDeltaMasses>& theoretical_patterns, std::vector<bool>& delta_mass_matched, int& index_label_set)
  {
    // loop over theoretical patterns
    for (std::vector<MultiplexDeltaMasses>::const_iterator it_pattern = theoretical_patterns.begin(); it_pattern != theoretical_patterns.end(); ++it_pattern)
    {
      std::vector<MultiplexDeltaMasses::DeltaMass> pattern = it_pattern->getDeltaMasses();
      
      double shift = matchLabelSet_(pattern, label_set, index_label_set);
      if (!std::isnan(shift))
      {        
        // reset boolean vector to false
        delta_mass_matched.assign(delta_mass_matched.size(), false);
        
        bool match = matchDeltaMasses_(consensus, pattern, shift, delta_mass_matched);
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
  double findNewMZ_(double mz, int charge, const std::vector<MultiplexDeltaMasses::DeltaMass>& pattern, const std::vector<bool>& delta_mass_matched)
  {
    // loop over delta masses in theoretical pattern
    std::vector<MultiplexDeltaMasses::DeltaMass>::const_iterator it_mass_shift;
    std::vector<bool>::const_iterator it_delta_mass_matched;
    for (it_mass_shift = pattern.begin(), it_delta_mass_matched = delta_mass_matched.begin();
         it_mass_shift != pattern.end() && it_delta_mass_matched != delta_mass_matched.end();
         ++it_mass_shift, ++it_delta_mass_matched)
    {
      // find the first match
      if (*it_delta_mass_matched)
      {
        return (mz * charge - it_mass_shift->delta_mass) / charge;
      }
    }
    
    // Should never happen.
    return mz;
  }
  
  /**
   * @brief check if this position is blacklisted
   * 
   * @param RT
   * @param mz
   * @param charge
   */
  bool isBlacklisted(double rt, double mz, size_t charge)
  {
    double mz_tolerance = mz_tolerance_ * mz / 1000000;    // m/z tolerance in Da
    
    MSExperiment::ConstIterator it_rt_begin = exp_blacklist_.RTBegin(rt - rt_tolerance_);
    MSExperiment::ConstIterator it_rt_end = exp_blacklist_.RTEnd(rt + rt_tolerance_);
    
    // loop over range of relevant spectra
    for (MSExperiment::ConstIterator it_rt = it_rt_begin; it_rt < it_rt_end; ++it_rt)
    {
      // Loop over first three isotopes in dummy feature (and check if one of them is blacklisted).
      for (size_t isotope = 0; isotope < 3; ++isotope)
      {
        double mz_isotope = mz + isotope * Constants::C13C12_MASSDIFF_U / charge;
        
        MSSpectrum::ConstIterator it_mz = it_rt->MZBegin(mz_isotope);
        
        if ((std::abs(it_mz->getMZ() - mz_isotope)) < mz_tolerance)
        {
          // There is a blacklisted peak close-by.
          return true;
        }
      }
    }
    
    // None of the first three isotopes has a blacklisted peak near-by.
    return false;
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
  ConsensusFeature completeConsensus_(const ConsensusFeature& consensus, const std::vector<MultiplexDeltaMasses::DeltaMass>& pattern, const std::vector<bool>& delta_mass_matched, int index_label_set)
  {
    // Nothing to do. Detected consensus is already complete.
    if (consensus.size() == pattern.size())
    {
      return ConsensusFeature(consensus);
    }
    
    if (pattern.size() != delta_mass_matched.size())
    {
       throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, delta_mass_matched.size());
    }
    
    // new complete consensus feature
    ConsensusFeature consensus_complete;
    
    int charge = consensus.getCharge();
    double RT = consensus.getRT();
    double mz = consensus.getMZ();
    
    // find m/z of the new complete consensus
    double mz_complete = findNewMZ_(mz, charge, pattern, delta_mass_matched);
    
    consensus_complete.setMZ(mz_complete);
    consensus_complete.setRT(consensus.getRT());
    consensus_complete.setCharge(consensus.getCharge());
    consensus_complete.setIntensity(consensus.getIntensity());    // Alternatively, reduce intensity due to new zero-intensity dummy features.
    consensus_complete.setQuality(consensus.getQuality());
    consensus_complete.setPeptideIdentifications(consensus.getPeptideIdentifications());
    consensus_complete.getPeptideIdentifications()[0].getHits()[0].setMetaValue("map_index", index_label_set);
    
    // loop over delta masses in theoretical pattern
    std::vector<MultiplexDeltaMasses::DeltaMass>::const_iterator it_mass_shift;
    std::vector<bool>::const_iterator it_delta_mass_matched;
    ConsensusFeature::HandleSetType::const_iterator it_feature;
    for (it_mass_shift = pattern.begin(), it_delta_mass_matched = delta_mass_matched.begin(), it_feature = consensus.getFeatures().begin();
         it_mass_shift != pattern.end() && it_delta_mass_matched != delta_mass_matched.end();
         ++it_mass_shift, ++it_delta_mass_matched)
    {
      
      //OPENMS_LOG_DEBUG << "    index = " << (it_mass_shift - pattern.begin()) << "    shift = " << it_mass_shift->delta_mass;
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
        feature_handle.setMZ(mz_complete + it_mass_shift->delta_mass / charge);
        feature_handle.setRT(RT);
        if (isBlacklisted(RT, mz_complete + it_mass_shift->delta_mass / charge, charge))
        {
          // Some peaks close-by were blacklisted during feature detection i.e. another peptide feature overlaps with the dummy feature.
          // Consequently, we better report NaN i.e. not quantifiable.
          feature_handle.setIntensity(std::numeric_limits<double>::quiet_NaN());
        }
        else
        {
          // There is no blacklisted peak near-by i.e. there is no peptide feature in the vicinity.
          // Consequently, we can confidently report zero i.e. the peptide is absent.
          feature_handle.setIntensity(0.0);
        }
        feature_handle.setCharge(charge);
        feature_handle.setMapIndex(it_mass_shift - pattern.begin());
        consensus_complete.insert(feature_handle);
        
        // debug output
        //std::cout << "dummy feature @ RT = " << RT << "   m/z = " << (mz_complete + it_mass_shift->delta_mass / charge) << "   blacklisted = " << isBlacklisted(RT, mz_complete + it_mass_shift->delta_mass / charge, charge) << "\n";
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
    // unsigned found_pattern_count = 0;
    std::vector<MultiplexDeltaMasses> theoretical_masses = generator.getDeltaMassesList();
    size_t multiplicity = theoretical_masses[0].getDeltaMasses().size();
    
    for (ConsensusMap::ConstIterator cit = map_in.begin(); cit != map_in.end(); ++cit)
    {
      //OPENMS_LOG_DEBUG << "consensus = " << (cit - map_in.begin());
      //OPENMS_LOG_DEBUG << "       RT = " << cit->getRT();
      //OPENMS_LOG_DEBUG << "       mz = " << cit->getMZ();
      
      // Consensus features without sequence annotations are written unchanged to the conflict output.
      if (cit->getPeptideIdentifications().empty())
      {
        //OPENMS_LOG_DEBUG << "  (no ID)\n\n";
        
        ConsensusFeature consensus(*cit);
        map_conflicts.push_back(consensus);
        
        continue;
      }
      
      // extract the label set from the attached peptide sequence (There should be only one, since IDConflictResolver was run first.)
      AASequence sequence = cit->getPeptideIdentifications()[0].getHits()[0].getSequence();      
      MultiplexDeltaMasses::LabelSet label_set = generator.extractLabelSet(sequence);
      std::vector<bool> delta_mass_matched(multiplicity, false);
      int index_label_set = -1;
      
      int index = findMatchingPattern_(cit, label_set, theoretical_masses, delta_mass_matched, index_label_set);
            
      if (index >= 0)
      {
        //OPENMS_LOG_DEBUG << "  (Ok)\n\n";
        // ++found_pattern_count;
        
        ConsensusFeature consensus = completeConsensus_(*cit, theoretical_masses[index].getDeltaMasses(), delta_mass_matched, index_label_set);
        map_out.push_back(consensus);
      }
      else
      {
        //OPENMS_LOG_DEBUG << "  (Conflict)\n\n";
        
        ConsensusFeature consensus(*cit);
        map_conflicts.push_back(consensus);
      }
    }
    
    // update map sizes
    for (unsigned map_index = 0; map_index < multiplicity; ++map_index)
    {
      map_out.getColumnHeaders()[map_index].size = map_out.size();
    }
    
    map_out.applyMemberFunction(&UniqueIdInterface::setUniqueId);
    map_conflicts.applyMemberFunction(&UniqueIdInterface::setUniqueId);
    
    /*OPENMS_LOG_DEBUG << "\n";
    OPENMS_LOG_DEBUG << "number of consensuses                   = " << map_in.size() << "\n";
    OPENMS_LOG_DEBUG << "number of consensuses without conflicts = " << found_pattern_count << "\n";
    OPENMS_LOG_DEBUG << "\n";*/
  }
  
public:
  TOPPMultiplexResolver() :
    TOPPBase("MultiplexResolver", "Completes peptide multiplets and resolves conflicts within them."),
    missed_cleavages_(0), mass_tolerance_(0.1)
  {
  }

  ExitCodes main_(int, const char**) override
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
    ConsensusMap map_in;
    FileHandler().loadConsensusFeatures(in_, map_in, {FileTypes::CONSENSUSXML});

    /**
     * load (optional) blacklist
     */
    if (!(in_blacklist_.empty()))
    {
      FileHandler().loadExperiment(in_blacklist_, exp_blacklist_, {FileTypes::MZML});
    }

    /**
     * generate patterns
     */
    MultiplexDeltaMassesGenerator generator = MultiplexDeltaMassesGenerator(labels_, missed_cleavages_, label_mass_shift_);
    #ifdef DEBUG
    generator.printSamplesLabelsList(std::cout);
    generator.printDeltaMassesList(std::cout);
    #endif
        
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
    FileHandler().storeConsensusFeatures(out_, map_out, {FileTypes::CONSENSUSXML});
    if (!out_conflicts_.empty())
    {
      FileHandler().storeConsensusFeatures(out_conflicts_, map_conflicts, {FileTypes::CONSENSUSXML});
    }
   
    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPMultiplexResolver tool;
  return tool.main(argc, argv);
}

///@endcond
