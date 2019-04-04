// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche, Fabian Kriegel, Frederic Lehnert $
// --------------------------------------------------------------------------

#include <OpenMS/SIMULATION/LABELING/SILACLabeler.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>

using std::vector;

namespace OpenMS
{
  const int SILACLabeler::LIGHT_FEATURE_MAPID_ = 0;
  const int SILACLabeler::MEDIUM_FEATURE_MAPID_ = 1;
  const int SILACLabeler::HEAVY_FEATURE_MAPID_ = 2;


  SILACLabeler::SILACLabeler() :
    BaseLabeler()
  {
    channel_description_ = "SILAC labeling on MS1 level with up to 3 channels and custom modifications.";

    defaults_.setValue("medium_channel:modification_lysine", "UniMod:481", "Modification of Lysine in the medium SILAC channel");
    defaults_.setValue("medium_channel:modification_arginine", "UniMod:188", "Modification of Arginine in the medium SILAC channel");
    defaults_.setSectionDescription("medium_channel", "Modifications for the medium SILAC channel.");

    defaults_.setValue("heavy_channel:modification_lysine", "UniMod:259", "Modification of Lysine in the heavy SILAC channel. If left empty, two channelSILAC is assumed.");
    defaults_.setValue("heavy_channel:modification_arginine", "UniMod:267", "Modification of Arginine in the heavy SILAC channel. If left empty, two-channel SILAC is assumed.");
    defaults_.setSectionDescription("heavy_channel", "Modifications for the heavy SILAC channel. If you want to use only 2 channels, just leave the Labels as they are and provide only 2 input files.");

    defaults_.setValue("fixed_rtshift", 0.0001, "Fixed retention time shift between labeled peptides. If set to 0.0 only the retention times computed by the RT model step are used.");
    defaults_.setMinFloat("fixed_rtshift", 0.0);

    defaultsToParam_();
  }

  void SILACLabeler::updateMembers_()
  {
    medium_channel_lysine_label_ = (String)param_.getValue("medium_channel:modification_lysine");
    medium_channel_arginine_label_ = (String)param_.getValue("medium_channel:modification_arginine");

    heavy_channel_lysine_label_ = (String)param_.getValue("heavy_channel:modification_lysine");
    heavy_channel_arginine_label_ = (String)param_.getValue("heavy_channel:modification_arginine");
  }

  bool SILACLabeler::canModificationBeApplied_(const String& modification_id, const String& aa) const
  {
    std::set<const ResidueModification*> modifications;
    ModificationsDB::getInstance()->searchModifications(modifications, modification_id, aa);
    
    if (modifications.empty())
    {
      String message = String("The modification '") + modification_id + "' could not be found in the local UniMod DB! Please check if you used the correct format (e.g. UniMod:Accession#)";
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, message);
    }

    return !modifications.empty();
  }

  SILACLabeler::~SILACLabeler()
  {

  }

  void SILACLabeler::preCheck(Param&) const
  {
    // check if the given modifications can be applied to the target
    // amino acids
    canModificationBeApplied_(medium_channel_lysine_label_, String("K"));
    canModificationBeApplied_(medium_channel_arginine_label_, String("R"));
    canModificationBeApplied_(heavy_channel_lysine_label_, String("K"));
    canModificationBeApplied_(heavy_channel_arginine_label_, String("R"));
  }

  void SILACLabeler::applyLabelToProteinHit_(SimTypes::FeatureMapSim& channel, const String& arginine_label, const String& lysine_label) const
  {
    for (std::vector<ProteinHit>::iterator protein_hit = channel.getProteinIdentifications()[0].getHits().begin();
         protein_hit != channel.getProteinIdentifications()[0].getHits().end();
         ++protein_hit)
    {
      AASequence aa = AASequence::fromString(protein_hit->getSequence());

      for (AASequence::Iterator residue = aa.begin(); residue != aa.end(); ++residue)
      {
        if (*residue == 'R')
        {
          aa.setModification(residue - aa.begin(), arginine_label);
        }
        else if (*residue == 'K')
        {
          aa.setModification(residue - aa.begin(), lysine_label);
        }
      }
      protein_hit->setSequence(aa.toString());
    }
  }

  void SILACLabeler::setUpHook(SimTypes::FeatureMapSimVector& features_to_simulate)
  {
    // check if we have the correct number of channels
    if (features_to_simulate.size() < 2 || features_to_simulate.size() > 3)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String(features_to_simulate.size()) + " channel(s) given. We currently support only 2-channel SILAC. Please provide two FASTA files!");
    }

    SimTypes::FeatureMapSim& medium_channel = features_to_simulate[1];
    if (medium_channel.getProteinIdentifications().size() > 0)
    {
      applyLabelToProteinHit_(medium_channel, medium_channel_arginine_label_, medium_channel_lysine_label_);
    }

    //check for third channel and label
    if (features_to_simulate.size() == 3)
    {
      SimTypes::FeatureMapSim& heavy_channel = features_to_simulate[2];
      if (heavy_channel.getProteinIdentifications().size() > 0)
      {
        applyLabelToProteinHit_(heavy_channel, heavy_channel_arginine_label_, heavy_channel_lysine_label_);
      }
    }
  }

  String SILACLabeler::getUnmodifiedSequence_(const Feature& feature, const String& arginine_label, const String& lysine_label) const
  {
    String unmodified_sequence = "";
    for (AASequence::ConstIterator residue = feature.getPeptideIdentifications()[0].getHits()[0].getSequence().begin();
         residue != feature.getPeptideIdentifications()[0].getHits()[0].getSequence().end();
         ++residue)
    {
      if (*residue == 'R' && residue->getModificationName() == arginine_label)
      {
        unmodified_sequence.append("R");
      }
      else if (*residue == 'K' && residue->getModificationName() == lysine_label)
      {
        unmodified_sequence.append("K");
      }
      else
      {
        unmodified_sequence.append(residue->getOneLetterCode());
      }
    }
    return unmodified_sequence;
  }

  void SILACLabeler::postDigestHook(SimTypes::FeatureMapSimVector& features_to_simulate)
  {

    SimTypes::FeatureMapSim& light_channel_features = features_to_simulate[0];
    SimTypes::FeatureMapSim& medium_channel_features = features_to_simulate[1];

    // merge the generated feature maps and create consensus
    SimTypes::FeatureMapSim final_feature_map = mergeProteinIdentificationsMaps_(features_to_simulate);

    if (features_to_simulate.size() == 2)
    {
      Map<String, Feature> unlabeled_features_index;
      for (SimTypes::FeatureMapSim::iterator unlabeled_features_iter = light_channel_features.begin();
           unlabeled_features_iter != light_channel_features.end();
           ++unlabeled_features_iter)
      {
        (*unlabeled_features_iter).ensureUniqueId();
        unlabeled_features_index.insert(std::make_pair(
                                          (*unlabeled_features_iter).getPeptideIdentifications()[0].getHits()[0].getSequence().toString()
                                                      ,
                                          *unlabeled_features_iter
                                          ));
      }

      // iterate over second map
      for (SimTypes::FeatureMapSim::iterator labeled_feature_iter = medium_channel_features.begin(); labeled_feature_iter != medium_channel_features.end(); ++labeled_feature_iter)
      {
        const String unmodified_sequence = getUnmodifiedSequence_(*labeled_feature_iter, medium_channel_arginine_label_, medium_channel_lysine_label_);

        // guarantee uniqueness
        (*labeled_feature_iter).ensureUniqueId();

        // check if we have a pair
        if (unlabeled_features_index.has(unmodified_sequence))
        {
          // own scope as we don't know what happens to 'f_modified' once we call erase() below
          Feature& unlabeled_feature = unlabeled_features_index[unmodified_sequence];
          // guarantee uniqueness
          unlabeled_feature.ensureUniqueId();

          // feature has a SILAC Label and is not equal to non-labeled
          if ((*labeled_feature_iter).getPeptideIdentifications()[0].getHits()[0].getSequence().isModified())
          {
            // add features to final map
            final_feature_map.push_back(*labeled_feature_iter);
            final_feature_map.push_back(unlabeled_feature);

            // create consensus feature
            ConsensusFeature cf;
            cf.insert(MEDIUM_FEATURE_MAPID_, *labeled_feature_iter);
            cf.insert(LIGHT_FEATURE_MAPID_, unlabeled_feature);
            cf.ensureUniqueId();
            consensus_.push_back(cf);

            // remove unlabeled feature
            unlabeled_features_index.erase(unmodified_sequence);
          }
          else
          {
            // merge features since they are equal
            Feature final_feature = mergeFeatures_(*labeled_feature_iter, unmodified_sequence, unlabeled_features_index, 1, 2);
            final_feature_map.push_back(final_feature);
          }
        }
        else // no SILAC pair, just add the labeled one
        {
          final_feature_map.push_back(*labeled_feature_iter);
        }
      }

      // add singletons from unlabeled channel
      // clean up unlabeled_index
      for (Map<String, Feature>::iterator unlabeled_index_iter = unlabeled_features_index.begin(); unlabeled_index_iter != unlabeled_features_index.end(); ++unlabeled_index_iter)
      {
        // the single ones from c0
        final_feature_map.push_back(unlabeled_index_iter->second);
      }
    }

    // merge three channels
    if (features_to_simulate.size() == 3)
    {

      // index of unlabeled channelunlabeled_feature
      Map<String, Feature> unlabeled_features_index;
      for (SimTypes::FeatureMapSim::iterator unlabeled_features_iter = light_channel_features.begin();
           unlabeled_features_iter != light_channel_features.end();
           ++unlabeled_features_iter)
      {
        (*unlabeled_features_iter).ensureUniqueId();
        unlabeled_features_index.insert(std::make_pair(
                                          (*unlabeled_features_iter).getPeptideIdentifications()[0].getHits()[0].getSequence().toString()
                                                      ,
                                          *unlabeled_features_iter
                                          ));
      }

      // index of labeled channel
      Map<String, Feature> medium_features_index;
      for (SimTypes::FeatureMapSim::iterator labeled_features_iter = medium_channel_features.begin();
           labeled_features_iter != medium_channel_features.end();
           ++labeled_features_iter)
      {
        (*labeled_features_iter).ensureUniqueId();
        medium_features_index.insert(std::make_pair(
                                       getUnmodifiedSequence_(*labeled_features_iter, medium_channel_arginine_label_, medium_channel_lysine_label_)
                                                   ,
                                       *labeled_features_iter
                                       ));
      }

      SimTypes::FeatureMapSim& heavy_labeled_features = features_to_simulate[2];
      for (SimTypes::FeatureMapSim::iterator heavy_labeled_feature_iter = heavy_labeled_features.begin();
           heavy_labeled_feature_iter != heavy_labeled_features.end();
           ++heavy_labeled_feature_iter)
      {

        Feature& heavy_feature = *heavy_labeled_feature_iter;
        heavy_feature.ensureUniqueId();

        String heavy_feature_unmodified_sequence = getUnmodifiedSequence_(heavy_feature, heavy_channel_arginine_label_, heavy_channel_lysine_label_);

        if (unlabeled_features_index.has(heavy_feature_unmodified_sequence) && medium_features_index.has(heavy_feature_unmodified_sequence))
        {
          // it is a triplet
          // c2 & c1 modified
          if (heavy_feature_unmodified_sequence != heavy_feature.getPeptideIdentifications()[0].getHits()[0].getSequence().toString())
          {
            // add features to final map
            final_feature_map.push_back(heavy_feature);
            final_feature_map.push_back(medium_features_index[heavy_feature_unmodified_sequence]);
            final_feature_map.push_back(unlabeled_features_index[heavy_feature_unmodified_sequence]);

            ConsensusFeature c_triplet;
            c_triplet.insert(HEAVY_FEATURE_MAPID_, heavy_feature);
            c_triplet.insert(LIGHT_FEATURE_MAPID_, unlabeled_features_index[heavy_feature_unmodified_sequence]);
            c_triplet.insert(MEDIUM_FEATURE_MAPID_, medium_features_index[heavy_feature_unmodified_sequence]);
            c_triplet.ensureUniqueId();

            consensus_.push_back(c_triplet);
          }
          else
          {
            // merge all three channels
            Feature completeMerge = mergeAllChannelFeatures_(heavy_feature, heavy_feature_unmodified_sequence, unlabeled_features_index, medium_features_index);
            final_feature_map.push_back(completeMerge);
          }
          // remove features from indices
          unlabeled_features_index.erase(heavy_feature_unmodified_sequence);
          medium_features_index.erase(heavy_feature_unmodified_sequence);
        }
        else if (unlabeled_features_index.has(heavy_feature_unmodified_sequence))
        {
          // 2nd case light and heavy pair
          if (heavy_feature_unmodified_sequence != heavy_feature.getPeptideIdentifications()[0].getHits()[0].getSequence().toString())
          {
            // add features to final map
            final_feature_map.push_back(heavy_feature);
            final_feature_map.push_back(unlabeled_features_index[heavy_feature_unmodified_sequence]);

            ConsensusFeature c_duplet;
            c_duplet.insert(HEAVY_FEATURE_MAPID_, heavy_feature);
            c_duplet.insert(LIGHT_FEATURE_MAPID_, unlabeled_features_index[heavy_feature_unmodified_sequence]);
            c_duplet.ensureUniqueId();

            consensus_.push_back(c_duplet);
          }
          else
          {
            // merge all three channels
            Feature completeMerge = mergeFeatures_(heavy_feature, heavy_feature_unmodified_sequence, unlabeled_features_index, 1, 3);
            final_feature_map.push_back(completeMerge);
          }
          // remove features from indices
          unlabeled_features_index.erase(heavy_feature_unmodified_sequence);
        }
        else if (medium_features_index.has(heavy_feature_unmodified_sequence))
        {
          // 3rd case medium and heavy pair
          if (heavy_feature_unmodified_sequence != heavy_feature.getPeptideIdentifications()[0].getHits()[0].getSequence().toString())
          {
            // add features to final map
            final_feature_map.push_back(heavy_feature);
            final_feature_map.push_back(medium_features_index[heavy_feature_unmodified_sequence]);

            ConsensusFeature c_duplet;
            c_duplet.insert(HEAVY_FEATURE_MAPID_, heavy_feature);
            c_duplet.insert(MEDIUM_FEATURE_MAPID_, medium_features_index[heavy_feature_unmodified_sequence]);
            c_duplet.ensureUniqueId();

            consensus_.push_back(c_duplet);
          }
          else
          {
            // merge all
            Feature completeMerge = mergeFeatures_(heavy_feature, heavy_feature_unmodified_sequence, medium_features_index, 2, 3);
            final_feature_map.push_back(completeMerge);
          }
          // remove features from indices
          medium_features_index.erase(heavy_feature_unmodified_sequence);
        }
        else
        {
          // heavy feature is a singleton
          final_feature_map.push_back(heavy_feature);
        }
      }

      // clean up labeled_index
      for (Map<String, Feature>::iterator medium_channle_index_iterator = medium_features_index.begin(); medium_channle_index_iterator != medium_features_index.end(); ++medium_channle_index_iterator)
      {
        Feature& medium_channel_feature = medium_channle_index_iterator->second;
        medium_channel_feature.ensureUniqueId();

        String medium_channel_feature_unmodified_sequence = getUnmodifiedSequence_(medium_channel_feature, medium_channel_arginine_label_, medium_channel_lysine_label_);

        if (unlabeled_features_index.has(medium_channel_feature_unmodified_sequence))
        {
          // 1. case: pair between c0 and c1
          if (medium_channel_feature.getPeptideIdentifications()[0].getHits()[0].getSequence().isModified())
          {
            // add features to final map
            final_feature_map.push_back(medium_channel_feature);
            final_feature_map.push_back(unlabeled_features_index[medium_channel_feature_unmodified_sequence]);

            ConsensusFeature c_duplet;
            c_duplet.insert(MEDIUM_FEATURE_MAPID_, medium_channel_feature);
            c_duplet.insert(LIGHT_FEATURE_MAPID_, unlabeled_features_index[medium_channel_feature_unmodified_sequence]);
            c_duplet.ensureUniqueId();
            consensus_.push_back(c_duplet);
          }
          else
          {
            // merge
            Feature completeMerge = mergeFeatures_(medium_channel_feature, medium_channel_feature_unmodified_sequence, unlabeled_features_index, 1, 2);
            final_feature_map.push_back(completeMerge);
          }
          // remove features from indices
          unlabeled_features_index.erase(medium_channel_feature_unmodified_sequence);
        }
        else
        {
          // c1 is alone
          final_feature_map.push_back(medium_channel_feature);
        }

      }

      // clean up unlabeled_index
      for (Map<String, Feature>::iterator unlabeled_index_iter = unlabeled_features_index.begin(); unlabeled_index_iter != unlabeled_features_index.end(); ++unlabeled_index_iter)
      {
        // the single ones from c0
        final_feature_map.push_back(unlabeled_index_iter->second);
      }
    }

    features_to_simulate.clear();
    features_to_simulate.push_back(final_feature_map);

    consensus_.setProteinIdentifications(final_feature_map.getProteinIdentifications());
    ConsensusMap::ColumnHeader map_description;
    map_description.label = "Simulation (Labeling Consensus)";
    map_description.size = features_to_simulate.size();
    consensus_.getColumnHeaders()[0] = map_description;
  }

  Feature SILACLabeler::mergeFeatures_(Feature& labeled_feature, const String& unmodified_feature_sequence, Map<String, Feature>& feature_index, Int index_channel_id, Int labeled_channel_id) const
  {
    // merge with feature from first map
    Feature merged_feature = feature_index[unmodified_feature_sequence];

    merged_feature.setMetaValue(getChannelIntensityName(index_channel_id), merged_feature.getIntensity());
    merged_feature.setMetaValue(getChannelIntensityName(labeled_channel_id), labeled_feature.getIntensity());

    merged_feature.setIntensity(merged_feature.getIntensity() + labeled_feature.getIntensity());

    mergeProteinAccessions_(merged_feature, labeled_feature);

    feature_index.erase(unmodified_feature_sequence);

    return merged_feature;
  }

  Feature SILACLabeler::mergeAllChannelFeatures_(Feature& heavy_channel_feature, const String& unmodified_feature_sequence, Map<String, Feature>& light_channel_feature_index, Map<String, Feature>& medium_channel_feature_index) const
  {
    // merge with feature from first map
    Feature merged_feature = light_channel_feature_index[unmodified_feature_sequence];

    merged_feature.setMetaValue(getChannelIntensityName(1), merged_feature.getIntensity());
    merged_feature.setMetaValue(getChannelIntensityName(1), medium_channel_feature_index[unmodified_feature_sequence].getIntensity());
    merged_feature.setMetaValue(getChannelIntensityName(3), heavy_channel_feature.getIntensity());

    merged_feature.setIntensity(merged_feature.getIntensity() + heavy_channel_feature.getIntensity() + medium_channel_feature_index[unmodified_feature_sequence].getIntensity());

    mergeProteinAccessions_(merged_feature, medium_channel_feature_index[unmodified_feature_sequence]);
    mergeProteinAccessions_(merged_feature, heavy_channel_feature);

    light_channel_feature_index.erase(unmodified_feature_sequence);
    medium_channel_feature_index.erase(unmodified_feature_sequence);

    return merged_feature;
  }

  bool weight_compare_less(Feature* f1, Feature* f2)
  {
    return (f1->getPeptideIdentifications())[0].getHits()[0].getSequence().getFormula().getMonoWeight()
           < (f2->getPeptideIdentifications())[0].getHits()[0].getSequence().getFormula().getMonoWeight();
  }

// TODO: rewrite
  void SILACLabeler::postRTHook(SimTypes::FeatureMapSimVector& features_to_simulate)
  {
    double rt_shift = param_.getValue("fixed_rtshift");

    // only adjust rt if we have a fixed shift
    if (rt_shift != 0.0)
    {
      // create map of all available features
      Map<UInt64, Feature*> id_map;
      SimTypes::FeatureMapSim& feature_map = features_to_simulate[0];
      for (SimTypes::FeatureMapSim::Iterator it = feature_map.begin(); it != feature_map.end(); ++it)
      {
        id_map.insert(std::make_pair<UInt64, Feature*>(it->getUniqueId(), &(*it)));
      }

      // recompute RT and set shape parameters for each consensus element
      for (ConsensusMap::Iterator consensus_it = consensus_.begin(); consensus_it != consensus_.end(); ++consensus_it)
      {
        vector<Feature*> original_features;

        // find all features that belong to this consensus element and adjust their rt
        ConsensusFeature& cf = *consensus_it;
        for (ConsensusFeature::iterator cfit = cf.begin(); cfit != cf.end(); ++cfit)
        {
          if (id_map.has(cfit->getUniqueId()))
          {
            original_features.push_back(id_map[cfit->getUniqueId()]);
          }
        }

        if (original_features.size() > 1)
        {
          std::sort(original_features.begin(), original_features.end(), weight_compare_less);

          // we use the shape parameters from the lightest fragment for all channels
          double variance = original_features[0]->getMetaValue("RT_egh_variance");
          double tau = original_features[0]->getMetaValue("RT_egh_tau");

          double startRT = original_features[0]->getRT();

          for (Size i = 0; i < original_features.size(); ++i)
          {
            original_features[i]->setRT(startRT + i * rt_shift);

            // copy shape parameters to features
            original_features[i]->setMetaValue("RT_egh_variance", variance);
            original_features[i]->setMetaValue("RT_egh_tau", tau);
          }
        }
      }
    }
  }

  void SILACLabeler::postDetectabilityHook(SimTypes::FeatureMapSimVector& /* features_to_simulate */)
  {
    // no changes to the features .. nothing todo here
  }

  void SILACLabeler::postIonizationHook(SimTypes::FeatureMapSimVector& /* features_to_simulate */)
  {
    // no changes to the features .. nothing todo here
  }

  void SILACLabeler::postRawMSHook(SimTypes::FeatureMapSimVector& features_to_simulate)
  {
    recomputeConsensus_(features_to_simulate[0]);
  }

  void SILACLabeler::postRawTandemMSHook(SimTypes::FeatureMapSimVector& /* features_to_simulate */, SimTypes::MSSimExperiment& /* simulated map */)
  {
    // no changes to the features .. nothing todo here
  }

} // namespace OpenMS
