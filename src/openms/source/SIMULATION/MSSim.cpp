// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Stephan Aiche$
// $Authors: Stephan Aiche, Chris Bielow$
// --------------------------------------------------------------------------

#include <OpenMS/SIMULATION/MSSim.h>

#include <OpenMS/SIMULATION/DigestSimulation.h>
#include <OpenMS/SIMULATION/DetectabilitySimulation.h>
#include <OpenMS/SIMULATION/RawMSSignalSimulation.h>
#include <OpenMS/SIMULATION/RawTandemMSSignalSimulation.h>
#include <OpenMS/SIMULATION/IonizationSimulation.h>
#include <OpenMS/SIMULATION/RTSimulation.h>

#include <OpenMS/SIMULATION/LABELING/BaseLabeler.h>

#include <OpenMS/DATASTRUCTURES/ListUtils.h>

//#define OPENMS_DEBUG_SIM_

using namespace std;

namespace OpenMS
{

  void verbosePrintFeatureMap(FeatureMapSimVector feature_maps, String stage)
  {
#ifdef OPENMS_DEBUG_SIM_
    cout << "############## DEBUG (" << stage << ") -- FEATURE MAPS ##############" << endl;

    Size map_count = 1;
    StringList keys;
    for (FeatureMapSimVector::iterator map_iter = feature_maps.begin(); map_iter != feature_maps.end(); ++map_iter)
    {
      cout << "FEATURE MAP #" << map_count << endl;

      cout << "contained proteins" << endl;
      ProteinIdentification protIdent = (*map_iter).getProteinIdentifications()[0];
      for (vector<ProteinHit>::iterator proteinHit = protIdent.getHits().begin();
           proteinHit != protIdent.getHits().end();
           ++proteinHit)
      {
        cout << "- " << proteinHit->getAccession() << endl;
      }
      cout << "----------------------------------------------" << endl;
      for (FeatureMapSim::const_iterator feat = (*map_iter).begin();
           feat != (*map_iter).end();
           ++feat)
      {
        (*feat).getKeys(keys);
        cout << " RT: " << (*feat).getRT()
             << " MZ: " << (*feat).getMZ()
             << " INT: " << (*feat).getIntensity()
             << " CHARGE: " << (*feat).getCharge()
             << " Det: " << (*feat).getMetaValue("detectibility") << endl
             << " Pep: " << feat->getPeptideIdentifications()[0].getHits()[0].getSequence().toString() << endl
             << " ID: " << feat->getUniqueId() << endl
             << " Meta: " << keys.concatenate(",") << endl;
        cout << "derived from protein(s): ";
        for (vector<String>::const_iterator it = (*feat).getPeptideIdentifications()[0].getHits()[0].getProteinAccessions().begin();
             it != (*feat).getPeptideIdentifications()[0].getHits()[0].getProteinAccessions().end();
             ++it)
        {
          cout << (*it) << " ";
        }
        cout << endl << "----------------------------------------------" << endl;
      }
      cout << endl;
      ++map_count;
    }
    cout << "############## END DEBUG -- FEATURE MAPS ##############" << endl;
#else
    if (feature_maps.empty())
      cout << stage;                             // just to avoid warnings of unused parameters
#endif
  }

  MSSim::MSSim() :
    DefaultParamHandler("MSSim"),
    experiment_(),
    feature_maps_(),
    consensus_map_(),
    labeler_(0)
  {
    // section params
    defaults_.insert("Digestion:", DigestSimulation().getDefaults());
    defaults_.insert("RT:", RTSimulation().getDefaults());
    defaults_.insert("Detectability:", DetectabilitySimulation().getDefaults());
    defaults_.insert("Ionization:", IonizationSimulation().getDefaults());
    defaults_.insert("RawSignal:", RawMSSignalSimulation().getDefaults());
    defaults_.insert("RawTandemSignal:", RawTandemMSSignalSimulation().getDefaults());

    subsections_.push_back("Labeling");

    //sync params (remove duplicates from modules and put them in a global module)
    syncParams_(defaults_, true);
    defaultsToParam_();
  }

  MSSim::~MSSim()
  {
    delete labeler_;
  }

  Param MSSim::getParameters() const
  {
    Param tmp;
    tmp.insert("", this->param_); // get non-labeling options

    vector<String> products = Factory<BaseLabeler>::registeredProducts();

    tmp.setValue("Labeling:type", "labelfree", "Select the labeling type you want for your experiment");
    tmp.setValidStrings("Labeling:type", products);

    for (vector<String>::iterator product_name = products.begin(); product_name != products.end(); ++product_name)
    {
      BaseLabeler * labeler = Factory<BaseLabeler>::create(*product_name);
      tmp.insert("Labeling:" + *product_name + ":", labeler->getDefaultParameters());
      if (!tmp.copy("Labeling:" + *product_name).empty())
      {
        // if parameters of labeler are empty, the section will not exist and
        // the command below would fail
        tmp.setSectionDescription("Labeling:" + *product_name, labeler->getDescription());
      }
      delete(labeler);
    }

    return tmp;
  }

  void MSSim::simulate(MutableSimRandomNumberGeneratorPtr rnd_gen, SampleChannels & channels)
  {
    /*todo: move to a global config file or into INI file */
    Log_fatal.setPrefix("%S: ");
    Log_error.setPrefix("%S: ");
    Log_warn.setPrefix("%S: ");
    Log_info.setPrefix("%S: ");
    Log_debug.setPrefix("%S: ");

    /*
      General progress should be
        1. digest Proteins
        2. predict retention times
        3. predict detectability
        4. simulate ionization
        5. simulate the ms signal
        6. select features for MS2
        7. generate MS2 signals for selected features
     */

    // re-distribute synced parameters:
    //param_.store("c:/mssim_param.ini"); // test reconstruction
    syncParams_(param_, false);

    // instantiate and pass params before doing any actual work
    // ... this way, each module can throw an Exception when the parameters
    // ... do not fit, and the users gets an immediate feedback
    DigestSimulation digest_sim;
    digest_sim.setParameters(param_.copy("Digestion:", true));
    RTSimulation rt_sim(rnd_gen);
    rt_sim.setParameters(param_.copy("RT:", true));
    DetectabilitySimulation dt_sim;
    dt_sim.setParameters(param_.copy("Detectability:", true));
    IonizationSimulation ion_sim(rnd_gen);
    ion_sim.setParameters(param_.copy("Ionization:", true));
    ion_sim.setLogType(this->getLogType());
    RawMSSignalSimulation raw_sim(rnd_gen);
    raw_sim.setParameters(param_.copy("RawSignal:", true));
    raw_sim.setLogType(this->getLogType());
    raw_sim.loadContaminants(); // check if the file is valid (if not, an error is raised here instead of half-way through simulation)


    String labeling = param_.getValue("Labeling:type");
    labeler_ = Factory<BaseLabeler>::create(labeling);
    Param labeling_parameters = param_.copy("Labeling:" + labeling + ":", true);
    labeler_->setParameters(labeling_parameters);
    labeler_->setRnd(rnd_gen);

    // check parameters ..
    labeler_->preCheck(param_);

    // convert sample proteins into an empty FeatureMap with ProteinHits
    for (SampleChannels::const_iterator channel_iterator = channels.begin(); channel_iterator != channels.end(); ++channel_iterator)
    {
      FeatureMapSim map;
      createFeatureMap_(*channel_iterator, map, feature_maps_.size());
      feature_maps_.push_back(map);
    }

    // Call setUpHook
    labeler_->setUpHook(feature_maps_);

    // digest
    for (FeatureMapSimVector::iterator map_iterator = feature_maps_.begin(); map_iterator != feature_maps_.end(); ++map_iterator)
    {
      digest_sim.digest(*map_iterator);
    }

    // post digest labeling
    labeler_->postDigestHook(feature_maps_);

    // debug
    verbosePrintFeatureMap(feature_maps_, "digested");

    // RT prediction
    for (FeatureMapSimVector::iterator map_iterator = feature_maps_.begin(); map_iterator != feature_maps_.end(); ++map_iterator)
    {
      rt_sim.predictRT(*map_iterator);
    }
    rt_sim.createExperiment(experiment_);

    peak_map_ = experiment_; // initial ground-truth for peak map is the same as for raw data

    // post rt sim labeling
    labeler_->postRTHook(feature_maps_);

    // debug
    verbosePrintFeatureMap(feature_maps_, "RT sim done");

    // Detectability prediction
    for (FeatureMapSimVector::iterator map_iterator = feature_maps_.begin(); map_iterator != feature_maps_.end(); ++map_iterator)
    {
      dt_sim.filterDetectability(*map_iterator);
    }

    // post detectability labeling
    labeler_->postDetectabilityHook(feature_maps_);

    // debug
    verbosePrintFeatureMap(feature_maps_, "DT sim done");

    // at this point all feature maps should be combined to one?
    ion_sim.ionize(feature_maps_.front(), consensus_map_, experiment_);

    // post ionization labeling
    labeler_->postIonizationHook(feature_maps_);

    // debug
    verbosePrintFeatureMap(feature_maps_, "ION sim done");

    raw_sim.generateRawSignals(feature_maps_.front(), experiment_, peak_map_, contaminants_map_);

    // post raw sim labeling
    labeler_->postRawMSHook(feature_maps_);

    // debug
    verbosePrintFeatureMap(feature_maps_, "RawSignal sim done");

    RawTandemMSSignalSimulation raw_tandemsim(rnd_gen);
    raw_tandemsim.setParameters(param_.copy("RawTandemSignal:", true));
    raw_tandemsim.generateRawTandemSignals(feature_maps_.front(), experiment_, peak_map_);

    labeler_->postRawTandemMSHook(feature_maps_, experiment_);


    // some last fixing of meta-values (this is impossible to do before as we do not know the final number of scans)
    for (Size i = 0; i < feature_maps_[0].size(); ++i)
    {
      Feature & f = feature_maps_[0][i];
      PeptideIdentification & pi = f.getPeptideIdentifications()[0];
      // search for closest scan index:
      MSSimExperiment::ConstIterator it_rt = experiment_.RTBegin(f.getRT());
      SignedSize scan_index = distance<MSSimExperiment::ConstIterator>(experiment_.begin(), it_rt);
      pi.setMetaValue("RT_index", scan_index);
      pi.setMetaValue("RT", f.getRT());
    }

    LOG_INFO << "Final number of simulated features: " << feature_maps_[0].size() << "\n";

    // re-index spectra to avoid naming conflicts
    Size id = 1;
    experiment_.sortSpectra();
    peak_map_.sortSpectra();
    if (experiment_.size() != peak_map_.size())
    {
      throw Exception::InvalidSize(__FILE__, __LINE__, __PRETTY_FUNCTION__, peak_map_.size() - experiment_.size());
    }
    for (MSSimExperiment::Iterator it_e = experiment_.begin(), it_ep = peak_map_.begin(); it_e != experiment_.end(); ++it_e, ++it_ep)
    {
      String spec_id = String("scan=") + id++;
      it_e->setNativeID(spec_id);
      it_ep->setNativeID(spec_id);
    }
  }

  void MSSim::createFeatureMap_(const SampleProteins & proteins, FeatureMapSim & feature_map, Size map_index)
  {
    // clear feature map
    feature_map.clear(true);
    ProteinIdentification protIdent;

    for (SampleProteins::const_iterator it = proteins.begin(); it != proteins.end(); ++it)
    {
      // add new ProteinHit to ProteinIdentification
      ProteinHit protHit(0.0, 1, (it->first).identifier, (it->first).sequence);
      // copy all meta values from FASTA file parsing
      protHit = (it->second);
      // additional meta values:
      protHit.setMetaValue("description", it->first.description);
      protHit.setMetaValue("map_index", map_index);
      protIdent.insertHit(protHit);
    }

    vector<ProteinIdentification> vec_protIdent;
    vec_protIdent.push_back(protIdent);
    feature_map.setProteinIdentifications(vec_protIdent);
  }

  void MSSim::syncParams_(Param & p, bool to_outer)
  {
    vector<StringList> globals;
    // here the globals params are listed that require to be in sync across several modules
    // - first the global param name and following that the module names where this param occurs
    // - Warning: the module params must have unchanged names and restrictions! (descriptions can differ though)
    globals.push_back(ListUtils::create<String>("ionization_type,Ionization,RawSignal,RawTandemSignal"));

    String global_prefix = "Global";
    // remove or add local params
    if (to_outer) // remove local params and merge to global
    {
      for (Size i = 0; i < globals.size(); ++i)
      {
        // set the global param:
        OPENMS_PRECONDITION(globals[i].size() >= 2, "Param synchronisation aborting due to missing local parameters!");
        p.insert(global_prefix + ":" + globals[i][0], p.copy(globals[i][1] + ":" + globals[i][0], true));
        // remove local params
        for (Size i_local = 1; i_local < globals[i].size(); ++i_local)
        {
          p.remove(globals[i][i_local] + ":" + globals[i][0]);
        }
      }
    }
    else     // restore local params from global one
    {
      for (Size i = 0; i < globals.size(); ++i)
      {
        // get the global param:
        OPENMS_PRECONDITION(globals[i].size() >= 2, "Param synchronisation aborting due to missing local parameters!");

        Param p_global = p.copy(global_prefix + ":" + globals[i][0], true);
        // insert into local params
        for (Size i_local = 1; i_local < globals[i].size(); ++i_local)
        {
          p.insert(globals[i][i_local] + ":" + globals[i][0], p_global);
        }
      }
    }

  }

  void MSSim::updateMembers_()
  {
  }

  const MSSimExperiment& MSSim::getExperiment() const
  {
    return experiment_;
  }

  const FeatureMapSim& MSSim::getSimulatedFeatures() const
  {
    OPENMS_PRECONDITION(feature_maps_.size() == 1, "More than one feature map remains after simulation. The channels should however be merged by now. Check!")
    return feature_maps_[0];
  }

  ConsensusMap& MSSim::getChargeConsensus()
  {
    return consensus_map_;
  }

  ConsensusMap& MSSim::getLabelingConsensus()
  {
    return labeler_->getConsensus();
  }

  const FeatureMapSim& MSSim::getContaminants() const
  {
    return contaminants_map_;
  }

  const MSSimExperiment& MSSim::getPeakMap() const
  {
    return peak_map_;
  }

  void MSSim::getMS2Identifications(vector<ProteinIdentification>& proteins,
                                    vector<PeptideIdentification>& peptides)
    const
  {
    proteins.clear();
    peptides.clear();
    set<String> accessions;
    for (MSSimExperiment::const_iterator ms_it = experiment_.begin();
         ms_it != experiment_.end(); ++ms_it)
    {
      if (ms_it->getMSLevel() != 2) continue;
      // "the" precursor is the one with highest intensity:
      Size index = 0;
      DoubleReal intensity = ms_it->getPrecursors()[0].getIntensity();
      for (Size i = 1; i < ms_it->getPrecursors().size(); ++i)
      {
        if (ms_it->getPrecursors()[i].getIntensity() > intensity)
        {
          intensity = ms_it->getPrecursors()[i].getIntensity();
          index = i;
        }
      }
      IntList feat_ids = ms_it->getMetaValue("parent_feature_ids");
      const Feature& feature = feature_maps_[0][feat_ids[index]];
      peptides.push_back(feature.getPeptideIdentifications()[0]);
      peptides.back().setMetaValue("RT", ms_it->getRT());
      peptides.back().setMetaValue("MZ", ms_it->getPrecursors()[index].getMZ());
      const PeptideHit& hit = peptides.back().getHits()[0];
      accessions.insert(hit.getProteinAccessions().begin(),
                        hit.getProteinAccessions().end());
    }
    const ProteinIdentification& protein =
      feature_maps_[0].getProteinIdentifications()[0];
    proteins.push_back(protein);
    proteins[0].getHits().clear();
    for (vector<ProteinHit>::const_iterator prot_it = protein.getHits().begin();
         prot_it != protein.getHits().end(); ++prot_it)
    {
      if (accessions.find(prot_it->getAccession()) != accessions.end())
      {
        proteins[0].insertHit(*prot_it);
      }
    }
  }

}
