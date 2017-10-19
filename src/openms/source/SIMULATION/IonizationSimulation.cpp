// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Chris Bielow$
// $Authors: Stephan Aiche, Chris Bielow$
// --------------------------------------------------------------------------

#include <OpenMS/SIMULATION/IonizationSimulation.h>

#include <OpenMS/DATASTRUCTURES/Adduct.h>
#include <OpenMS/DATASTRUCTURES/Compomer.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <OpenMS/CONCEPT/Constants.h>

#include <cmath>
#include <algorithm>

#include <boost/bind.hpp>
#include <boost/random/binomial_distribution.hpp>
#include <boost/random/discrete_distribution.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace OpenMS
{

  IonizationSimulation::IonizationSimulation() :
    DefaultParamHandler("IonizationSimulation"),
    ProgressLogger(),
    ionization_type_(),
    basic_residues_(),
    esi_probability_(),
    esi_impurity_probabilities_(),
    esi_adducts_(),
    max_adduct_charge_(),
    maldi_probabilities_(),
    rnd_gen_(new SimTypes::SimRandomNumberGenerator())
  {
    setDefaultParams_();
    updateMembers_();
  }

  IonizationSimulation::IonizationSimulation(SimTypes::MutableSimRandomNumberGeneratorPtr random_generator) :
    DefaultParamHandler("IonizationSimulation"),
    ProgressLogger(),
    ionization_type_(),
    basic_residues_(),
    esi_probability_(),
    esi_impurity_probabilities_(),
    esi_adducts_(),
    max_adduct_charge_(),
    maldi_probabilities_(),
    rnd_gen_(random_generator)
  {
    setDefaultParams_();
    updateMembers_();
  }

  IonizationSimulation::IonizationSimulation(const IonizationSimulation& source) :
    DefaultParamHandler(source),
    ProgressLogger(source),
    ionization_type_(source.ionization_type_),
    basic_residues_(source.basic_residues_),
    esi_probability_(source.esi_probability_),
    esi_impurity_probabilities_(source.esi_impurity_probabilities_),
    esi_adducts_(source.esi_adducts_),
    max_adduct_charge_(source.max_adduct_charge_),
    maldi_probabilities_(source.maldi_probabilities_),
    rnd_gen_(source.rnd_gen_)
  {
    //updateMembers_();
  }

  IonizationSimulation& IonizationSimulation::operator=(const IonizationSimulation& source)
  {
    DefaultParamHandler::operator=(source);
    ionization_type_ = source.ionization_type_;
    basic_residues_ = source.basic_residues_;
    esi_probability_ = source.esi_probability_;
    esi_impurity_probabilities_ = source.esi_impurity_probabilities_;
    esi_adducts_ = source.esi_adducts_;
    max_adduct_charge_ = source.max_adduct_charge_;
    maldi_probabilities_ = source.maldi_probabilities_;
    rnd_gen_ = source.rnd_gen_;
    //updateMembers_();
    return *this;
  }

  IonizationSimulation::~IonizationSimulation()
  {
  }

  void IonizationSimulation::ionize(SimTypes::FeatureMapSim& features, ConsensusMap& charge_consensus, SimTypes::MSSimExperiment& experiment)
  {
    LOG_INFO << "Ionization Simulation ... started" << std::endl;

    // clear the consensus map
    charge_consensus = ConsensusMap();
    charge_consensus.setProteinIdentifications(features.getProteinIdentifications());

    switch (ionization_type_)
    {
    case MALDI:
      ionizeMaldi_(features, charge_consensus);
      break;

    case ESI:
      ionizeEsi_(features, charge_consensus);
      break;
    }

    // add params for subsequent modules
    ScanWindow sw;
    sw.begin = minimal_mz_measurement_limit_;
    sw.end = maximal_mz_measurement_limit_;
    for (Size i = 0; i < experiment.size(); ++i)
    {
      experiment[i].getInstrumentSettings().getScanWindows().push_back(sw);
    }

    ConsensusMap::FileDescription map_description;
    map_description.label = "Simulation (Charge Consensus)";
    map_description.size = features.size();
    charge_consensus.getFileDescriptions()[0] = map_description;
  }

  void IonizationSimulation::setDefaultParams_()
  {
    defaults_.setValue("ionization_type", "ESI", "Type of Ionization (MALDI or ESI)");
    defaults_.setValidStrings("ionization_type", ListUtils::create<String>("MALDI,ESI"));

    defaults_.setValue("esi:ionized_residues", ListUtils::create<String>("Arg,Lys,His"), "List of residues (as three letter code) that will be considered during ES ionization. The N-term is always assumed to carry a charge. This parameter will be ignored during MALDI ionization.");
    StringList valid_ionized_residues = ListUtils::create<String>("Ala,Cys,Asp,Glu,Phe,Gly,His,Ile,Lys,Leu,Met,Asn,Pro,Gln,Arg,Sec,Ser,Thr,Val,Trp,Tyr");
    defaults_.setValidStrings("esi:ionized_residues", valid_ionized_residues);
    defaults_.setValue("esi:charge_impurity", ListUtils::create<String>("H+:1"), "List of charged ions that contribute to charge with weight of occurrence (their sum is scaled to 1 internally), e.g. ['H:1'] or ['H:0.7' 'Na:0.3'], ['H:4' 'Na:1'] (which internally translates to ['H:0.8' 'Na:0.2'])");

    defaults_.setValue("esi:max_impurity_set_size", 3, "Maximal #combinations of charge impurities allowed (each generating one feature) per charge state. E.g. assuming charge=3 and this parameter is 2, then we could choose to allow '3H+, 2H+Na+' features (given a certain 'charge_impurity' constraints), but no '3H+, 2H+Na+, 3Na+'", ListUtils::create<String>("advanced"));

    // ionization probabilities
    defaults_.setValue("esi:ionization_probability", 0.8, "Probability for the binomial distribution of the ESI charge states");
    defaults_.setValue("maldi:ionization_probabilities", ListUtils::create<double>("0.9,0.1"), "List of probabilities for the different charge states during MALDI ionization (the list must sum up to 1.0)");

    // maximal size of map in mz dimension
    defaults_.setValue("mz:lower_measurement_limit", 200.0, "Lower m/z detector limit.");
    defaults_.setMinFloat("mz:lower_measurement_limit", 0.0);
    defaults_.setValue("mz:upper_measurement_limit", 2500.0, "Upper m/z detector limit.");
    defaults_.setMinFloat("mz:upper_measurement_limit", 0.0);

    defaultsToParam_();
  }

  void IonizationSimulation::updateMembers_()
  {
    String type = param_.getValue("ionization_type");
    if (type == "ESI")
    {
      ionization_type_ = ESI;
    }
    else if (type == "MALDI")
    {
      ionization_type_ = MALDI;
    }
    else
    {
      /// unsupported ionization model
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "IonizationSimulation got invalid Ionization type '" + type + "'");
    }

    // get basic residues from params
    basic_residues_.clear();
    StringList basic_residues = param_.getValue("esi:ionized_residues");
    for (StringList::const_iterator it = basic_residues.begin(); it != basic_residues.end(); ++it)
    {
      basic_residues_.insert(*it);
    }

    // parse possible ESI adducts
    StringList esi_charge_impurity = param_.getValue("esi:charge_impurity");
    if (esi_charge_impurity.empty())
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("IonizationSimulation got empty esi:charge_impurity! You need to specify at least one adduct (usually 'H+:1')"));
    StringList components;
    max_adduct_charge_ = 0;
    // reset internal state:
    esi_impurity_probabilities_.clear();
    esi_adducts_.clear();
    // accumulate probabilities in list
    double summed_probability(0);
    for (Size i = 0; i < esi_charge_impurity.size(); ++i)
    {
      esi_charge_impurity[i].split(':', components);
      if (components.size() != 2)
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("IonizationSimulation got invalid esi:charge_impurity (") + esi_charge_impurity[i] + ") with " + components.size() + " components instead of 2.");
      // determine charge of adduct (by # of '+')
      Size l_charge = components[0].size();
      l_charge -= components[0].remove('+').size();
      EmpiricalFormula ef(components[0].remove('+'));
      // effectively subtract electrons
      ef.setCharge(l_charge); ef -= EmpiricalFormula(String("H") + String(l_charge));
      // create adduct
      Adduct a((Int)l_charge, 1, ef.getMonoWeight(), components[0].remove('+'), log(components[1].toDouble()), 0);
      esi_adducts_.push_back(a);
      esi_impurity_probabilities_.push_back(components[1].toDouble());
      summed_probability += esi_impurity_probabilities_.back();

      max_adduct_charge_ = std::max(max_adduct_charge_, l_charge);
    }

    // scale probability to 1
    for (Size i = 0; i < esi_charge_impurity.size(); ++i)
    {
      esi_impurity_probabilities_[i] /= summed_probability;
    }

    // MALDI charge distribution
    maldi_probabilities_ = param_.getValue("maldi:ionization_probabilities");

    esi_probability_ = param_.getValue("esi:ionization_probability");

    // detector ranges
    maximal_mz_measurement_limit_ = param_.getValue("mz:upper_measurement_limit");
    minimal_mz_measurement_limit_ = param_.getValue("mz:lower_measurement_limit");

    if (minimal_mz_measurement_limit_ > maximal_mz_measurement_limit_)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "m/z measurement limits do not define a valid interval!");
    }

  }

  class IonizationSimulation::CompareCmpByEF_
  {
public:
    bool operator()(const Compomer& x, const Compomer& y) const { return x.getAdductsAsString() < y.getAdductsAsString(); }
  };

  void IonizationSimulation::ionizeEsi_(SimTypes::FeatureMapSim& features, ConsensusMap& charge_consensus)
  {
    for (size_t i = 0; i < esi_impurity_probabilities_.size(); ++i)
      std::cout << "esi_impurity_probabilities_[" << i << "]: " << esi_impurity_probabilities_.at(i) << std::endl;

    std::vector<double> weights;
    std::transform(esi_impurity_probabilities_.begin(),
                   esi_impurity_probabilities_.end(),
                   std::back_inserter(weights),
                   boost::bind(std::multiplies<double>(), _1, 10));
    for (size_t i = 0; i < weights.size(); ++i)
      std::cout << "weights[" << i << "]: " << weights.at(i) << std::endl;
    boost::random::discrete_distribution<Size, double> ddist(weights.begin(), weights.end());

    try
    {
      // map for charged features
      SimTypes::FeatureMapSim copy_map = features;
      // but leave meta information & other stuff intact
      copy_map.clear(false);

      // features which are not ionized
      Size uncharged_feature_count = 0;
      // features discarded - out of mz detection range
      Size undetected_features_count = 0;

      LOG_INFO << "Simulating " << features.size() << " features" << std::endl;

      this->startProgress(0, features.size(), "Ionization");
      Size progress(0);

      // iterate over all features
#pragma omp parallel for reduction(+: uncharged_feature_count, undetected_features_count)
      for (SignedSize index = 0; index < (SignedSize)features.size(); ++index)
      {
        // no barrier here .. only an atomic update of progress value
#pragma omp atomic
        ++progress;

#ifdef _OPENMP
        // progress logger, only master thread sets progress (no barrier here)
        if (omp_get_thread_num() == 0)
          this->setProgress(progress);
#else
        this->setProgress(progress);
#endif

        ConsensusFeature cf;

        // iterate on abundance
        Int abundance = (Int) ceil(features[index].getIntensity());
        UInt basic_residues_c = countIonizedResidues_(features[index].getPeptideIdentifications()[0].getHits()[0].getSequence());


        /// shortcut: if abundance is >1000, we 1) downsize by power of 2 until 1000 < abundance_ < 2000
        ///                                     2) dice distribution
        ///                                     3) blow abundance up to original level  (to save A LOT of computation time)
        Int power_factor_2(0);
        while (abundance > 1000)
        {
          ++power_factor_2;
          abundance /= 2;
        }

        if (basic_residues_c == 0)
        {
          ++uncharged_feature_count; // OMP
          continue;
        }

        // precompute random numbers:
        std::vector<UInt> prec_rndbin(abundance);
        {
          boost::random::binomial_distribution<Int, double> bdist(basic_residues_c, esi_probability_);
          for (Int j = 0; j < abundance; ++j)
          {
            Int rnd_no = bdist(rnd_gen_->getTechnicalRng());
            prec_rndbin[j] = (UInt) rnd_no; //cast is save because random dist should give result in the intervall [0, basic_residues_c]
          }
        }

        std::vector<Size> prec_rnduni(50); // uniform numbers container
        Size prec_rnduni_remaining(0);

        // assumption: each basic residue can hold one charged adduct
        // , we need a custom comparator, as building Compomers step by step can lead to
        // numeric diffs (and thus distinct compomers) - we only use EF to discern, that's sufficient here
        std::map<Compomer, UInt, CompareCmpByEF_> charge_states;
        Size adduct_index;
        UInt charge;

        // sample different charge states (dice for each peptide molecule separately)
        for (Int j = 0; j < abundance; ++j)
        {
          // currently we might also loose some molecules here (which is ok?)
          // sample charge state from binomial

          charge = prec_rndbin[j]; // get precomputed rnd

          if (charge == 0)
          {
            continue;
          }

          /////
          // distribute charges across adduct types
          /////
          Compomer cmp;
          // if there is only one adduct allowed (usually H+), this is easy
          if (esi_adducts_.size() == 1)
          {
            cmp.add(esi_adducts_[0] * charge, Compomer::RIGHT);
          }
          else // for more elaborate adducts
          {
            for (UInt charge_site = 0; charge_site < charge; ++charge_site)
            {
              if (prec_rnduni_remaining == 0)
              {
                // refill discrete rnd numbers if container is depleted
                {
                  for (Size i_rnd = 0; i_rnd < prec_rnduni.size(); ++i_rnd)
                  {
                    prec_rnduni[i_rnd] = ddist(rnd_gen_->getTechnicalRng());
                  }
                  prec_rnduni_remaining = prec_rnduni.size();
                }
              }
              adduct_index = prec_rnduni[--prec_rnduni_remaining];
              cmp.add(esi_adducts_[adduct_index], Compomer::RIGHT);
            }
          }

          // add 1 to abundance of sampled charge state
          ++charge_states[cmp];
        }

        // no charges > 0 selected (this should be really rare)
        if (charge_states.empty())
        {
          ++uncharged_feature_count; // OMP!
          continue;
        }

        // re-scale abundance to original value if it was below 1000
        //   -> this might lead to small numerical differences to original abundance
        UInt factor = pow(2.0, power_factor_2);
        for (std::map<Compomer, UInt, CompareCmpByEF_>::const_iterator it_m = charge_states.begin(); it_m != charge_states.end(); ++it_m)
        {
          charge_states[it_m->first] *= factor;
        }

        // transform into a set (for sorting by abundance)
        Int max_observed_charge(0);
        std::set<std::pair<UInt, Compomer> > charge_states_sorted;
        for (std::map<Compomer, UInt, CompareCmpByEF_>::const_iterator it_m = charge_states.begin(); it_m != charge_states.end(); ++it_m) // create set of pair(abundance, Compomer)
        {
          charge_states_sorted.insert(charge_states_sorted.begin(), std::make_pair(it_m->second, it_m->first));
          // update maximal observed charge
          max_observed_charge = std::max(max_observed_charge, it_m->first.getNetCharge());
        }

        Int max_compomer_types = param_.getValue("esi:max_impurity_set_size");
        std::vector<Int> allowed_entities_of_charge(max_observed_charge + 1, max_compomer_types);
        // start at highest abundant ions
        for (std::set<std::pair<UInt, Compomer> >::reverse_iterator it_s = charge_states_sorted.rbegin();
             it_s != charge_states_sorted.rend();
             ++it_s)
        {
          Int lcharge = it_s->second.getNetCharge();
          if (allowed_entities_of_charge[lcharge] > 0)
          {
            Feature charged_feature(features[index]);

            setFeatureProperties_(charged_feature, it_s->second.getMass(), it_s->second.getAdductsAsString(1), lcharge, it_s->first, index);

            // remember the original feature as parent feature (needed for labeling consensus)
            charged_feature.setMetaValue("parent_feature", String(features[index].getUniqueId()));

            if (!isFeatureValid_(charged_feature))
            {
              ++undetected_features_count; // OMP!
              continue;
            }

#pragma omp critical (OPENMS_copy_map)
            {
              copy_map.push_back(charged_feature);
            }
            // add to consensus
            cf.insert(0, charged_feature);

            // decrease # of allowed compomers of current compomer's charge
            --allowed_entities_of_charge[lcharge];
          }
        }

        // add consensus element containing all charge variants just created
#pragma omp critical (OPENMS_charge_consensus)
        {
          charge_consensus.push_back(cf);
        }

      } // ! for feature  (parallel)

      this->endProgress();

      for (Size i = 0; i < charge_consensus.size(); ++i) // this cannot be done inside the parallel-for as the copy_map might be populated meanwhile, which changes the internal uniqueid-map (used in below function)
      {
        charge_consensus[i].computeDechargeConsensus(copy_map);
      }

      // swap feature maps
      features.swap(copy_map);

      LOG_INFO << "#Peptides not ionized: " << uncharged_feature_count << std::endl;
      LOG_INFO << "#Peptides outside mz range: " << undetected_features_count << std::endl;
    }
    catch (std::exception& e)
    {
      // before leaving: free
      LOG_WARN << "Exception (" << e.what() << ") caught in " << __FILE__ << "\n";
      throw;
    }

    features.applyMemberFunction(&UniqueIdInterface::ensureUniqueId);
    charge_consensus.applyMemberFunction(&UniqueIdInterface::ensureUniqueId);
  }

  UInt IonizationSimulation::countIonizedResidues_(const AASequence& seq) const
  {
    UInt count = 1; // +1 for N-term
    for (Size i = 0; i < seq.size(); ++i)
    {
      // check for basic residues
      if (basic_residues_.count(seq[i].getShortName()) == 1)
      {
        ++count;
      }
    }

    return count;
  }

  void IonizationSimulation::ionizeMaldi_(SimTypes::FeatureMapSim& features, ConsensusMap& charge_consensus)
  {
    std::vector<double> weights;
    std::transform(maldi_probabilities_.begin(),
                   maldi_probabilities_.end(),
                   std::back_inserter(weights),
                   boost::bind(std::multiplies<double>(), _1, 10));
    boost::random::discrete_distribution<Size, double> ddist(weights.begin(), weights.end());

    try
    {
      // features discarded - out of mz detection range
      Size undetected_features_count = 0;
      Size feature_index = 0;

      SimTypes::FeatureMapSim copy_map(features);
      copy_map.clear(false);
      double h_mono_weight = Constants::PROTON_MASS_U;

      this->startProgress(0, features.size(), "Ionization");
      Size progress = 0;

      for (SignedSize index = 0; index < (SignedSize)features.size(); ++index)
      {
        Int abundance = (Int) ceil(features[index].getIntensity());
        std::vector<UInt> charge_states(maldi_probabilities_.size() + 1);
        // sample different charge states
        for (Int j = 0; j < abundance; ++j)
        {
          // sample charge from discrete distribution
          Size charge = ddist(rnd_gen_->getTechnicalRng()) + 1;

          // add 1 to abundance of sampled charge state
          ++charge_states[charge];
        }

        ConsensusFeature cf;
        // only consider charged (charge >= 1) ions
        for (UInt c = 1; c < charge_states.size(); ++c)
        {
          // empty charge states won't be generated
          if (charge_states[c] == 0)
          {
            continue;
          }
          else
          {
            Feature charged_feature(features[index]);

            setFeatureProperties_(charged_feature, h_mono_weight * c, String("H") + String(c), c, charge_states[c], feature_index);

            // remember the original feature as parent feature (needed for labeling consensus)
            charged_feature.setMetaValue("parent_feature", String(features[index].getUniqueId()));

            if (!isFeatureValid_(charged_feature))
            {
              ++undetected_features_count;
              continue;
            }

            copy_map.push_back(charged_feature);

            cf.insert(0, charged_feature);
          }
        }
        // add consensus element containing all charge variants just created
        cf.computeDechargeConsensus(copy_map);
        charge_consensus.push_back(cf);

        this->setProgress(progress);
        ++feature_index;
      } // ! feature loop (parallel)

      this->endProgress();

      // swap feature maps
      features.swap(copy_map);

      LOG_INFO << "#Peptides outside mz range: " << undetected_features_count << std::endl;
    }
    catch (std::exception& e)
    {
      LOG_WARN << "Exception (" << e.what() << ") caught in " << __FILE__ << "\n";
      throw;
    }

    features.applyMemberFunction(&UniqueIdInterface::ensureUniqueId);
    charge_consensus.applyMemberFunction(&UniqueIdInterface::ensureUniqueId);
  }

  void IonizationSimulation::setFeatureProperties_(Feature& f,
                                                   const double& adduct_mass,
                                                   const String& adduct_formula,
                                                   const SimTypes::SimChargeType charge,
                                                   const SimTypes::SimIntensityType new_intensity,
                                                   const Size parent_index)
  {
    EmpiricalFormula feature_ef = f.getPeptideIdentifications()[0].getHits()[0].getSequence().getFormula();

    f.setMZ((feature_ef.getMonoWeight() + adduct_mass) / charge);
    f.setCharge(charge);
    std::vector<PeptideHit> hits = f.getPeptideIdentifications()[0].getHits();
    hits[0].setCharge(charge);
    f.getPeptideIdentifications()[0].setHits(hits);
    // set "main" intensity
    SimTypes::SimIntensityType old_intensity = f.getIntensity();
    f.setIntensity(new_intensity);
    double factor = new_intensity / old_intensity;

#pragma omp critical (OPENMS_setfeatureprop)
    {
      // ensure uniqueness
      f.setUniqueId();
      // add meta information on compomer (mass)
      f.setMetaValue("charge_adduct_mass", adduct_mass);
      f.setMetaValue("charge_adducts", adduct_formula);
      f.setMetaValue("parent_feature_number", parent_index);

      // adapt "other" intensities (iTRAQ...) by the factor we just decreased real abundance
      StringList keys;
      f.getKeys(keys);
      for (StringList::const_iterator it_key = keys.begin(); it_key != keys.end(); ++it_key)
      {
        if (it_key->hasPrefix("intensity"))
        {
          f.setMetaValue(*it_key, SimTypes::SimIntensityType(f.getMetaValue(*it_key)) * factor);
        }
      }
    } // ! pragma
  }

  bool IonizationSimulation::isFeatureValid_(const Feature& feature)
  {
    if (feature.getMZ() > maximal_mz_measurement_limit_ || feature.getMZ() < minimal_mz_measurement_limit_) // remove feature
    {
      return false;
    }
    else
    {
      return true;
    }
  }

}
