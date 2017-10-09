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

#include <OpenMS/SIMULATION/DigestSimulation.h>

#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestionLogModel.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>

namespace OpenMS
{

  DigestSimulation::DigestSimulation() :
    DefaultParamHandler("DigestSimulation")
  {
    setDefaultParams_();
  }

  DigestSimulation::DigestSimulation(const DigestSimulation& source) :
    DefaultParamHandler(source)
  {
  }

  DigestSimulation& DigestSimulation::operator=(const DigestSimulation& source)
  {
    if (this != &source)
    {
      DefaultParamHandler::operator=(source);
    }
    return *this;
  }

  DigestSimulation::~DigestSimulation()
  {
  }

  void DigestSimulation::setDefaultParams_()
  {
    // supported enzymes
    StringList enzymes;
    ProteaseDB::getInstance()->getAllNames(enzymes);
    defaults_.setValue("enzyme", "Trypsin", "Enzyme to use for digestion (select 'no cleavage' to skip digestion)");
    defaults_.setValidStrings("enzyme", enzymes);

    // cleavages
    defaults_.setValue("model", "naive", "The cleavage model to use for digestion. 'Trained' is based on a log likelihood model (see DOI:10.1021/pr060507u).");
    defaults_.setValidStrings("model", ListUtils::create<String>("trained,naive"));

    defaults_.setValue("model_trained:threshold", 0.50, "Model threshold for calling a cleavage. Higher values increase the number of cleavages. -2 will give no cleavages, +4 almost full cleavage.");
    defaults_.setMinFloat("model_trained:threshold", -2);
    defaults_.setMaxFloat("model_trained:threshold", 4);

    defaults_.setValue("model_naive:missed_cleavages", 1, "Maximum number of missed cleavages considered. All possible resulting peptides will be created.");
    defaults_.setMinInt("model_naive:missed_cleavages", 0);

    // pep length
    defaults_.setValue("min_peptide_length", 3, "Minimum peptide length after digestion (shorter ones will be discarded)");
    defaults_.setMinInt("min_peptide_length", 1);

    defaultsToParam_();
  }

  void DigestSimulation::digest(SimTypes::FeatureMapSim& feature_map)
  {
    LOG_INFO << "Digest Simulation ... started" << std::endl;

    if ((String)param_.getValue("enzyme") == String("no cleavage"))
    {
      //peptides = proteins;
      // convert all proteins into peptides

      // for each protein_hit in the FeatureMap
      for (std::vector<ProteinHit>::iterator protein_hit = feature_map.getProteinIdentifications()[0].getHits().begin();
           protein_hit != feature_map.getProteinIdentifications()[0].getHits().end();
           ++protein_hit)
      {
        // generate a PeptideHit hit with the correct link to the protein
        PeptideHit pep_hit(1.0, 1, 0, AASequence::fromString(protein_hit->getSequence()));
        PeptideEvidence pe;
        pe.setProteinAccession(protein_hit->getAccession());
        pep_hit.addPeptideEvidence(pe);

        // add the PeptideHit to the PeptideIdentification
        PeptideIdentification pep_id;
        pep_id.insertHit(pep_hit);

        // generate Feature with correct Intensity and corresponding PeptideIdentification
        Feature f;
        f.getPeptideIdentifications().push_back(pep_id);
        f.setIntensity(protein_hit->getMetaValue("intensity"));

        // copy intensity meta-values and additional annotations from Protein to Feature
        StringList keys;
        protein_hit->getKeys(keys);
        for (StringList::const_iterator it_key = keys.begin(); it_key != keys.end(); ++it_key)
        {
          f.setMetaValue(*it_key, protein_hit->getMetaValue(*it_key));
        }

        // add Feature to SimTypes::FeatureMapSim
        feature_map.push_back(f);
      }

      return;
    }


    UInt min_peptide_length = param_.getValue("min_peptide_length");
    bool use_log_model = param_.getValue("model") == "trained" ? true : false;
    UInt missed_cleavages = param_.getValue("model_naive:missed_cleavages");
    double cleave_threshold = param_.getValue("model_trained:threshold");
    if (use_log_model)
    {
      EnzymaticDigestionLogModel digestion;
      digestion.setLogThreshold(cleave_threshold);
    }
    else
    {
      ProteaseDigestion digestion;
      digestion.setEnzyme((String)param_.getValue("enzyme"));
    }

    std::vector<AASequence> digestion_products;

    // keep track of generated features
    std::map<AASequence, Feature> generated_features;

    // Iterate through ProteinHits in the FeatureMap and digest them
    for (std::vector<ProteinHit>::iterator protein_hit = feature_map.getProteinIdentifications()[0].getHits().begin();
         protein_hit != feature_map.getProteinIdentifications()[0].getHits().end();
         ++protein_hit)
    {
      // determine abundance of each digestion product (this is quite long now...)
      // we assume that each digestion product will have the same abundance
      // note: missed cleavages reduce overall abundance as they combine two (or more) single peptides

      // how many "atomic"(i.e. non-cleavable) peptides are created?
      Size complete_digest_count;
      if (use_log_model)
      {
        EnzymaticDigestionLogModel digestion;
        digestion.setLogThreshold(cleave_threshold);
        complete_digest_count = digestion.peptideCount(AASequence::fromString(protein_hit->getSequence()));
      }
      else
      {
        ProteaseDigestion digestion;
        digestion.setEnzyme((String)param_.getValue("enzyme"));
        digestion.setMissedCleavages(0);
        complete_digest_count = digestion.peptideCount(AASequence::fromString(protein_hit->getSequence()));
      }

      // compute average number of "atomic" peptides summed from all digestion products
      Size number_atomic_whole = 0;
      Size number_of_digestion_products = 0;
      for (Size i = 0; (i <= missed_cleavages) && (i < complete_digest_count); ++i)
      {
        number_atomic_whole += (complete_digest_count - i) * (i + 1);
        number_of_digestion_products += (complete_digest_count - i);
      }

      // mean number of "atomic" peptides per digestion product is now: number_atomic_whole / number_of_digestion_products
      // -> thus abundance of a digestion product is: #proteins / avg#of"atomic"peptides
      // i.e.: protein->second / (number_atomic_whole / number_of_digestion_products)

      Map<String, SimTypes::SimIntensityType> intensities;
      StringList keys;
      protein_hit->getKeys(keys);
      for (StringList::const_iterator it_key = keys.begin(); it_key != keys.end(); ++it_key)
      {
        if (!it_key->hasPrefix("intensity"))
          continue;
        intensities[*it_key] = std::max(SimTypes::SimIntensityType(1), SimTypes::SimIntensityType(protein_hit->getMetaValue(*it_key))
                                        * SimTypes::SimIntensityType(number_of_digestion_products)
                                        / SimTypes::SimIntensityType(number_atomic_whole)); // order changed for numeric stability
      }

      // do real digest
      if (use_log_model)
      {
        EnzymaticDigestionLogModel digestion;
        digestion.setLogThreshold(cleave_threshold);
        digestion.digest(AASequence::fromString(protein_hit->getSequence()), digestion_products);
      }
      else
      {
        ProteaseDigestion digestion;
        digestion.setEnzyme((String)param_.getValue("enzyme"));
        digestion.setMissedCleavages(missed_cleavages);
        digestion.digest(AASequence::fromString(protein_hit->getSequence()), digestion_products);
      }

      for (std::vector<AASequence>::const_iterator dp_it = digestion_products.begin();
           dp_it != digestion_products.end();
           ++dp_it)
      {
        if (dp_it->size() < min_peptide_length)
          continue;

        // sum equal peptide's intensities
        // *dp_it -> peptide
        // If we see this Peptide the first time -> generate corresponding feature
        if (generated_features.count(*dp_it) == 0)
        {
          PeptideHit pep_hit(1.0, 1, 0, *dp_it);

          PeptideIdentification pep_id;
          pep_id.insertHit(pep_hit);

          // create feature
          Feature f;
          f.getPeptideIdentifications().push_back(pep_id);
          // set intensity to 0 to avoid problems when summing up
          f.setIntensity(0.0);

          // copy all non-intensity meta values
          StringList lkeys;
          protein_hit->getKeys(lkeys);
          for (StringList::iterator key = lkeys.begin(); key != lkeys.end(); ++key)
          {
            if (!key->hasPrefix("intensity"))
            {
              f.setMetaValue(*key, protein_hit->getMetaValue(*key));
            }
          }

          // insert into map
          generated_features.insert(std::make_pair(*dp_it, f));
        }

        // sum up intensity values
        generated_features[*dp_it].setIntensity(generated_features[*dp_it].getIntensity() + intensities["intensity"]);
        // ... same for other intensities (iTRAQ...)
        for (Map<String, SimTypes::SimIntensityType>::const_iterator it_other = intensities.begin(); it_other != intensities.end(); ++it_other)
        {
          if (!generated_features[*dp_it].metaValueExists(it_other->first))
          {
            generated_features[*dp_it].setMetaValue(it_other->first, it_other->second);
          }
          else
          {
            generated_features[*dp_it].setMetaValue(it_other->first, SimTypes::SimIntensityType(generated_features[*dp_it].getMetaValue(it_other->first)) + it_other->second);
          }
        }

        // add current protein accession
        // existing proteins accessions ...
        std::set<String> protein_accessions = generated_features[*dp_it].getPeptideIdentifications()[0].getHits()[0].extractProteinAccessionsSet();

        // ... add accession of current protein
        protein_accessions.insert(protein_hit->getAccession());

        std::vector<PeptideIdentification> pep_idents = generated_features[*dp_it].getPeptideIdentifications();
        std::vector<PeptideHit> pep_hits = pep_idents[0].getHits();

        for (std::set<String>::const_iterator s_it = protein_accessions.begin(); s_it != protein_accessions.end(); ++s_it)
        {
          PeptideEvidence pe;
          pe.setProteinAccession(*s_it);
          pep_hits[0].addPeptideEvidence(pe);
        }
        pep_idents[0].setHits(pep_hits);
        generated_features[*dp_it].setPeptideIdentifications(pep_idents);
      }
    }

    // add generated_features to FeatureMap
    for (std::map<AASequence, Feature>::iterator it_gf = generated_features.begin();
         it_gf != generated_features.end();
         ++it_gf)
    {
      // round up intensity
      (it_gf->second).setIntensity(ceil((it_gf->second).getIntensity()));
      feature_map.push_back(it_gf->second);
    }

  }

}
