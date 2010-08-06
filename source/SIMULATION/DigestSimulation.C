// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow$
// $Authors: Stephan Aiche, Chris Bielow$
// --------------------------------------------------------------------------

#include<OpenMS/SIMULATION/DigestSimulation.h>

#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/KERNEL/Feature.h>

namespace OpenMS {

  DigestSimulation::DigestSimulation()
    : DefaultParamHandler("DigestSimulation")
  {
    setDefaultParams_();
	}

  DigestSimulation::DigestSimulation(const DigestSimulation& source)
    : DefaultParamHandler(source)
  {
  }

  DigestSimulation& DigestSimulation::operator = (const DigestSimulation& source)
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
		enzymes.resize(EnzymaticDigestion::SIZE_OF_ENZYMES + 1);
		for (UInt i=0;i<EnzymaticDigestion::SIZE_OF_ENZYMES;++i) enzymes[i] = EnzymaticDigestion::NamesOfEnzymes[i];
		enzymes[EnzymaticDigestion::SIZE_OF_ENZYMES] = "none";
		defaults_.setValue("enzyme", enzymes[0], "Enzyme to use for digestion (select 'none' to skip digestion)");
		defaults_.setValidStrings("enzyme", enzymes);
		
		// cleavages
    defaults_.setValue("model", "trained", "The cleavage model to use for digestion. 'Trained' is based on a log likelihood model (see DOI:10.1021/pr060507u).");
    defaults_.setValidStrings("model", StringList::create("trained,naive"));

    defaults_.setValue("model_trained:threshold",0.50,"Model threshold for calling a cleavage. Higher values increase the number of cleavages. -2 will give no cleavages, +4 almost full cleavage.");
    defaults_.setMinFloat("model_trained:threshold", -2);
    defaults_.setMaxFloat("model_trained:threshold",  4);

    defaults_.setValue("model_naive:missed_cleavages",1,"Maximum number of missed cleavages considered. All possible resulting peptides will be created.");
    defaults_.setMinInt("model_naive:missed_cleavages",0);

		// pep length
		defaults_.setValue("min_peptide_length",3,"Minimum peptide length after digestion (shorter ones will be discarded)");
    defaults_.setMinInt("min_peptide_length",1);

		defaultsToParam_();		    
  }

  void DigestSimulation::digest(FeatureMapSim & feature_map)
  {
		if ((String)param_.getValue("enzyme") == String("none"))
		{
      //peptides = proteins;
      // convert all proteins into peptides
      
      // for each protein_hit in the FeatureMap
      for(std::vector<ProteinHit>::iterator protein_hit = feature_map.getProteinIdentifications()[0].getHits().begin();
          protein_hit != feature_map.getProteinIdentifications()[0].getHits().end();
          ++protein_hit)
      {
        // generate a PeptideHit hit with the correct link to the protein      
        PeptideHit pep_hit(1.0, 1, 0, protein_hit->getSequence());
        std::vector<String> prot_accessions;
        prot_accessions.push_back(protein_hit->getAccession());
        pep_hit.setProteinAccessions(prot_accessions);

        // add the PeptideHit to the PeptideIdentification
        PeptideIdentification pep_id;
        pep_id.insertHit(pep_hit);
        
        // generate Feature with correct Intensity and corresponding PeptideIdentification
        Feature f;        
        f.getPeptideIdentifications().push_back(pep_id);
        f.setIntensity(protein_hit->getMetaValue("intensity"));
				// copy intensity meta-values from Protein to Feature
				StringList keys;
				protein_hit->getKeys(keys);
				for (StringList::const_iterator it_key = keys.begin(); it_key != keys.end(); it_key++)
				{
					if (!it_key->hasPrefix("intensity")) continue;
					f.setMetaValue(*it_key, protein_hit->getMetaValue(*it_key)); 
				}
			        
        // add Feature to FeatureMapSim
        feature_map.push_back(f);     
      }
      
			return;
		}


    UInt min_peptide_length     = param_.getValue("min_peptide_length");
    bool use_log_model          = param_.getValue("model") == "trained" ? true : false;
		UInt missed_cleavages       = param_.getValue("model_naive:missed_cleavages");
    DoubleReal cleave_threshold = param_.getValue("model_trained:threshold");
    
    EnzymaticDigestion digestion;
    digestion.setEnzyme(digestion.getEnzymeByName((String)param_.getValue("enzyme")));
    digestion.setLogModelEnabled(use_log_model);
    digestion.setLogThreshold(cleave_threshold);

		std::vector<AASequence> digestion_products;

    // keep track of generated features
    std::map<AASequence, Feature> generated_features;
    
    // Iterate through ProteinHits in the FeatureMap and digest them
    for(std::vector<ProteinHit>::iterator protein_hit = feature_map.getProteinIdentifications()[0].getHits().begin();
        protein_hit != feature_map.getProteinIdentifications()[0].getHits().end();
        ++protein_hit)
    {
			// determine abundance of each digestion product (this is quite long now...)
			// we assume that each digestion product will have the same abundance
			// note: missed cleavages reduce overall abundance as they combine two (or more) single peptides

			// how many "atomic"(i.e. non-cleavable) peptides are created?
			digestion.setMissedCleavages( 0 );
			Size complete_digest_count = digestion.peptideCount(protein_hit->getSequence());
			// compute average number of "atomic" peptides summed from all digestion products
			Size number_atomic_whole = 0;
			Size number_of_digestion_products = 0;
			for (Size i=0; (i<=missed_cleavages) && (i<complete_digest_count);++i)
			{
				number_atomic_whole += (complete_digest_count-i)*(i+1);
				number_of_digestion_products += (complete_digest_count-i);
			}

			// mean number of "atomic" peptides per digestion product is now: number_atomic_whole / number_of_digestion_products
			// -> thus abundance of a digestion product is: #proteins / avg#of"atomic"peptides
			//																				i.e.: protein->second / (number_atomic_whole / number_of_digestion_products)
			
			Map<String, SimIntensityType> intensities;
			StringList keys;
			protein_hit->getKeys(keys);
			for (StringList::const_iterator it_key = keys.begin(); it_key != keys.end(); it_key++)
			{
				if (!it_key->hasPrefix("intensity")) continue;
				intensities[*it_key] = std::max(SimIntensityType(1), SimIntensityType(protein_hit->getMetaValue(*it_key)) 
																																*SimIntensityType(number_of_digestion_products) 
																																/SimIntensityType(number_atomic_whole) ); // order changed for numeric stability
			}
			
			// do real digest
			digestion.setMissedCleavages( missed_cleavages );
      digestion.digest(AASequence(protein_hit->getSequence()), digestion_products);
			
			for (std::vector<AASequence>::const_iterator dp_it = digestion_products.begin();
           dp_it != digestion_products.end();
           ++dp_it)
      {
        if (dp_it->size() < min_peptide_length) continue;
				
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
          
          // insert into map
          generated_features.insert(std::make_pair(*dp_it, f));
        }
        
        // sum up intensity values
        generated_features[*dp_it].setIntensity(generated_features[*dp_it].getIntensity() + intensities["intensity"]);
        // ... same for other intensities (iTRAQ...)
        for (Map<String, SimIntensityType>::const_iterator it_other=intensities.begin(); it_other!=intensities.end(); ++it_other)
        {
					if (!generated_features[*dp_it].metaValueExists(it_other->first))
					{
						generated_features[*dp_it].setMetaValue(it_other->first, it_other->second);
					}
					else
					{
						generated_features[*dp_it].setMetaValue(it_other->first, SimIntensityType(generated_features[*dp_it].getMetaValue(it_other->first)) + it_other->second);
					}
        }

        // add current protein accession
        // existing proteins accessions ...
        std::vector<String> prot_accessions(generated_features[*dp_it].getPeptideIdentifications()[0].getHits()[0].getProteinAccessions());

        // ... add accession of current protein
        prot_accessions.push_back(protein_hit->getAccession());

        std::vector<PeptideIdentification> pep_idents = generated_features[*dp_it].getPeptideIdentifications();
        std::vector<PeptideHit> pep_hits = pep_idents[0].getHits();
        pep_hits[0].setProteinAccessions(prot_accessions);
        pep_idents[0].setHits(pep_hits);
        generated_features[*dp_it].setPeptideIdentifications(pep_idents);
      }
		}
    
    // add generated_features to FeatureMap
    for(std::map<AASequence, Feature>::iterator it_gf = generated_features.begin();
        it_gf != generated_features.end();
        ++it_gf)
    {
      // round up intensity
      (it_gf->second).setIntensity(ceil((it_gf->second).getIntensity()));
      feature_map.push_back(it_gf->second);
    }

  }

}
