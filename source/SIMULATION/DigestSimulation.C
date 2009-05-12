// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
  
  void DigestSimulation::digest(FeatureMapSim & feature_map)
  {
		if ((String)param_.getValue("enzyme") == String("none"))
		{
      ProteinIdentification protIdent = feature_map.getProteinIdentifications()[0];
      // for each ProteinHit in the FeatureMap
      for(std::vector<ProteinHit>::iterator proteinHit = protIdent.getHits().begin();
          proteinHit != protIdent.getHits().end();
          ++proteinHit)
      {
        // generate a PeptideHit hit with the correct link to the protein      
        PeptideHit pep_hit(1.0, 1, 1, proteinHit->getSequence());
        std::vector<String> prot_accessions(1);
        prot_accessions.push_back(proteinHit->getAccession());
        pep_hit.setProteinAccessions(prot_accessions);

        // add the PeptideHit to the PeptideIdentification
        PeptideIdentification pep_id;
        pep_id.insertHit(pep_hit);
        
        // generate Feature with correct Intensity and corresponding PeptideIdentification
        Feature f;        
        f.getPeptideIdentifications().push_back(pep_id);
        f.setIntensity(proteinHit->getMetaValue("intesity"));
        
        // add Feature to FeatureMapSim
        feature_map.push_back(f);     
      }
      
			//peptides = proteins;
      // convert all proteins into peptides
			return;
		}

		EnzymaticDigestion digestion;
		digestion.setEnzyme(digestion.getEnzymeByName((String)param_.getValue("enzyme")));

    UInt min_peptide_length = param_.getValue("min_peptide_length");
		UInt missed_cleavages   = param_.getValue("missed_cleavages");
		
		std::vector<AASequence> digestion_products;

    // keep track of generated features
    std::map<AASequence, Feature> generated_features;
    
    // Iterate through ProteinHits in the FeatureMap and digest them
		/*
    for (SampleProteins::const_iterator protein = proteins.begin();
         protein != proteins.end();
         ++protein)*/
    for(std::vector<ProteinHit>::iterator proteinHit = feature_map.getProteinIdentifications()[0].getHits().begin();
        proteinHit != feature_map.getProteinIdentifications()[0].getHits().end();
        ++proteinHit)
    {
			// determine abundance of each digestion product (this is quite long now...)
			// we assume that each digestion product will have the same abundance
			// note: missed cleavages reduce overall abundance as they combine two (or more) single peptides

			// how many "atomic"(i.e. non-cleavable) peptides are created?
			digestion.setMissedCleavages( 0 );
			Size complete_digest_count = digestion.peptideCount(proteinHit->getSequence());
			// compute average numer of "atomic" peptides summed from all digestion products
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
			SimIntensityType abundance = std::max(SimIntensityType(1), SimIntensityType(proteinHit->getMetaValue("intensity")) 
																																*SimIntensityType(number_of_digestion_products) 
																																/SimIntensityType(number_atomic_whole) ); // order changed for numeric stability
			
			// do real digest
			digestion.setMissedCleavages( missed_cleavages );
      digestion.digest(AASequence(proteinHit->getSequence()), digestion_products);
			
			for (std::vector<AASequence>::const_iterator dp_it = digestion_products.begin();
           dp_it != digestion_products.end();
           ++dp_it)
      {
        if (dp_it->size() < min_peptide_length) continue;
				
				// sum equal peptide's intensities
        // *dp_it -> peptide
        // If we see this Peptide the first time -> generate corresponding feature
        if(generated_features.count(*dp_it) == 0) {
          PeptideHit pep_hit(1.0, 1, 1, *dp_it);

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
        
        // sum up intesity values
        generated_features[*dp_it].setIntensity(generated_features[*dp_it].getIntensity() + abundance);

        // add current protein accession
        // generate new vector of accessions
        std::vector<String> prot_accessions;

        // copy old values into new vector
        prot_accessions.assign(generated_features[*dp_it].getPeptideIdentifications()[0].getHits()[0].getProteinAccessions().begin(), 
                               generated_features[*dp_it].getPeptideIdentifications()[0].getHits()[0].getProteinAccessions().end());
        prot_accessions.push_back(proteinHit->getAccession());

        std::vector<PeptideIdentification> pepIdents = generated_features[*dp_it].getPeptideIdentifications();
        std::vector<PeptideHit> pepHits = pepIdents[0].getHits();
        pepHits[0].setProteinAccessions(prot_accessions);
        pepIdents[0].setHits(pepHits);
        generated_features[*dp_it].setPeptideIdentifications(pepIdents);
      }
		}
    
    // add generated_features to FeatureMap
    for(std::map<AASequence, Feature>::iterator gfIt = generated_features.begin();
        gfIt != generated_features.end();
        ++gfIt)
    {
      // TODO: document the rounding
      (gfIt->second).setIntensity(ceil((gfIt->second).getIntensity()));
      feature_map.push_back(gfIt->second);
    }

  }
  
  void DigestSimulation::setDefaultParams_()
  {
		// supported enzymes
		StringList enzymes;
		enzymes.resize(EnzymaticDigestion::SIZE_OF_ENZYMES + 1);
		for (UInt i=0;i<EnzymaticDigestion::SIZE_OF_ENZYMES;++i) enzymes[i] = EnzymaticDigestion::NamesOfEnzymes[i];
		enzymes[EnzymaticDigestion::SIZE_OF_ENZYMES] = "none";
		defaults_.setValue("enzyme", enzymes[0], "Enzyme to use for digestion");
		defaults_.setValidStrings("enzyme", enzymes);
		
		// cleavages
		defaults_.setValue("missed_cleavages",1,"maximum number of missed cleavages");
    
		// pep length
		defaults_.setValue("min_peptide_length",3,"minimum peptide length after digestion");
    
		defaultsToParam_();		    
  }
}
