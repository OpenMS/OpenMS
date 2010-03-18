#// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: David Wojnar $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------


#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/MATH/STATISTICS/PosteriorErrorProbabilityModel.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <vector>
#include <iostream>

using namespace OpenMS;
using namespace Math;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_IDPosteriorErrorProbability IDPosteriorErrorProbability
	
	@brief  Estimates posterior error probabilities.
	
	So far an estimation of the false score distribution with a gumbel distribution
	and the correct score distribution with a gaussian distribution is performed. 
	The probabilities are calculated using bayes law, similar to PeptideProphet.
	
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPIDPosteriorErrorProbability
	: public TOPPBase
{
 public:
	TOPPIDPosteriorErrorProbability()
		: TOPPBase("IDPosteriorErrorProbability","Estimates posterior error probabilities using a mixture model.")
	{
		
	}
	
 protected:
	void registerOptionsAndFlags_()
	{
		registerInputFile_("in","<file>","","input file ");
		setValidFormats_("in",StringList::create("idXML"));
		registerOutputFile_("out","<file>","","output file ");
	  setValidFormats_("out",StringList::create("idXML"));
	  registerDoubleOption_("smallest_e_value","<value>",10e-20,"This value gives a lower bound to E-Values. It should not be 0, as transformation in a real number (log of E-value) is not possible for certain values then.",false,true);
	  
	  registerSubsection_("fit_algorithm", "Algorithm parameter subsection");
		addEmptyLine_();	
	}
	
	//there is only one parameter at the moment
	Param getSubsectionDefaults_(const String& /*section*/) const
	{
		PosteriorErrorProbabilityModel pepm;
		return pepm.getParameters();
	}
	

	ExitCodes main_(int , const char**)
	{
		//-------------------------------------------------------------
		// parsing parameters
		//-------------------------------------------------------------
			
		String inputfile_name = getStringOption_("in");			
		String outputfile_name = getStringOption_("out");
		DoubleReal smallest_e_value = getDoubleOption_("smallest_e_value");
	
		
		//-------------------------------------------------------------
		// reading input
		//-------------------------------------------------------------
		IdXMLFile file;
		vector< ProteinIdentification > protein_ids;
		vector< PeptideIdentification > peptide_ids;
		file.load(inputfile_name, protein_ids, peptide_ids);
		
		//-------------------------------------------------------------
		// calculations
		//-------------------------------------------------------------
		vector<double> score_vector_mascot;
		vector<double> score_vector_omssa;
		vector<double> score_vector_xtandem;
		for(vector< ProteinIdentification >::iterator prot_iter = protein_ids.begin(); prot_iter < protein_ids.end(); ++prot_iter)
		{
			String engine  = prot_iter->getSearchEngine();
			for(vector< PeptideIdentification >::iterator it = peptide_ids.begin();it < peptide_ids.end(); ++it)
			{
				if(prot_iter->getIdentifier().compare(it->getIdentifier()) == 0)
				{
					vector<PeptideHit> hits = it->getHits();
					for(std::vector<PeptideHit>::iterator  hit  = hits.begin(); hit < hits.end(); ++hit)
					{
						if( engine == "OMSSA" )
						{
							score_vector_omssa.push_back((-1)* log10(max(hit->getScore(),smallest_e_value)));
						}
						else if(engine.compare("XTandem") == 0)
						{
							score_vector_xtandem.push_back((-1)* log10(max((DoubleReal)hit->getMetaValue("E-Value"),smallest_e_value)));
						}
						else if (engine.toUpper() == "MASCOT")
						{
							score_vector_mascot.push_back(hit->getScore() - (DoubleReal)hit->getMetaValue("identity_threshold"));
						}
						else
						{
							
							throw Exception::UnableToFit(__FILE__,__LINE__,__PRETTY_FUNCTION__,"No parameters for choosen search engine","The choosen search engine is currently not supported");
						}
					}
				}
			}
		}
		PosteriorErrorProbabilityModel PEP_model;
		vector<double> probabilities_mascot;
		PEP_model.fit(score_vector_mascot, probabilities_mascot);
		vector<double> probabilities_omssa;
		PEP_model.fit(score_vector_omssa, probabilities_omssa);
		vector<double> probabilities_xtandem;
		PEP_model.fit(score_vector_xtandem, probabilities_xtandem);		
		
		vector<double>::iterator prob_mascot = probabilities_mascot.begin();
		vector<double>::iterator prob_omssa = probabilities_omssa.begin();
		vector<double>::iterator prob_xtandem = probabilities_xtandem.begin();
		for(vector< ProteinIdentification >::iterator prot_iter = protein_ids.begin(); prot_iter < protein_ids.end(); ++prot_iter)
		{
			String engine  = prot_iter->getSearchEngine();		
			for(vector< PeptideIdentification >::iterator it = peptide_ids.begin();it < peptide_ids.end(); ++it)
			{
			if(prot_iter->getIdentifier().compare(it->getIdentifier()) == 0)
				{
					if( engine == "OMSSA" )
					{
						vector<PeptideHit> hits = it->getHits();
						for(std::vector<PeptideHit>::iterator  hit  = hits.begin(); hit < hits.end(); ++hit)
						{	
							hit->setMetaValue("PEP",*prob_omssa);
							++prob_omssa;
						}
						it->setHits(hits);						
					}
					else if(engine == "XTandem")
					{
						vector<PeptideHit> hits = it->getHits();
						for(std::vector<PeptideHit>::iterator  hit  = hits.begin(); hit < hits.end(); ++hit)
						{	
							hit->setMetaValue("PEP",*prob_xtandem);
							++prob_xtandem;
						}
						it->setHits(hits);
					}
					else if (engine.toUpper() == "MASCOT")
					{
						vector<PeptideHit> hits = it->getHits();
						for(std::vector<PeptideHit>::iterator  hit  = hits.begin(); hit < hits.end(); ++hit)
						{	
							hit->setMetaValue("PEP",*prob_mascot);
							++prob_mascot;
						}
						it->setHits(hits);
					}
				}
			}
		}
		//-------------------------------------------------------------
		// writing output
		//-------------------------------------------------------------
	
		file.store(outputfile_name, protein_ids, peptide_ids);
		return EXECUTION_OK;
	}
};


int main( int argc, const char** argv )
{
	TOPPIDPosteriorErrorProbability tool;

	return tool.main(argc,argv);
}

/// @endcond
