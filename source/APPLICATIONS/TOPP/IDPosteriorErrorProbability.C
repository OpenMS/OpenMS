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

using namespace OpenMS;
using namespace Math; //PosteriorErrorProbabilityModel
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
	At the moment, it is able to handle Xtandem, Mascot and OMSSA scores.
	
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
	  registerFlag_("split_charge", "The search enginge scores are splitted by charge if this flag is set. Thus, for each charge state a new model will be computed.");
	  
	  registerSubsection_("fit_algorithm", "Algorithm parameter subsection");
		addEmptyLine_();	
	}	
	
	//there is only one parameter at the moment
	Param getSubsectionDefaults_(const String& /*section*/) const
	{
		PosteriorErrorProbabilityModel pepm;
		return pepm.getParameters();
	}
	
	double get_score_(String& engine, const PeptideHit& hit)
	{
		if( engine == "OMSSA" )
		{
			return ((-1)* log10(max(hit.getScore(),smallest_e_value_)));
		}
		else if(engine.compare("XTandem") == 0)
		{
			return((-1)* log10(max((DoubleReal)hit.getMetaValue("E-Value"),smallest_e_value_)));
		}
		else if (engine == "MASCOT")
		{
			return((-1)* log10(max((DoubleReal)hit.getMetaValue("EValue"),smallest_e_value_)));
		}
		else
		{		
			throw Exception::UnableToFit(__FILE__,__LINE__,__PRETTY_FUNCTION__,"No parameters for choosen search engine","The choosen search engine is currently not supported");
		}
	}

	ExitCodes main_(int , const char**)
	{
		//-------------------------------------------------------------
		// parsing parameters
		//-------------------------------------------------------------
			
		String inputfile_name = getStringOption_("in");			
		String outputfile_name = getStringOption_("out");
	  smallest_e_value_ = getDoubleOption_("smallest_e_value");
		Param fit_algorithm = getParam_().copy("fit_algorithm:",true);
		bool split_charge = getFlag_("split_charge");
		
		//-------------------------------------------------------------
		// reading input
		//-------------------------------------------------------------
		IdXMLFile file;
		vector< ProteinIdentification > protein_ids;
		vector< PeptideIdentification > peptide_ids;
		file.load(inputfile_name, protein_ids, peptide_ids);
		vector<double> scores;
		vector<double>probabilities;
		PosteriorErrorProbabilityModel PEP_model;
		PEP_model.setParameters(fit_algorithm);
		StringList search_engines = StringList::create("XTandem,OMSSA,MASCOT");
		//-------------------------------------------------------------
		// calculations
		//-------------------------------------------------------------
		
		if(split_charge)
		{
			vector<Int> charges;
			for(vector< PeptideIdentification >::iterator it = peptide_ids.begin();it < peptide_ids.end(); ++it)
			{
				vector<PeptideHit> hits = it->getHits();
				for(std::vector<PeptideHit>::iterator  hit  = hits.begin(); hit < hits.end(); ++hit)
				{
					if(charges.end() == find(charges.begin(), charges.end(), hit->getCharge()))
					{
						charges.push_back(hit->getCharge());
					}
				}
			}
			for(vector<Int>::iterator charge = charges.begin(); charge < charges.end(); ++charge)
			{
				for(StringList::iterator engine = search_engines.begin(); engine < search_engines.end(); ++engine)
				{
					for(vector< ProteinIdentification >::iterator prot_iter = protein_ids.begin(); prot_iter < protein_ids.end(); ++prot_iter)
					{
						String searchengine_toUpper =  prot_iter->getSearchEngine();
						searchengine_toUpper.toUpper();
						if(*engine == prot_iter->getSearchEngine() || *engine == searchengine_toUpper)
						{
							for(vector< PeptideIdentification >::iterator it = peptide_ids.begin();it < peptide_ids.end(); ++it)
							{
								if(prot_iter->getIdentifier().compare(it->getIdentifier()) == 0)
								{
									vector<PeptideHit> hits = it->getHits();
									for(std::vector<PeptideHit>::iterator  hit  = hits.begin(); hit < hits.end(); ++hit)
									{
										if(hit->getCharge() == *charge)
										{
											scores.push_back(get_score_(*engine, *hit));
										}
									}
								}
							}
						}
					}
					PEP_model.fit(scores, probabilities);
					for(vector< ProteinIdentification >::iterator prot_iter = protein_ids.begin(); prot_iter < protein_ids.end(); ++prot_iter)
					{
						String searchengine_toUpper =  prot_iter->getSearchEngine();
						searchengine_toUpper.toUpper();
						if(*engine == prot_iter->getSearchEngine() || *engine == searchengine_toUpper)
						{
							for(vector< PeptideIdentification >::iterator it = peptide_ids.begin();it < peptide_ids.end(); ++it)
							{
								if(prot_iter->getIdentifier().compare(it->getIdentifier()) == 0)
								{
									vector<PeptideHit> hits = it->getHits();
										for(std::vector<PeptideHit>::iterator  hit  = hits.begin(); hit < hits.end(); ++hit)
									{
										if(hit->getCharge() == *charge)
										{
											hit->setMetaValue("Search engine score",hit->getScore());
											hit->setScore(PEP_model.computeProbability(get_score_(*engine, *hit)));
										}
									}
									it->setHits(hits);	
								}
								it->setScoreType("Posterior Error Probability");
								it->setHigherScoreBetter(false);
							}
						}
					}				
					scores.clear();	
				}				
			}
		}
		else
		{
			for(StringList::iterator engine = search_engines.begin(); engine < search_engines.end(); ++engine)
			{
				for(vector< ProteinIdentification >::iterator prot_iter = protein_ids.begin(); prot_iter < protein_ids.end(); ++prot_iter)
				{
						String searchengine_toUpper =  prot_iter->getSearchEngine();
						searchengine_toUpper.toUpper();
						if(*engine == prot_iter->getSearchEngine() || *engine == searchengine_toUpper)
					{
						for(vector< PeptideIdentification >::iterator it = peptide_ids.begin();it < peptide_ids.end(); ++it)
						{
							if(prot_iter->getIdentifier().compare(it->getIdentifier()) == 0)
							{
								vector<PeptideHit> hits = it->getHits();
								for(std::vector<PeptideHit>::iterator  hit  = hits.begin(); hit < hits.end(); ++hit)
								{
									scores.push_back(get_score_(*engine, *hit));
								}
							}
						}
					}
				}
				PEP_model.fit(scores, probabilities);
				for(vector< ProteinIdentification >::iterator prot_iter = protein_ids.begin(); prot_iter < protein_ids.end(); ++prot_iter)
				{	
						String searchengine_toUpper =  prot_iter->getSearchEngine();
						searchengine_toUpper.toUpper();
						if(*engine == prot_iter->getSearchEngine() || *engine == searchengine_toUpper)
					{
						for(vector< PeptideIdentification >::iterator it = peptide_ids.begin();it < peptide_ids.end(); ++it)
						{
							if(prot_iter->getIdentifier().compare(it->getIdentifier()) == 0)
							{
								vector<PeptideHit> hits = it->getHits();
								for(std::vector<PeptideHit>::iterator  hit  = hits.begin(); hit < hits.end(); ++hit)
								{
									hit->setMetaValue("Search engine score",hit->getScore());
									hit->setScore(PEP_model.computeProbability(get_score_(*engine, *hit)));	
								}
								it->setHits(hits);
							}
							it->setScoreType("Posterior Error Probability");
							it->setHigherScoreBetter(false);
						}
					}
				}				
				scores.clear();	
			}
		}
		//-------------------------------------------------------------
		// writing output
		//-------------------------------------------------------------
	
		file.store(outputfile_name, protein_ids, peptide_ids);
		return EXECUTION_OK;
	}
	
	DoubleReal smallest_e_value_;
};


int main( int argc, const char** argv )
{
	TOPPIDPosteriorErrorProbability tool;

	return tool.main(argc,argv);
}

/// @endcond
