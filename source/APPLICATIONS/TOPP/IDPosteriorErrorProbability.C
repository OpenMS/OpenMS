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
	
	@brief  Tool to estimate the probability of peptide hit to be incorrectly assigned.
	
	By default an estimation is performed using the (inverse) gumbel distribution for incorrectly assigned sequences 
	and a gaussian distribution for correctly assigned sequences. The probabilities are calculated using bayes law, similar to PeptideProphet.
	Alternatively, a second gaussian distribution can be used for incorreclty assigned sequences.
	At the moment, it is able to handle Xtandem, Mascot and OMSSA scores.
	
	In order to validate the computed probabilities one can adjust the fit_algorithm subsection. 
	The easiest way to do this, is to create a default ini file with the parameter -write_ini <file_name>.
	Secondly, it is suggested to open the created ini-file with the INIFileEditor. 
	There are three parameters for the plot: 
	The parameter output_plots is by default false. If set to true the plots will be created.
	The scores are plotted in form of bins. Each bin represents a set of scores in a range of (highest_score - smallest_score)/number_of_bins (if all scores have positive values).
	The middle point of the bin is the mean of the scores it represents.
	Finally, the parameter output_name should be used to give the plot a unique name. The algorithm then creates two files. One with the binned scores and one with all steps of the estimation.
	
	The plots are saved in a gnuplot file. Therefore, to visualize the plots as a pdf one has to use gnuplot. e.g. gnuplot <filename>. This should output a file which contains all steps of the estimation.
	
	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_IDPosteriorErrorProbability.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPIDPosteriorErrorProbability
	: public TOPPBase
{
 public:
	TOPPIDPosteriorErrorProbability()
		: TOPPBase("IDPosteriorErrorProbability","Estimates peptide probabilities for a set of search engine scores using a mixture model.")
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
