// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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

	@brief  Tool to estimate the probability of peptide hits to be incorrectly assigned.

	<CENTER>
	<table>
		<tr>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
			<td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ IDPosteriorErrorProbability \f$ \longrightarrow \f$</td>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_MascotAdapter (or other ID engines) </td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_ConsensusID </td>
		</tr>
	</table>
	</CENTER>

	@experimental This tool has not been tested thoroughly and might behave not as expected!

	By default an estimation is performed using the (inverse) Gumbel distribution for incorrectly assigned sequences
	and a Gaussian distribution for correctly assigned sequences. The probabilities are calculated by using Bayes' law, similar to PeptideProphet.
	Alternatively, a second Gaussian distribution can be used for incorrectly assigned sequences.
	At the moment, IDPosteriorErrorProbability is able to handle X!Tandem, Mascot, MyriMatch and OMSSA scores.
  
  No target/decoy information needs to be provided, since the model fits are done on the mixed distribution.

	In order to validate the computed probabilities one can adjust the fit_algorithm subsection.

	There are three parameters for the plot:
	The parameter 'output_plots' is by default false. If set to true the plot will be created.
	The scores are plotted in form of bins. Each bin represents a set of scores in a range of (highest_score - smallest_score)/number_of_bins (if all scores have positive values).
	The midpoint of the bin is the mean of the scores it represents.
	Finally, the parameter output_name should be used to give the plot a unique name. Two files are created. One with the binned scores and one with all steps of the estimation.
	If top_hits_only is set, only the top hits of each PeptideIndentification are used for the estimation process.
	Additionally, if 'top_hits_only' is set, target_decoy information are available and a False Discovery Rate run was performed before, an additional plot will be plotted with target and decoy bins(output_plot must be true in fit_algorithm subsection).
	A peptide hit is assumed to be a target if its q-value is smaller than fdr_for_targets_smaller.
	
	Actually, the plots are saved as a gnuplot file. Therefore, to visualize the plots one has to use gnuplot, e.g. gnuplot file_name. This should output a postscript file which contains all steps of the estimation.

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_IDPosteriorErrorProbability.cli
	<B>INI file documentation of this tool:</B>
	@htmlinclude TOPP_IDPosteriorErrorProbability.html

	For the parameters of the algorithm section see the algorithms documentation: @n
		@ref OpenMS::Math::PosteriorErrorProbabilityModel "fit_algorithm" @n

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPIDPosteriorErrorProbability
	: public TOPPBase
{
 public:
	TOPPIDPosteriorErrorProbability()
		: TOPPBase("IDPosteriorErrorProbability","Estimates probabilities for incorrectly assigned peptide sequences and a set of search engine scores using a mixture model.")
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
	  registerFlag_("split_charge", "The search engine scores are split by charge if this flag is set. Thus, for each charge state a new model will be computed.");
		registerFlag_("top_hits_only","If set only the top hits of every PeptideIdentification will be used");
		registerDoubleOption_("fdr_for_targets_smaller","<value>",0.05,"Only used, when top_hits_only set. Additionally, target_decoy information should be available. The score_type must be q-value from an previous False Discovery Rate run.",false,true);
	  registerFlag_("ignore_bad_data","If set errors will be written but ignored. Useful for pipelines with many datasets where only a few are bad, but the pipeline should run through.");
	  registerFlag_("prob_correct","If set scores will be calculated as 1-ErrorProbabilities and can be interpreted as probabilities for correct identifications.");
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
    else if( engine == "MyriMatch" )
    {
      //double e_val = exp(-hit.getScore());
      //double score_val = ((-1)* log10(max(e_val,smallest_e_value_)));
      //printf("myri score: %e ; e_val: %e ; score_val: %e\n",hit.getScore(),e_val,score_val);
      //return score_val;
      return hit.getScore();
    }
		else if(engine.compare("XTandem") == 0)
		{
			return((-1)* log10(max((DoubleReal)hit.getMetaValue("E-Value"),smallest_e_value_)));
		}
		else if (engine == "MASCOT")
		{
			if(hit.metaValueExists("EValue"))
			{
				return((-1)* log10(max((DoubleReal)hit.getMetaValue("EValue"),smallest_e_value_)));
			}
			if(hit.metaValueExists("expect"))
			{
				return((-1)* log10(max((DoubleReal)hit.getMetaValue("expect"),smallest_e_value_)));
			}
    } 
		else if (engine == "SpectraST")
    {
      return (100*hit.getScore());  // SpectraST f-val
    }
		else
		{
			throw Exception::UnableToFit(__FILE__,__LINE__,__PRETTY_FUNCTION__,"No parameters for chosen search engine","The chosen search engine is currently not supported");
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
		bool top_hits_only = getFlag_("top_hits_only");
		DoubleReal fdr_for_targets_smaller = getDoubleOption_("fdr_for_targets_smaller");
		bool target_decoy_available = false;
		bool ignore_bad_data = getFlag_("ignore_bad_data");
		bool prob_correct = getFlag_("prob_correct");
		//-------------------------------------------------------------
		// reading input
		//-------------------------------------------------------------
		bool return_value;
		IdXMLFile file;
		vector< ProteinIdentification > protein_ids;
		vector< PeptideIdentification > peptide_ids;
		file.load(inputfile_name, protein_ids, peptide_ids);
		vector<double> scores;
		vector<double> decoy;
		vector<double> target;
		vector<double> probabilities;
		vector<Int> charges;
		PosteriorErrorProbabilityModel PEP_model;
		PEP_model.setParameters(fit_algorithm);
    StringList search_engines = StringList::create("XTandem,OMSSA,MASCOT,SpectraST,MyriMatch");
		//-------------------------------------------------------------
		// calculations
		//-------------------------------------------------------------
		if(split_charge)
		{
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
			if(charges.empty())
			{
				throw Exception::ElementNotFound(__FILE__,	__LINE__,__PRETTY_FUNCTION__,"no charges found!");
			}
		}
		for(vector< PeptideIdentification >::iterator it = peptide_ids.begin();it < peptide_ids.end(); ++it)
		{
			if(!it->getHits().empty())
			{
				target_decoy_available = (it->getScoreType() == "q-value" &&  it->getHits()[0].getMetaValue("target_decoy") != DataValue::EMPTY);
				break;
			}
		}		
    
    vector<Int>::iterator charge = charges.begin(); // charges can be empty, no problem if split_charge is not set
    if (split_charge && charges.empty())
    {
      throw Exception::Precondition(__FILE__,__LINE__, __PRETTY_FUNCTION__, "split_charge is set and the list of charge states is empty but should not be!");
    }
		map<String,vector<vector<double> > > all_scores;
		char splitter = ',';//to split the engine from the charge state later on
		do
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
								if(top_hits_only)
								{
									if(!hits.empty() && (!split_charge || hits[0].getCharge() == *charge))
									{
										scores.push_back(get_score_(*engine, hits[0]));
										if(target_decoy_available)
										{
											if(hits[0].getScore() < fdr_for_targets_smaller)
											{
												target.push_back(get_score_(*engine, hits[0]));
										  }
								 			else
										 	{
												decoy.push_back(get_score_(*engine, hits[0]));
								 			}
										}
									}									
								}
								else
								{
									for(std::vector<PeptideHit>::iterator  hit  = hits.begin(); hit < hits.end(); ++hit)
									{
										if(!split_charge || hit->getCharge() == *charge)
										{
											scores.push_back(get_score_(*engine, *hit));
										}
									}
								}
							}
						}
					}
				}
				if(scores.size() > 2)
				{
					vector< vector <double > > tmp;
					tmp.push_back(scores);
					tmp.push_back(target);
					tmp.push_back(decoy);
					if(split_charge)
					{
						String engine_with_charge_state = *engine +String(splitter)+ String(*charge);
						all_scores.insert(make_pair(engine_with_charge_state,tmp ));
					}
					else
					{
						all_scores.insert(make_pair(*engine,tmp));
					}
				}
			
				scores.clear();
				target.clear();
				decoy.clear();
			}	

      if(split_charge) ++charge;

    } while(charge < charges.end());

		if(all_scores.empty() )
		{
			writeLog_("No data collected. Check whether search engine is supported.");
			if(!ignore_bad_data)return INPUT_FILE_EMPTY;
		}
		for(map<String, vector< vector<double> > >::iterator it =all_scores.begin(); it != all_scores.end();++it)
		{	
			vector<String> engine_info;
		  it->first.split(splitter, engine_info);
			String engine = engine_info[0];
			Int charge= -1;
			if(engine_info.size() == 2)
			{
				charge = engine_info[1].toInt();
			}
			if(split_charge) 
			{
				String output_name  = fit_algorithm.getValue("output_name");
				fit_algorithm.setValue("output_name", output_name + "_charge_" + String(charge) , "...",StringList::create("advanced,output file"));	
				PEP_model.setParameters(fit_algorithm);
			}
				
			return_value = PEP_model.fit(it->second[0]);
			if(!return_value )	 writeLog_("unable to fit data. Algorithm did not run through for the following search engine: " + engine);
			if(!return_value && !ignore_bad_data) return UNEXPECTED_RESULT;
				//plot target_decoy
			if(target_decoy_available && it->second[0].size() > 0 && return_value)
			{
			  PEP_model.plotTargetDecoyEstimation(it->second[1], it->second[2]);//target, decoy
			}					
			if(return_value)
			{
				bool unable_to_fit_data =true;
				bool data_might_not_be_well_fit = true;
				for(vector< ProteinIdentification >::iterator prot_iter = protein_ids.begin(); prot_iter < protein_ids.end(); ++prot_iter)
				{
					String searchengine_toUpper =  prot_iter->getSearchEngine();
					searchengine_toUpper.toUpper();
					
					if(engine == prot_iter->getSearchEngine() || engine == searchengine_toUpper)
					{
						for(vector< PeptideIdentification >::iterator it = peptide_ids.begin();it < peptide_ids.end(); ++it)
						{
							if(prot_iter->getIdentifier().compare(it->getIdentifier()) == 0)
							{
                String score_type = it->getScoreType() + "_score";
								vector<PeptideHit> hits = it->getHits();
									for(std::vector<PeptideHit>::iterator  hit  = hits.begin(); hit < hits.end(); ++hit)
								{
		 							if(!split_charge || hit->getCharge() == charge)
									{
                    DoubleReal score;
                    hit->setMetaValue(score_type, hit->getScore());
										score = PEP_model.computeProbability(get_score_(engine, *hit));
										if(score > 0 && score < 1) unable_to_fit_data = false; //only if all it->second[0] are 0 or 1 unable_to_fit_data stays true
										if(score > 0.2 && score < 0.8) data_might_not_be_well_fit = false;//same as above
										hit->setScore(score);
										if(prob_correct)
										{
											hit->setScore(1-score);
										}
										else
										{
											hit->setScore(score);
										}
									}
								}
								it->setHits(hits);
							}
							it->setScoreType("Posterior Error Probability");
							it->setHigherScoreBetter(false);
						}
					}
				}
				if(unable_to_fit_data) writeLog_(String("unable to fit data for search engine: ")+ engine);
				if(unable_to_fit_data && !ignore_bad_data) return UNEXPECTED_RESULT;
				if(data_might_not_be_well_fit) writeLog_(String("data might not be well fitted for search engine: ")+engine);
			}
 		}
		//-------------------------------------------------------------
		// writing output
		//-------------------------------------------------------------
		file.store(outputfile_name, protein_ids, peptide_ids);
		return EXECUTION_OK;
	}
	//Used in several functions
	DoubleReal smallest_e_value_;
};


int main( int argc, const char** argv )
{
	TOPPIDPosteriorErrorProbability tool;

	return tool.main(argc,argv);
}

/// @endcond
