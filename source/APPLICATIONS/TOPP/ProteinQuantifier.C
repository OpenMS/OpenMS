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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/SVOutStream.h>

#include <fstream>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_ProteinQuantifier ProteinQuantifier
	
	@brief Application to compute protein abundances from annotated feature/consensus maps

	Peptide/protein IDs from multiple identification runs can be handled, but will not be differentiated (i.e. protein accessions for a peptide will be accumulated over all identification runs).


	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_ProteinQuantifier.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

namespace OpenMS
{

  class TOPPProteinQuantifier 
  	: public TOPPBase
  {
	public:
		TOPPProteinQuantifier() :
			TOPPBase("ProteinQuantifier", "Compute protein abundances")
      {
      }

	protected:

		struct peptide_data
		{
			map<Int, DoubleReal> abundances; // abundance by charge
			set<String> accessions; // protein accessions
		};
		typedef map<AASequence, peptide_data> peptide_quant;
		

		bool isAnnotationAmbiguous_(Feature& feature)
			{
				// only the best peptide hit in each peptide ID is taken into account!
				vector<PeptideIdentification>& peptides = 
					feature.getPeptideIdentifications();
				if (peptides.empty()) return true;
				peptides.front().sort();
				const AASequence& seq = 
					peptides.front().getHits().front().getSequence();
				for (vector<PeptideIdentification>::iterator pep_it = 
							 ++peptides.begin(); pep_it != peptides.end(); ++pep_it)
				{
					pep_it->sort();
					if (pep_it->getHits().front().getSequence() != seq) return true;
				}
				return false;
			}


		void quantifyPeptides_(FeatureMap<>& features, peptide_quant& quant)
			{
				for (FeatureMap<>::iterator feat_it = features.begin(); 
						 feat_it != features.end(); ++feat_it)
				{
					// skip features without or with ambiguous annotation:
					if (isAnnotationAmbiguous_(*feat_it)) continue;

					vector<PeptideIdentification>& peptides = 
						feat_it->getPeptideIdentifications();
					const AASequence& seq = 
						peptides.front().getHits().front().getSequence();
					quant[seq].abundances[feat_it->getCharge()] += 
						feat_it->getIntensity(); // new map element is initialized with 0
					for (vector<PeptideIdentification>::iterator pep_it = 
								 peptides.begin(); pep_it != peptides.end(); ++pep_it)
					{
						// check charge of peptide hit against charge of feature?
						quant[seq].accessions.insert(
							pep_it->getHits().front().getProteinAccessions().begin(),
							pep_it->getHits().front().getProteinAccessions().end());
					}
				}
			}

		
		void writePeptideTable_(peptide_quant& quant, SVOutStream& out)
			{
				out << "peptide" << "charge" << "abundance" << "protein"
						<< "n_diff_protein" << endl;
				bool filter_charge = getFlag_("filter_charge");
				for (peptide_quant::iterator q_it = quant.begin(); q_it != quant.end();
						 ++q_it)
				{
					StringList accessions;
					for (set<String>::iterator acc_it = q_it->second.accessions.begin();
							 acc_it != q_it->second.accessions.end(); ++acc_it)
					{
						String acc = *acc_it;
						accessions << acc.substitute('/', '_');
					}
					String protein = accessions.concatenate("/");
					if (!filter_charge)
					{ // replace individual abundances by sum (stored in charge 0):
						DoubleReal sum = 0;
						for (map<Int, DoubleReal>::iterator ab_it = 
									 q_it->second.abundances.begin(); ab_it != 
									 q_it->second.abundances.end(); ++ab_it)
						{
							sum += ab_it->second;
						}
						q_it->second.abundances.clear();
						q_it->second.abundances[0] = sum;
					}
					for (map<Int, DoubleReal>::iterator ab_it = 
								 q_it->second.abundances.begin(); ab_it != 
								 q_it->second.abundances.end(); ++ab_it)
					{
						out << q_it->first.toString() << ab_it->first << ab_it->second
								<< protein << accessions.size() << endl;
					}
				}
			}

	
		void registerOptionsAndFlags_()
      {
				registerInputFile_("in", "<file>", "", "Input file");
				setValidFormats_("in", StringList::create("featureXML,consensusXML"));
        registerOutputFile_("out", "<file>", "", "Output file for protein abundances (CSV format)");
				registerOutputFile_("peptide_out", "<file>", "", "Output file for peptide abundances (CSV format)", false);
				registerIntOption_("top", "<number>", 5, "Calculate protein abundance from this number of proteotypic peptides (best first; '0' for all)", false);
				setMinInt_("top", 0);
				registerFlag_("include_fewer", "Include results for proteins with fewer than 'top' proteotypic peptides");
				registerFlag_("filter_charge", "Set this flag to distinguish between charge states of the same peptide. For peptides, abundances will be reported separately for each charge;\nfor proteins, abundances will be computed based only on the most prevalent charge of each peptide.\n(By default, abundances are summed over all charge states.)");
				registerStringOption_("average", "<method>", "median", "Averaging method for computing protein abundances from peptide abundances", false);
				setValidStrings_("average", StringList::create("median,mean"));
				addEmptyLine_();
        addText_("Options for consensusXML files:");
				registerFlag_("normalize", "Scale peptide abundances so that medians of all samples are equal");
      }

		
		ExitCodes main_(int, const char**)
      {
				String in = getStringOption_("in");
				FileTypes::Type in_type = FileHandler::getType(in);

				String peptide_out = getStringOption_("peptide_out");
				
				if (in_type == FileTypes::FEATUREXML)
				{
					FeatureMap<> features;
          FeatureXMLFile().load(in, features);
					peptide_quant quant;
					quantifyPeptides_(features, quant);

					if (!peptide_out.empty())
					{
						ofstream outstr(peptide_out.c_str());
						SVOutStream output(outstr, ",", "_", String::DOUBLE);
						writePeptideTable_(quant, output);
						outstr.close();
					}
				}

				else // consensusXML
				{
					throw Exception::NotImplemented(__FILE__, __LINE__, 
																					__PRETTY_FUNCTION__);
				}

				return EXECUTION_OK;
			}
	};
	
} 


int main(int argc, const char** argv)
{
  TOPPProteinQuantifier t;
  return t.main(argc, argv);
}

/// @endcond
