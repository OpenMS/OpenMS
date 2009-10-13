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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/ANALYSIS/ID/ConsensusID.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_ConsensusID ConsensusID

	@brief Computes a consensus identification from peptide identification engines.

	The input file can contain several searches, e.g. from several identification engines.
	You can combine several searches with @ref TOPP_IDMerger. Identification runs can be mapped
	to featureXML and consensusXML with the @ref TOPP_IDMapper tool.

	For a detailed description of the algorithms and parameters see the documentation of
	the %OpenMS ConsensusID class.

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_ConsensusID.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

//Helper class
struct IDData
{
	DoubleReal mz;
	DoubleReal rt;
	vector<PeptideIdentification> ids;
};

class TOPPConsensusID
	: public TOPPBase
{
	public:
		TOPPConsensusID()
			: TOPPBase("ConsensusID","Computes a consensus identification from peptide identifications of several identification engines.")
		{
		}

	protected:

		Param getSubsectionDefaults_(const String& /*section*/) const
		{
			return ConsensusID().getDefaults();
		}

		void registerOptionsAndFlags_()
		{
			registerInputFile_("in","<file>","","input file");
			setValidFormats_("in",StringList::create("idXML,featureXML,consensusXML"));
			registerOutputFile_("out","<file>","","output file");
			setValidFormats_("out",StringList::create("idXML,featureXML,consensusXML"));

			addEmptyLine_();
			registerDoubleOption_("rt_delta","<value>",0.1, "Maximum allowed precursor RT deviation between identifications.", false);
			setMinFloat_("rt_delta",0.0);
			registerDoubleOption_("mz_delta","<value>",0.1, "Maximum allowed precursor m/z deviation between identifications.", false);
			setMinFloat_("mz_delta",0.0);

			registerSubsection_("algorithm","Consensus algorithm section");
		}

		ExitCodes main_(int , const char**)
		{
			String in = getStringOption_("in");
			FileTypes::Type in_type = FileHandler::getType(in);
			String out = getStringOption_("out");

			DoubleReal rt_delta = getDoubleOption_("rt_delta");
			DoubleReal mz_delta = getDoubleOption_("mz_delta");

			//----------------------------------------------------------------
			//set up ConsensusID
			//----------------------------------------------------------------
			ConsensusID consensus;
			Param alg_param = getParam_().copy("algorithm:",true);
			if (alg_param.empty())
			{
				writeLog_("No parameters for ConsensusID given. Aborting!");
				return ILLEGAL_PARAMETERS;
			}
			writeDebug_("Parameters passed to ConsensusID (without number of runs)", alg_param, 3);

			//----------------------------------------------------------------
			// idXML
			//----------------------------------------------------------------
			if (in_type == FileTypes::IDXML)
			{
				vector<ProteinIdentification> prot_ids;
				vector<PeptideIdentification> pep_ids;
				String document_id;
				IdXMLFile().load(in,prot_ids, pep_ids, document_id);

				//merge peptide ids by precursor position
				vector<IDData> prec_data;
				for (vector<PeptideIdentification>::iterator pep_id_it = pep_ids.begin(); pep_id_it != pep_ids.end(); ++pep_id_it)
				{
					DoubleReal rt = (DoubleReal)(pep_id_it->getMetaValue("RT"));
					DoubleReal mz = (DoubleReal)(pep_id_it->getMetaValue("MZ"));
					writeDebug_(String("  ID: ") + rt + " / " + mz, 4);
					vector<IDData>::iterator pos = prec_data.begin();
					while (pos != prec_data.end())
					{
						if (fabs(pos->rt - rt) < rt_delta && fabs(pos->mz - mz) < mz_delta )
						{
							break;
						}
						++pos;
					}
					//right position was found => append ids
					if (pos != prec_data.end())
					{
						writeDebug_(String("    Appending IDs to precursor: ") + pos->rt + " / " + pos->mz, 4);
						pos->ids.push_back(*pep_id_it);
					}
					//insert new entry
					else
					{
						IDData tmp;
						tmp.mz = mz;
						tmp.rt = rt;
						tmp.ids.push_back(*pep_id_it);
						prec_data.push_back(tmp);
						writeDebug_(String("    Inserting new precursor: ") + tmp.rt + " / " + tmp.mz, 4);
					}
				}

				//compute consensus
				alg_param.setValue("number_of_runs",(UInt)prot_ids.size());
				consensus.setParameters(alg_param);
				for (vector<IDData>::iterator it = prec_data.begin(); it!=prec_data.end(); ++it)
				{
					writeDebug_(String("Calculating consensus for : ") + it->rt + " / " + it->mz + " #peptide ids: " + it->ids.size(), 4);
					consensus.apply(it->ids);
				}

				// writing output
				pep_ids.clear();
				for (vector<IDData>::iterator it = prec_data.begin(); it!=prec_data.end(); ++it)
				{
					pep_ids.push_back(it->ids[0]);
					pep_ids.back().setMetaValue("RT",it->rt);
					pep_ids.back().setMetaValue("MZ",it->mz);
				}

				//create new identification run
				vector< ProteinIdentification > prot_id_out(1);
				prot_id_out[0].setDateTime(DateTime::now());
				prot_id_out[0].setSearchEngine("OpenMS/ConsensusID");
				prot_id_out[0].setSearchEngineVersion(VersionInfo::getVersion());

				//store consensus
				IdXMLFile().store(out,prot_id_out,pep_ids);
			}

			//----------------------------------------------------------------
			// featureXML
			//----------------------------------------------------------------
			if (in_type == FileTypes::FEATUREXML)
			{
				//load map
				FeatureMap<> map;
				FeatureXMLFile().load(in,map);

				//compute consensus
				alg_param.setValue("number_of_runs",(UInt)map.getProteinIdentifications().size());
				consensus.setParameters(alg_param);
				for (Size i = 0; i < map.size(); ++i)
				{
					consensus.apply(map[i].getPeptideIdentifications());
				}

				//create new identification run
				map.getProteinIdentifications().clear();
				map.getProteinIdentifications().resize(1);
				map.getProteinIdentifications()[0].setDateTime(DateTime::now());
				map.getProteinIdentifications()[0].setSearchEngine("OpenMS/ConsensusID");
				map.getProteinIdentifications()[0].setSearchEngineVersion(VersionInfo::getVersion());

				//store consensus
				FeatureXMLFile().store(out,map);
			}

			//----------------------------------------------------------------
			// consensusXML
			//----------------------------------------------------------------
			if (in_type == FileTypes::CONSENSUSXML)
			{
				//load map
				ConsensusMap map;
				ConsensusXMLFile().load(in,map);

				//compute consensus
				alg_param.setValue("number_of_runs",(UInt)map.getProteinIdentifications().size());
				consensus.setParameters(alg_param);
				for (Size i = 0; i < map.size(); ++i)
				{
					consensus.apply(map[i].getPeptideIdentifications());
				}

				//create new identification run
				map.getProteinIdentifications().clear();
				map.getProteinIdentifications().resize(1);
				map.getProteinIdentifications()[0].setDateTime(DateTime::now());
				map.getProteinIdentifications()[0].setSearchEngine("OpenMS/ConsensusID");
				map.getProteinIdentifications()[0].setSearchEngineVersion(VersionInfo::getVersion());

				//store consensus
				ConsensusXMLFile().store(out,map);
			}

			return EXECUTION_OK;
		}
};


int main( int argc, const char** argv )
{
	TOPPConsensusID tool;
	return tool.main(argc,argv);
}

/// @endcond
