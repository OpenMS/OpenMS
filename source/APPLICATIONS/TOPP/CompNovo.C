// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------


#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/ANALYSIS/DENOVO/CompNovoIdentification.h>
#include <OpenMS/ANALYSIS/DENOVO/CompNovoIdentificationCID.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>

using namespace OpenMS;
using namespace std;

/**
	@page TOPP_CompNovo CompNovo
	
	@brief Performs a peptide/protein identification with the CompNovo engine.

	The CompNovo engine can operate in CID/ETD mode were the ETD spectra of the
	same precursor are used to support the de novo identification. In CID mode
	all spectra are assumed to be CID spectra only.

	The details are described in the publication:
 
	Andreas Bertsch, Andreas Leinenbach, Anton Pervukhin, Markus Lubeck, 
	Ralf Hartmer, 	Carsten Baessmann, Yasser A Elnakady, Rolf M&uuml;ller, 
	Sebastian B&ouml;cker, Christian G Huber and Oliver Kohlbacher  (2009)
	"De novo peptide sequencing by tandem MS using complementary CID and 
	electron transfer dissociation"
	Electrophoresis, 30(21):3736-3747. (PubMed ID: 19862751)

	@experimental This implementation may contain bugs!
	
	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_CompNovo.cli
*/


// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPCompNovo
	: public TOPPBase
{
	public:
		TOPPCompNovo()
			: TOPPBase("CompNovo", "Performs a de novo peptide identification using the CompNovo engine.")
		{
		}
	
	protected:

		Param getSubsectionDefaults_(const String& /*section*/) const
    {
			String type = getStringOption_("type");
			if (type == "CompNovo")
			{
      	return CompNovoIdentification().getDefaults();
			}
			if (type == "CompNovoCID")
			{
				return CompNovoIdentificationCID().getDefaults();
			}
			return CompNovoIdentification().getDefaults();
    }
		
		void registerOptionsAndFlags_()
		{
			registerInputFile_("in", "<file>", "", "input file in mzML format", true);
			setValidFormats_("in", StringList::create("mzML"));

			registerOutputFile_("out", "<file>", "", "output file in IdXML format", true);
			setValidFormats_("out", StringList::create("idXML"));

			registerStringOption_("type","<name>","","type of algorithm",true);
			setValidStrings_("type", StringList::create("CompNovo,CompNovoCID"));

			registerSubsection_("algorithm","Algorithm section");
			addEmptyLine_();
		}
		
		ExitCodes main_(int , const char**)
		{
			//-------------------------------------------------------------
			// parameter handling
			//-------------------------------------------------------------
	
			//input/output files
			String in(getStringOption_("in"));
			String out(getStringOption_("out"));

			String type(getStringOption_("type"));

      //-------------------------------------------------------------
      // loading input
      //-------------------------------------------------------------

      PeakMap exp;
			MzMLFile f;
      f.setLogType(log_type_);

	    PeakFileOptions options;
      options.clearMSLevels();
			options.addMSLevel(2);
      f.getOptions() = options;
      f.load(in, exp);

			writeDebug_("Data set contains " + String(exp.size()) + " spectra", 1);
			
      //-------------------------------------------------------------
      // calculations
      //-------------------------------------------------------------
		
			vector<PeptideIdentification> pep_ids;
			Param type_param = getParam_().copy("algorithm:",true);
		
			if (type == "CompNovo")
			{
		  	CompNovoIdentification comp_novo_id;

				// set the options
				comp_novo_id.setParameters(type_param);
				comp_novo_id.getIdentifications(pep_ids, exp);
			}
			else
			{
				if (type == "CompNovoCID")
				{
					CompNovoIdentificationCID comp_novo_id;
					comp_novo_id.setParameters(type_param);
					comp_novo_id.getIdentifications(pep_ids, exp);
				}
			}

			//-------------------------------------------------------------
			// writing output
			//-------------------------------------------------------------

  		DateTime now = DateTime::now();
  		String date_string = now.get();
  		String identifier(type + "_" + date_string);

  		for (vector<PeptideIdentification>::iterator it = pep_ids.begin(); it != pep_ids.end(); ++it)
  		{
    		it->assignRanks();
    		it->setIdentifier(identifier);
 			}

		  vector<ProteinIdentification> prot_ids;
		  ProteinIdentification prot_id;
		  prot_id.setIdentifier(identifier);
			prot_id.setDateTime(now);

			ProteinIdentification::SearchParameters search_parameters;
			search_parameters.charges = "+2-+3";
			if (type_param.getValue("tryptic_only").toBool())
			{
				search_parameters.enzyme = ProteinIdentification::TRYPSIN;
			}
			else
			{
				search_parameters.enzyme = ProteinIdentification::NO_ENZYME;
			}
      search_parameters.mass_type = ProteinIdentification::MONOISOTOPIC;
      search_parameters.fixed_modifications = (StringList)type_param.getValue("fixed_modifications");
      search_parameters.variable_modifications = (StringList)type_param.getValue("variable_modifications");

      search_parameters.missed_cleavages = (UInt)type_param.getValue("missed_cleavages");
      search_parameters.peak_mass_tolerance = (DoubleReal)type_param.getValue("fragment_mass_tolerance");
      search_parameters.precursor_tolerance = (DoubleReal)type_param.getValue("precursor_mass_tolerance");
			prot_id.setSearchParameters(search_parameters);
			prot_id.setSearchEngineVersion("0.9beta");
      prot_id.setSearchEngine("CompNovo");
			prot_ids.push_back(prot_id);

  		IdXMLFile().store(out, prot_ids, pep_ids);
	
			return EXECUTION_OK;
		}
};

/// @endcond


int main( int argc, const char** argv )
{
	TOPPCompNovo tool;
	return tool.main(argc,argv);
}

