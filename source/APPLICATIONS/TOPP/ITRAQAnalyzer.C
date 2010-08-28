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
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqChannelExtractor.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqQuantifier.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/SYSTEM/File.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_ITRAQAnalyzer ITRAQAnalyzer

	@brief Extracts and normalizes iTRAQ information from an MS experiment.

<CENTER>
	<table>
		<tr>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
			<td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ ITRAQAnalyzer \f$ \longrightarrow \f$</td>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FileConverter </td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_IDMapper</td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FileFilter </td>
		</tr>
	</table>
</CENTER>

	@experimental This tool has not been tested thoroughly and might behave not as expected!

	Provide an idXML file that you obtained from the same data (e.g. by using InspectAdapter)
	to have protein ratios reported, instead of peptide ratios.

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_ITRAQAnalyzer.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPITRAQAnalyzer
	: public TOPPBase
{
 public:
	TOPPITRAQAnalyzer()
		: TOPPBase("ITRAQAnalyzer","\n\nWARNING: EXPERIMENTAL\n\n Calculates iTRAQ quantitative values for peptides or proteins (when idXML available)", true, true)
	{
	}

 protected:
	void registerOptionsAndFlags_()
	{
		registerStringOption_("type","<name>","","iTRAQ experiment type\n",true);
		setValidStrings_("type", getToolList()[toolName_()] );

		registerInputFile_("in","<file>","","input raw/picked data file ");
		setValidFormats_("in",StringList::create("mzML"));
		registerOutputFile_("out","<file>","","output consensusXML file with quantitative information");
		setValidFormats_("out",StringList::create("consensusXML"));

		registerInputFile_ ("idxml", "<file>", "", "!not supported yet! idXML file with peptide identifications from tandemMS of the -in file", false, false);

		addEmptyLine_();
		//addText_("Note: We highly recommend providing an idXML file with identifications. This enables ITRAQAnalyzer to report protein ratios!");

  	registerSubsection_("algorithm","Algorithm parameters section");

			// report ProteinIDs for Peptides: Mascot, OpenSource: XTandem, OMSA
			//--> filter for search engine!
			// to-check: SEQUEST?
			//or check for Protein-Candidates manually! suffix-array
	}

	Param getSubsectionDefaults_(const String& /*section*/) const
	{
		String type = getStringOption_("type");
		Int t = (type=="4plex" ?  ItraqQuantifier::FOURPLEX : ItraqQuantifier::EIGHTPLEX );
	  Param tmp;
		tmp.insert("Extraction:",ItraqChannelExtractor(t).getParameters());
	  tmp.insert("Quantification:",ItraqQuantifier(t).getParameters());
		tmp.setValue ("MetaInformation:Program", "OpenMS::ITRAQAnalyzer", "", StringList::create("advanced"));
	  return tmp;
	}

 	ExitCodes main_(int , const char**)
	{
		//-------------------------------------------------------------
		// parameter handling
		//-------------------------------------------------------------
		String in = getStringOption_("in");
		String out = getStringOption_("out");
		String idxml = getStringOption_("idxml");

		Int itraq_type = (getStringOption_("type")=="4plex" ?  ItraqQuantifier::FOURPLEX : ItraqQuantifier::EIGHTPLEX );
		//-------------------------------------------------------------
		// loading input
		//-------------------------------------------------------------

		MzMLFile mz_data_file;
		MSExperiment<Peak1D > exp;
		mz_data_file.setLogType(log_type_);
		mz_data_file.load(in,exp);

		//-------------------------------------------------------------
		// calculations
		//-------------------------------------------------------------
		Param extract_param(getParam_().copy("algorithm:Extraction:",true));
		ItraqChannelExtractor itraq_ce(itraq_type, extract_param);

		ConsensusMap consensus_map_raw, consensus_map_quant;
		// extract raw signals
		itraq_ce.run(exp, consensus_map_raw);

		// do normalization
		Param quant_param(getParam_().copy("algorithm:Quantification:",true));
		ItraqQuantifier itraq_quant(itraq_type, quant_param);

		// TODO: allow multiple input files and do statistics on labelling efficiency etc...

		if (File::readable(idxml))
		{
			IdXMLFile f;
			std::vector< ProteinIdentification > protein_ids;
			std::vector< PeptideIdentification > peptide_ids;
			String document_id;
			f.load (idxml, protein_ids, peptide_ids, document_id);
			itraq_quant.run(consensus_map_raw, peptide_ids, protein_ids, consensus_map_quant);
		}
		else
		{
			itraq_quant.run(consensus_map_raw, consensus_map_quant);
		}


		// assign unique ID to output file (this might throw an exception.. but thats ok, as we want the programm to quit then)
		if (getStringOption_("id_pool").trim().length()>0) getDocumentIDTagger_().tag(consensus_map_raw);

		// annotate output file with MetaInformation
		Param metainfo_param(getParam_().copy("algorithm:MetaInformation:",true));
		for (Param::ParamIterator it = metainfo_param.begin(); it!=metainfo_param.end(); ++it)
		{
			consensus_map_raw.setMetaValue (it->name, it->value);
		}


		//-------------------------------------------------------------
		// writing output
		//-------------------------------------------------------------
		ConsensusXMLFile cm_file;

		//annotate output with data processing info
		addDataProcessing_(consensus_map_raw, getProcessingInfo_(DataProcessing::QUANTITATION));
		addDataProcessing_(consensus_map_quant, getProcessingInfo_(DataProcessing::QUANTITATION));

		cm_file.store(out,consensus_map_raw);
		cm_file.store(out+"_quant",consensus_map_quant);

		return EXECUTION_OK;
	}

};

int main( int argc, const char** argv )
{
    TOPPITRAQAnalyzer tool;
    return tool.main(argc,argv);
}

/// @endcond
