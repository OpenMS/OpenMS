// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqChannelExtractor.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqQuantifier.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
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

  Extract the iTRAQ reporter ion intensities (4plex or 8plex) from
  raw MS2 data, does isotope corrections and stores the resulting quantitation as 
  consensusXML, where each consensus centroid corresponds to one iTRAQ MS2 scan (e.g., HCD).
  The position of the centroid is the precursor position, its sub-elements are the channels (thus having
  m/z's of 113-121).

  Isotope correction is done using non-negative least squares (NNLS), i.e.,

  Minimize ||Ax - b||, subject to x >= 0, where b is the vector of observed reporter intensities (with 'contaminating'
  isotope species), A is a correction matrix (as supplied by the manufacturer AB Sciex) and x is the desired vector of corrected (real)
  reporter intensities.
  Other software solves this problem using an inverse matrix multiplication, but this can yield entries in x which are negative. In a real sample,
  this solution cannot possibly be true, so usually negative values (= negative reporter intensities) are set to 0.
  However, a negative result usually means, that noise was not accounted for thus we use NNLS to get a non-negative solution, without the need to
  truncate negative values. In (the usual) case that inverse matrix multiplication yields only positive values, our NNLS will give 
  the exact same optimal solution.

  The correction matrices can be found (and changed) in the INI file. However, these matrices for both 4plex and 8plex are now stable, and every
  kit delivered should have the same isotope correction values. Thus, there should be no need to change them, but feel free to compare the values in the
  INI file with your kit's Certificate.
  
  After this quantitation step, you might want to annotate the consensus elements with the respective identifications, obtained from
  an identification pipeline.

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
		: TOPPBase("ITRAQAnalyzer","Calculates iTRAQ quantitative values for peptides", true, true)
	{
	}

 protected:
	void registerOptionsAndFlags_()
	{
		registerStringOption_("type", "<mode>", "4plex", "iTRAQ experiment type\n", false);
		setValidStrings_("type", StringList::create("4plex,8plex"));

		registerInputFile_("in", "<file>", "", "input raw/picked data file ");
		setValidFormats_("in", StringList::create("mzML"));
		registerOutputFile_("out", "<file>", "", "output consensusXML file with quantitative information");
		setValidFormats_("out", StringList::create("consensusXML"));

    registerOutputFile_("out_stats", "<file>", "", "output statistics as tab-separated file (readable by R or Excel or ...)", false);
    setValidFormats_("out_stats", StringList::create("tsv"));

    addEmptyLine_();

    registerSubsection_("algorithm","Algorithm parameters section");
	}

	Param getSubsectionDefaults_(const String& /*section*/) const
	{
	  Param tmp;
		tmp.insert("Extraction:", ItraqChannelExtractor(ItraqQuantifier::FOURPLEX).getParameters()); // type is irrelevant - ini is the same
	  tmp.insert("Quantification:", ItraqQuantifier(ItraqQuantifier::FOURPLEX).getParameters());   // type is irrelevant - ini is the same
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
    String out_stats = getStringOption_("out_stats");

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

		itraq_quant.run(consensus_map_raw, consensus_map_quant);

		// assign unique ID to output file (this might throw an exception.. but thats ok, as we want the program to quit then)
		if (getStringOption_("id_pool").trim().length()>0) getDocumentIDTagger_().tag(consensus_map_quant);

		// annotate output file with MetaInformation
		Param metainfo_param(getParam_().copy("algorithm:MetaInformation:",true));
		for (Param::ParamIterator it = metainfo_param.begin(); it!=metainfo_param.end(); ++it)
		{
			consensus_map_quant.setMetaValue (it->name, it->value);
		}


		//-------------------------------------------------------------
		// writing output
		//-------------------------------------------------------------

		//annotate output with data processing info
		addDataProcessing_(consensus_map_quant, getProcessingInfo_(DataProcessing::QUANTITATION));

    // add filename references
    for (ConsensusMap::FileDescriptions::iterator it = consensus_map_quant.getFileDescriptions().begin();
                                                  it != consensus_map_quant.getFileDescriptions().end();
                                                  ++it)
    {
      it->second.filename = in;
    }

		ConsensusXMLFile cm_file;
		cm_file.store(out, consensus_map_quant);
    
    std::cout << itraq_quant.getStats();
    if (!out_stats.trim().empty())
    {
      ofstream f;
      f.open(out_stats.c_str(), ios_base::out);
      f << itraq_quant.getStats();
      f.close();
    }

		return EXECUTION_OK;
	}

};

int main( int argc, const char** argv )
{
    TOPPITRAQAnalyzer tool;
    return tool.main(argc,argv);
}

/// @endcond
