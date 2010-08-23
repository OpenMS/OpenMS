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
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/ANALYSIS/DECHARGING/FeatureDeconvolution.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
   @page TOPP_Decharger Decharger

   @brief Decharges a feature map by clustering charge variants of a peptide to zero-charge entities.
<CENTER>
	<table>
		<tr>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
			<td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ Decharger \f$ \longrightarrow \f$</td>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureFinder </td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_ProteinQuantifier</td>
		</tr>
	</table>
</CENTER>

   The Decharger uses an ILP approach to group charge variants of the same peptide, which
   usually occur in ESI ionization mode. The resulting zero-charge peptides, which are defined by RT and mass,
   are written to consensusXML. Intensities of charge variants are summed up. The position of the zero charge
   variant is the average of all clustered peptides in each dimension (m/z and RT).
   It is also possible to include adducted species to the charge ladders (see 'potential_adducts' parameter).
   Via this mechanism it is also possible to use this tool to find pairs/triples/quadruples/... in labeled data (by specifing the mass
   tag weight as an adduct). If mass tags induce an RT shift (e.g. deuterium labeled data) you can also specify this also in the adduct list.
   This will allow to tighten the RT search window, thus reducing false positive results.

	 <B>The command line parameters of this tool are:</B>
   @verbinclude TOPP_Decharger.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPDecharger
  : virtual public TOPPBase
{
 public:
  TOPPDecharger()
    : TOPPBase("Decharger","Decharges and merges different feature charge variants of the same peptide.")
  {
  }

 protected:
  void registerOptionsAndFlags_()
  {
    registerInputFile_("in","<file>","","input file ");
		setValidFormats_("in",StringList::create("FeatureXML"));
    registerOutputFile_("out_fm","<file>","","output feature map");
    registerOutputFile_("out_cm","<file>","","output consensus map");
    registerOutputFile_("outpairs","<file>","","output file");
	  setValidFormats_("out_fm",StringList::create("FeatureXML"));
	  setValidFormats_("out_cm",StringList::create("ConsensusXML"));
	  setValidFormats_("outpairs",StringList::create("ConsensusXML"));

    addEmptyLine_();
    addText_("All other options of the Decharger depend on the FeatureDeconvolution class.\n"
             "They can be given only in the 'algorithm' section  of the INI file.");
    
    registerSubsection_("algorithm","Feature decharging algorithm section");
  }

	Param getSubsectionDefaults_(const String& /*section*/) const
	{
	  // there is only one subsection: 'algorithm' (s.a) .. and in it belongs the FeatureDecharger param
	  FeatureDeconvolution fdc;
	  Param tmp;
	  tmp.insert("FeatureDeconvolution:",fdc.getParameters());
	  return tmp;
	}

  ExitCodes main_(int , const char**)
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------
    String infile = getStringOption_("in");
    String outfile_fm = getStringOption_("out_fm");
    String outfile_cm = getStringOption_("out_cm");
    String outfile_p = getStringOption_("outpairs");

    FeatureDeconvolution fdc;
    Param const& dc_param = getParam_().copy("algorithm:FeatureDeconvolution:",true);

    writeDebug_("Parameters passed to Decharger", dc_param, 3);
    
    fdc.setParameters(dc_param);

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------

    writeDebug_("Loading input file", 1);
    
    typedef FeatureMap<> FeatureMapType;
    FeatureMapType map_in, map_out;
    FeatureXMLFile().load(infile, map_in);

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    ConsensusMap cm, cm2;
    fdc.compute(map_in, map_out, cm, cm2);
    
    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    
    writeDebug_("Saving output files", 1);

    cm.getFileDescriptions()[0].filename = infile;
    cm2.getFileDescriptions()[0].filename = infile;

		//annotate output with data processing info
		addDataProcessing_(map_out, getProcessingInfo_(DataProcessing::CHARGE_DECONVOLUTION));
		addDataProcessing_(cm, getProcessingInfo_(DataProcessing::CHARGE_DECONVOLUTION));
		addDataProcessing_(cm2, getProcessingInfo_(DataProcessing::CHARGE_DECONVOLUTION));

		FeatureXMLFile().store(outfile_fm, map_out);
    ConsensusXMLFile f;
    f.store(outfile_cm, cm);
    f.store(outfile_p, cm2);
    

    return EXECUTION_OK;
  }
};


int main( int argc, const char** argv )
{
    TOPPDecharger tool;
    return tool.main(argc,argv);
}

/// @endcond
