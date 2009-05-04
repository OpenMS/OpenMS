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
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/DECHARGING/FeatureDecharger.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/ANALYSIS/DECHARGING/FeatureDeconvolution.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <coin/OsiClpSolverInterface.hpp>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
   @page TOPP_Decharger Decharger

   @brief Decharges a feature map by clustering charge variants of a peptide to zero-charge entities.

   The Decharger uses a hierarchical clustering (complete linkage) to group charge variants of the same peptide, which
   usually occur in ESI ionization mode. The resulting zero-charge peptides, which are defined by RT and mass,
   are written to a featureXML file. Intensities of charge variants are summed up. The position of the zero charge
   variant is the average of all clustered peptides in each dimension.
   If several peptides with the same charge variant are grouped (which is clearly not allowed), a heuristic is used:
   <ul>
   <li>cluster consists of only one charge variant (but several peptides) -> split cluster into single elements</li>
   <li>cluster consists of several charge variants -> dispose cluster</li>
   </ul>
   

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_DBImporter.cli
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
    registerOutputFile_("out","<file>","","output file");
    registerOutputFile_("outpairs","<file>","","output file");
	  setValidFormats_("out",StringList::create("ConsensusXML"));
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
    String outfile = getStringOption_("out");
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
    FeatureMapType map;
    FeatureXMLFile().load(infile, map);
    //map.sortByPosition();

		OsiClpSolverInterface solver;
    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    ConsensusMap cm, cm2;
    fdc.compute(map, cm,cm2);
    
    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    
    writeDebug_("Saving output files", 1);

    cm.getFileDescriptions()[0].filename = infile;
    cm2.getFileDescriptions()[0].filename = infile;

    ConsensusXMLFile f;
    f.store(outfile, cm);
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
