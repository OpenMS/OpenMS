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
// $Maintainer: Chris Bielow, Andreas Bertsch $
// $Authors: Chris Bielow, Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FILTERING/TRANSFORMERS/SpectraMerger.h>

#include <algorithm>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_SpectraMerger SpectraMerger
	
	@brief Allows to add up several spectra.

	@experimental This TOPP-tool is not well tested and not all features might be properly implemented and tested!

  This tool can add several consecutive scans, increasing S/N ratio (for MS1 and above)<br>
  or<br>
  merge scans which stem from similar precursors (for MS2 and above).

  In any case, the number of scans will be reduced.

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_SpectraMerger.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPSpectraMerger
	: public TOPPBase
{
	public:
		TOPPSpectraMerger()
			: TOPPBase("SpectraMerger","Merges spectra (each MS level separately), increasing S/N ratios.", false)
		{
			
		}
	
	protected:
		void registerOptionsAndFlags_()
		{
			registerInputFile_("in","<file>","","Input mzML file containing the spectra.");
			setValidFormats_("in", StringList::create("mzML"));
			registerOutputFile_("out","<file>","","Output mzML file.");
			setValidFormats_("in", StringList::create("mzML"));

			registerStringOption_("merging_method", "<method>", "precursor_method", "Method of merging which should be used.");
			setValidStrings_("merging_method", StringList::create("precursor_method,block_method"));

			registerSubsection_("algorithm","Algorithm section for merging spectra");
		}

	 	Param getSubsectionDefaults_(const String& /*section*/) const
  	{
    	return SpectraMerger().getParameters();
  	}	

		ExitCodes main_(int , const char**)
		{
			//-------------------------------------------------------------
			// parsing parameters
			//-------------------------------------------------------------
			String in(getStringOption_("in"));
			String out(getStringOption_("out"));
			String merging_method(getStringOption_("merging_method"));			

			//-------------------------------------------------------------
			// reading input
			//-------------------------------------------------------------

			FileHandler fh;
      FileTypes::Type in_type = fh.getType(in);

      PeakMap exp;
      fh.loadExperiment(in, exp, in_type, log_type_);
			exp.sortSpectra();

			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------					

			SpectraMerger merger;
			merger.setParameters(getParam_().copy("algorithm:",true));
			if (merging_method == "precursor_method")
			{
				merger.mergeSpectraPrecursors(exp);				
			}
			else if (merging_method == "block_method")
			{
				merger.mergeSpectraBlockWise(exp);
			}

			//-------------------------------------------------------------
      // writing output
      //-------------------------------------------------------------
		

			fh.storeExperiment(out, exp, log_type_);
	
			return EXECUTION_OK;
		}
};


int main( int argc, const char** argv )
{
	TOPPSpectraMerger tool;
	return tool.main(argc,argv);
}
  
/// @endcond





