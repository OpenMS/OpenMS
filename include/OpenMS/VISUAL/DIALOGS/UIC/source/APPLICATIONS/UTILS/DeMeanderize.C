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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page UTILS_DeMeanderize DeMeanderize
	
	@brief Repairs MALDI experiments which were spotted line by line.

	

	<B>The command line parameters of this tool are:</B>
	@verbinclude UTILS_DeMeanderize.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPDeMeanderize
	: public TOPPBase
{
	public:
		TOPPDeMeanderize()
			: TOPPBase("DeMeanderize","Orders the spectra of MALDI spotting plates correctly.", false)
		{
			
		}
	
	protected:
		void registerOptionsAndFlags_()
		{
			registerInputFile_("in","<mzML-file>","","Input experiment file, containing the wrongly sorted spectra.");
			setValidFormats_("in", StringList::create("mzML"));
			registerOutputFile_("out","<mzML-file>","","Output experiment file with correctly sorted spectra.");
			setValidFormats_("out", StringList::create("mzML"));
			registerIntOption_("num_spots_per_row", "<integer>", 48, "Number of spots in one column, until next row is spotted.", false);
			setMinInt_("num_spots_per_row", 1);
			registerDoubleOption_("RT_distance", "<integer>", 1.0, "RT distance between two spots which is used to calculated pseudo RT.", false, true);
			setMinFloat_("RT_distance", 0.0);
		}

		ExitCodes main_(int , const char**)
		{
			//-------------------------------------------------------------
			// parsing parameters
			//-------------------------------------------------------------
			String in(getStringOption_("in"));
			String out(getStringOption_("out"));
			Size num_spots_per_row(getIntOption_("num_spots_per_row"));
			DoubleReal RT_distance(getDoubleOption_("RT_distance"));
			
			//-------------------------------------------------------------
			// reading input
			//-------------------------------------------------------------

			PeakMap exp;
			MzMLFile f;
	    f.setLogType(log_type_);
			f.load(in, exp);

			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------					

			ProgressLogger pl;
			pl.setLogType(log_type_);
			pl.startProgress(0, exp.size(), "Assigning pseudo RTs.");
			Size num_ms1(0), num_ms1_base(0), row_counter(0);
			bool row_to_reverse(false);
			DoubleReal actual_RT(0);
			for (Size i = 0; i != exp.size(); ++i)
			{
				pl.setProgress(i);
				if (row_to_reverse)
				{
					actual_RT = (DoubleReal)(num_ms1_base + (num_spots_per_row - row_counter)) * RT_distance;
					writeDebug_("RT=" + String(actual_RT) + " (modified, row_counter=" + String(row_counter) + ")", 1);
				}
				else
				{
					actual_RT = (DoubleReal)num_ms1 * RT_distance;
					writeDebug_("RT=" + String(actual_RT), 1);
				}

				exp[i].setRT(actual_RT);

				if (exp[i].getMSLevel() == 1)
				{
					if (++row_counter >= num_spots_per_row)
					{
						row_counter = 0;
						if (row_to_reverse)
						{
							row_to_reverse = false;
						}
						else
						{
							row_to_reverse = true;
						}
					}
					++num_ms1;
					if (!row_to_reverse)
					{
						num_ms1_base = num_ms1;
					}
				}
			}
			pl.endProgress();

			// sort the spectra according to their new RT
			exp.sortSpectra();

			//-------------------------------------------------------------
      // writing output
      //-------------------------------------------------------------
		
			f.store(out, exp);
	
			return EXECUTION_OK;
		}
};


int main( int argc, const char** argv )
{
	TOPPDeMeanderize tool;
	return tool.main(argc,argv);
}
  
/// @endcond





