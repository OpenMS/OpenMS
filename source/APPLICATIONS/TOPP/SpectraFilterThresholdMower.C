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
// $Maintainer: Mathias Walzer$
// $Authors: $
// --------------------------------------------------------------------------


#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/FILTERING/TRANSFORMERS/ThresholdMower.h>

#include <OpenMS/FORMAT/MzMLFile.h>

#include <typeinfo>

using namespace OpenMS;
using namespace std;

/**
	@page TOPP_SpectraFilterThresholdMower SpectraFilterThresholdMower

	@brief Filters the top Peaks in the given spectra according to a given schema/thresholdset
	
<CENTER>
	<table>
		<tr>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
			<td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ SpectraFilter \f$ \longrightarrow \f$</td>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPicker </td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any tool operating on MS peak data @n (in mzML format)</td>
		</tr>
	</table>
</CENTER>



	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_SpectraFilterThresholdMower.cli

	For the parameters of the algorithm section see the class documentation: @n
		@ref OpenMS::ThresholdMower @n
*/


// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPSpectraFilterThresholdMower
	: public TOPPBase
{
	public:
    TOPPSpectraFilterThresholdMower()
      : TOPPBase("SpectraFilterThresholdMower", "Applies thresholdfilter to peak spectra.")
		{
		}

	protected:

		void registerOptionsAndFlags_()
		{
			registerInputFile_("in", "<file>", "", "input file ");
			setValidFormats_("in",StringList::create("mzML"));
			registerOutputFile_("out", "<file>", "", "output file ");
	  	setValidFormats_("out",StringList::create("mzML"));

			// register one section for each algorithm
			registerSubsection_("algorithm","Algorithm parameter subsection.");
			
		}
		
		Param getSubsectionDefaults_(const String& /*section*/) const
		{
			return ThresholdMower().getParameters();
		}

		ExitCodes main_(int , const char**)
		{
			//-------------------------------------------------------------
			// parameter handling
			//-------------------------------------------------------------

			//input/output files
			String in(getStringOption_("in"));
			String out(getStringOption_("out"));

      //-------------------------------------------------------------
      // loading input
      //-------------------------------------------------------------

      MSExperiment<> exp;
      MzMLFile f;
      f.setLogType(log_type_);
      f.load(in, exp);

      //-------------------------------------------------------------
      // if meta data arrays are present, remove them and warn
      //-------------------------------------------------------------
			if (exp.clearMetaDataArrays())
			{
				writeLog_("Warning: Spectrum meta data arrays cannot be sorted. They are deleted.");
			}

      //-------------------------------------------------------------
      // filter
      //-------------------------------------------------------------
			Param filter_param = getParam_().copy("algorithm:", true);
			writeDebug_("Used filter parameters", filter_param, 3);
			
			ThresholdMower filter;
			filter.setParameters(filter_param);
			filter.filterPeakMap(exp);
			
			//-------------------------------------------------------------
			// writing output
			//-------------------------------------------------------------

			//annotate output with data processing info
			addDataProcessing_(exp, getProcessingInfo_(DataProcessing::FILTERING));

			f.store(out, exp);

			return EXECUTION_OK;
		}
};

/// @endcond


int main( int argc, const char** argv )
{
  TOPPSpectraFilterThresholdMower tool;
	return tool.main(argc,argv);
}

