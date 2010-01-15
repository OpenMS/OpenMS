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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FILTERING/BASELINE/MorphologicalFilter.h>
#include <OpenMS/FORMAT/PeakTypeEstimator.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_BaselineFilter BaselineFilter
	
	@brief Executes the top-hat filter to remove the baseline of an MS experiment.
	
	This nonlinear filter, known as the top-hat operator in morphological
	mathematics (see Soille, ''Morphological Image Analysis''), is independent
	of the underlying baseline shape.  It is able to detect an over brightness
	even if the environment is not uniform.  The principle is based on the
	subtraction of a signal from its opening (erosion followed by a dilation).
	The size the structuring element (here a flat line) being conditioned by the
	width of the lineament (in our case the maximum width of a mass
	spectrometric peak) to be detected.
	
	Before baseline filtering the @ref TOPP_NoiseFilter is often applied.
		
	@note The length (given in Thomson) of the structuring element should be wider than the
	maximum peak width in the raw data.

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_BaselineFilter.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPBaselineFilter
	: public TOPPBase
{
 public:
	TOPPBaselineFilter()
		: TOPPBase("BaselineFilter","Removes the baseline from profile spectra using a top-hat filter.")
	{
	}

 protected:
	void registerOptionsAndFlags_()
	{
	  	registerInputFile_("in","<file>","","input raw data file ");
			setValidFormats_("in",StringList::create("mzML"));
			registerOutputFile_("out","<file>","","output raw data file ");
	  	setValidFormats_("out",StringList::create("mzML"));
      registerDoubleOption_("struc_elem_length","<size>",3,"Length of the structuring element.",false);
      registerStringOption_("struc_elem_unit","<unit>","Thomson","Unit of 'struc_elem_length' parameter.",false);
			setValidStrings_("struc_elem_unit",StringList::create("Thomson,DataPoints"));
      registerStringOption_("method","<string>","tophat","The name of the morphological filter to be applied. If you are unsure, use the default.",false);
			setValidStrings_("method",StringList::create("identity,erosion,dilation,opening,closing,gradient,tophat,bothat,erosion_simple,dilation_simple"));
      addEmptyLine_();
			addText_("Note: The top-hat filter works only on roughly uniform data (to generate equally-spaced data you can use the Resampler tool!)");
	}
  
 	ExitCodes main_(int , const char**)
	{
		//-------------------------------------------------------------
		// parameter handling
		//-------------------------------------------------------------
		String in = getStringOption_("in");
		String out = getStringOption_("out");

		//-------------------------------------------------------------
		// loading input
		//-------------------------------------------------------------

		MzMLFile mz_data_file;
		MSExperiment<Peak1D > ms_exp;
		mz_data_file.setLogType(log_type_);
		mz_data_file.load(in,ms_exp);

		// check for peak type (raw data required)
		if (PeakTypeEstimator().estimateType(ms_exp[0].begin(),ms_exp[0].end())==SpectrumSettings::PEAKS)
		{
			writeLog_("Warning: OpenMS peak type estimation indicates that this is not raw data!");
		}

		//check if spectra are sorted
		for (Size i=0; i< ms_exp.size(); ++i)
		{
			if (!ms_exp[i].isSorted())
			{
				writeLog_("Error: Not all spectra are sorted according to peak m/z positions. Use FileFilter to sort the input!");
				return INCOMPATIBLE_INPUT_DATA;
			}
		}

		//-------------------------------------------------------------
		// calculations
		//-------------------------------------------------------------
		MorphologicalFilter morph_filter;
    morph_filter.setLogType(log_type_);
    
    Param parameters;
    parameters.setValue("struc_elem_length",getDoubleOption_("struc_elem_length"));
    parameters.setValue("struc_elem_unit",getStringOption_("struc_elem_unit"));
    parameters.setValue("method",getStringOption_("method"));
    
    morph_filter.setParameters(parameters);
		morph_filter.filterExperiment( ms_exp );

		//-------------------------------------------------------------
		// writing output
		//-------------------------------------------------------------
		
		//annotate output with data processing info
		addDataProcessing_(ms_exp, getProcessingInfo_(DataProcessing::BASELINE_REDUCTION));

		mz_data_file.store(out,ms_exp);

		return EXECUTION_OK;
	}

};




int main( int argc, const char** argv )
{
    TOPPBaselineFilter tool;
    return tool.main(argc,argv);
}

/// @endcond
