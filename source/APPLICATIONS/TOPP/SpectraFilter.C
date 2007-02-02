// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch$
// --------------------------------------------------------------------------


#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/FILTERING/TRANSFORMERS/PreprocessingFunctor.h>

#include <OpenMS/FORMAT/MzDataFile.h>

#include <typeinfo>

using namespace OpenMS;
using namespace std;

/**
	@page SpectraFilter SpectraFilter
	
	@brief Applies different spectrum modification filters to the data.
	
	Examples of filters are:
	<UL>
		<LI> NLargest -- keeps the n most intensive peaks of each spectrum
		<LI> ParentPeakMower -- reduces the intensity of the parent peak
		<LI> SqrtMower -- set each intensity to the square root of the original intensity
		<LI> WindowMower -- keeps the biggest peaks in a sliding window
		<LI> Normalizer -- normalizes the peaks in the spectrum with different modes (to_one, to_TIC)
		<LI> Scaler -- scales the peaks according to their rank
		<LI> BernNorm -- does the Bern et al. normalization
	</UL>

	Parameters of the different filters are documented at the class documentation of each filter
	respectively. The options can be set using the ini file. Each filter might have its own section
	named by the filter name with the parameters which should be used. An example section might look
	like:

  @code
  <NODE name="NLargest"> 
    <ITEM name="n" value="100" type="float"/>
  </NODE>
	@endcode

*/


// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPSpectraFilter
	: public TOPPBase
{
	public:
		TOPPSpectraFilter()
			: TOPPBase("SpectraFilter", "can apply several spectra filters to the spectra")
		{
		}
	
	protected:

		void registerOptionsAndFlags_()
		{
			registerStringOption_("in", "<file>", "", "input file in MzData format");
			registerStringOption_("out", "<file>", "", "output file in MzData format");
			registerStringOption_("filters", "<filter1>[,<filter2>]", "", "filter to be applied");
			addEmptyLine_();
			addText_("Available filters and their parameters are:\n"
							 "  - NLargest: keeps the n most intensive peaks of each spectrum\n"
							 "    - n: the numer of peaks to keep [200]\n"
							 "  - ParentPeakMower: reduces the intensity of the unfragmented precursor peak ions\n"
							 "    - window_size: the size of the m/z window where the peaks are removed, +/- window_size [2.0]\n"
							 "    - default_charge: if the precursor has no charge set, the default charge is assumed [2]\n"
							 "    - clean_all_charge_states: set to 1 if precursor ions of all possible charge states should be removed [1]\n"
							 "    - set_to_zero: reduce the intensities of the precursor and related ions to zero [1]\n"
							 "    - reduce_by_factor: reduce the intensities by a given factor (set 'set_to_zero' to 0) [0]\n"
							 "    - factor: factor which is used to reduce the intensities if \"reduce_by_factor\" is selected [1000.0]\n"
							 "    - consider_NH3_loss: whether NH3 loss peaks from the precursor should be removed [1]\n"
							 "    - consider_H2O_loss: whether H2O loss peaks from the precursor should be removed [1]\n"
							 "  - SqrtMower: set each intensity to the square root of the original intensity\n"
							 "  - WindowMower: keeps the most abundand peaks in a sliding window\n"
							 "    - windowsize: the size of the sliding window along the m/z axis [50]\n"
							 "    - peakcount: the number of peaks that should be kept [2]\n"
							 "  - Normalizer: normalizes the peaks to a maximum of '1'\n"
							 "   - method: normalize to TIC (\"to_TIC\") or normalize to max intensity of one (\"to_one\") [to_TIC]\n"
							 "  - Scaler: scales the peaks according to their rank in terms of intensity\n"
							 "  - BernNorm: does the Bern et al. normalization\n"
							 "    - C1 - C1 value of the normalization [48.0]\n"
							 "    - C2 - C2 value of the normalization [400.0]\n"							 
							 "    - threshold - threshold of the Bern et al. normalization [0.1]");
			addEmptyLine_();
			addText_("Parameters for the filter can only be fiven in the INI file.\n"
							 "Example parameters section for the 'NLargest':\n"
							 "  <NODE name=\"NLargest\">\n"
							 "    <ITEM name=\"n\" value=\"100\" type=\"float\"/>\n"
							 "  </NODE>");
			// register one section for each algorithm
			registerSubsection_("NLargest");
			registerSubsection_("ParentPeakMower");
			registerSubsection_("SqrtMower");
			registerSubsection_("WindowMower");
			registerSubsection_("Normalizer");
			registerSubsection_("Scaler");
			registerSubsection_("BernNorm");
		}
		
		ExitCodes main_(int , char**)
		{
			//-------------------------------------------------------------
			// parameter handling
			//-------------------------------------------------------------
	
			//input/output files
			String in(getStringOption_("in"));
			String out(getStringOption_("out"));
						
			// get the filternames
			vector<String> filter_names;
			String filter_command = getStringOption_("filters");
			filter_command.split(',', filter_names);
			if (filter_names.size() == 0)
			{	
				filter_names.push_back(filter_command);
			}

			Factory<PreprocessingFunctor>* factory = Factory<PreprocessingFunctor>::instance();

			// get the FilterFunctor pointers from the names
			vector<PreprocessingFunctor*> functors;
			for (vector<String>::const_iterator it = filter_names.begin(); it != filter_names.end(); ++it)
			{
				try 
				{
					writeDebug_("Trying to get filter '" + *it + "' from factory ", 3);
					functors.push_back(factory->create(*it));
				}
				catch (Exception::Base& e)
				{
					writeLog_("Unkown filter: '" + *it + "'");
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
			}

      //-------------------------------------------------------------
      // loading input
      //-------------------------------------------------------------

      MSExperiment<> exp;
      MzDataFile f;
      f.load(in, exp);

      //-------------------------------------------------------------
      // calculations
      //-------------------------------------------------------------
			
			// for every filter
			for (vector<PreprocessingFunctor*>::iterator it = functors.begin(); it != functors.end(); ++it)
			{
				Param filter_param = getParam_().copy((*it)->getName()+":", true);
				writeDebug_("Used filter parameters", filter_param, 3);
				writeDebug_("Applying filter: " +  (*it)->getName(), 1);
				(*it)->setParameters(filter_param);
				(*it)->filterPeakMap(exp);
			}

			//-------------------------------------------------------------
			// writing output
			//-------------------------------------------------------------
			
			f.store(out, exp);
			
			return EXECUTION_OK;
		}
};

/// @endcond


int main( int argc, char ** argv )
{
	TOPPSpectraFilter tool;
	return tool.main(argc,argv);
}

