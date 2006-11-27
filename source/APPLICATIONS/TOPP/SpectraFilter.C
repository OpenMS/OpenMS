// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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


#include <OpenMS/APPLICATIONS/TOPPBase2.h>

#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/FILTERING/TRANSFORMERS/PreprocessingFunctor.h>
#include <OpenMS/FILTERING/DATAREDUCTION/DataReducer.h>

#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/KERNEL/RangeUtils.h>

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
		<LI> ThresholdMower -- removes peaks lower than a threshold intensity
		<LI> WindowMower -- keeps the biggest peaks in a sliding window
		<LI> Normalizer -- Normalizes the peaks in the spectrum with different modes
		<LI> Scaler -- Scales the peaks according to their rank
		<LI> BernNorm -- Does the Bern et al. normalization
	</UL>
	
	@ingroup TOPP
*/


// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPSpectraFilter
	: public TOPPBase2
{
	public:
		TOPPSpectraFilter()
			: TOPPBase2("SpectraFilter", "can apply several spectra filters to the spectra")
		{
		}
	
	protected:

		void registerOptionsAndFlags_()
		{
			registerStringOption_("in", "<file>", "", "input file in MzData format");
			registerStringOption_("out", "<file>", "", "output file in MzData format");
			registerStringOption_("filters", "<filter1>[,<filter2>]", "NLargest, Scaler, BernNorm, ParentPeakMower, Normalizer, SqrtMower, ThresholdMower, WindowMower", "filter to be applied");

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
				Param filter_param = getParamCopy_(getIniLocation_()+(*it)->getName()+":");
				writeDebug_("Used filter parameters", filter_param, 3);
				writeDebug_("Applying filter: " +  (*it)->getName(), 1);
				(*it)->setParam(filter_param);
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

