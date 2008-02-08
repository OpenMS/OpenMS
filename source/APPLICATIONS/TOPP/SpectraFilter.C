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
	respectively. The options can be set using the ini file. Each filter has its own section
	named by the filter name with the parameters which should be used.
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
			addText_("Parameters for the filter can only be fiven in the INI file.");
			
			// register one section for each algorithm
			registerSubsection_("NLargest","Keeps the n most intensive peaks of each spectrum.");
			registerSubsection_("ParentPeakMower","Reduces the intensity of the unfragmented precursor peak ions.");
			registerSubsection_("WindowMower","Keeps the most abundand peaks in a sliding window.");
			registerSubsection_("Normalizer","Normalizes the peaks to a maximum of '1'.");
			registerSubsection_("BernNorm","Does the Bern et al. normalization.");
		}
		
		Param getSubsectionDefaults_(const String& section) const
		{
			return Factory<PreprocessingFunctor>::create(section)->getDefaults();
		}
		
		ExitCodes main_(int , const char**)
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

			// get the FilterFunctor pointers from the names
			vector<PreprocessingFunctor*> functors;
			for (vector<String>::const_iterator it = filter_names.begin(); it != filter_names.end(); ++it)
			{
				try 
				{
					writeDebug_("Trying to get filter '" + *it + "' from factory ", 3);
					functors.push_back(Factory<PreprocessingFunctor>::create(*it));
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
      f.setLogType(log_type_);
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


int main( int argc, const char** argv )
{
	TOPPSpectraFilter tool;
	return tool.main(argc,argv);
}

