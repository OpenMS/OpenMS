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


#include "TOPPBase.h"

#include <OpenMS/COMPARISON/CLUSTERING/ClusterFactory.h>
#include <OpenMS/FILTERING/TRANSFORMERS/NLargest.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>
#include <OpenMS/FILTERING/TRANSFORMERS/BernNorm.h>
#include <OpenMS/FILTERING/TRANSFORMERS/ParentPeakMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Scaler.h>
#include <OpenMS/FILTERING/TRANSFORMERS/SqrtMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/ThresholdMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/WindowMower.h>

#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/KERNEL/RangeUtils.h>

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
	</UL>
	
	@ingroup TOPP
*/


// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPSpectraFilter
	: public TOPPBase
{
	public:
		TOPPSpectraFilter()
			: TOPPBase("SpectraFilter")
		{
		}
	
	protected:
		void printToolUsage_()
		{
			cerr  << endl
						<< tool_name_ << " -- applies different spectrum modification filters to the data." << endl
						<< "Version: " << VersionInfo::getVersion() << endl
						<< endl
						<< "Usage:" << endl
						<< " " << tool_name_ << " [options]" << endl
						<< endl
						<< "Options are:" << endl
						<< "  -in <file>                   input mzData file name" << endl
						<< "  -out <file>                  output mzData file name" << endl
						<< "  -rt [min]:[max]              retention time range to extract" << endl
						<< "  -filters <name>[,<name>,...] filters to apply (see --help-opt for complete list)" << endl;
		}
	
		void printToolHelpOpt_()
		{
			cerr << endl
		       << tool_name_ << endl
		       << endl
		       << "INI options:" << endl
					 << "  in        input mzData file name" << endl
					 << "  out       output mzData file name" << endl
					 << "  rt        retention time range to extract" << endl
					 << "  filters   possible spectra filters are: " << endl
					 << "            - NLargest, keeps the n most intensive peaks of each spectrum" << endl
					 << "            - Normalizer, normalizes the intensity" << endl
					 << "            - BernNorm, normalizes due to method of Bern et. al" << endl
					 << "            - ParentPeakMower, reduces the intensity of the parent peak" << endl
					 << "            - Scaler, scales the intensities" << endl
					 << "            - SqrtMower, set each intensity to the square root of the original intensity" << endl
					 << "            - ThresholdMower, removes peaks lower than a threshold intensity" << endl
					 << "            - WindowMower, keeps the biggest peaks in a sliding window" << endl
					 << "            to specify options of the filters (different from the defaults) a TOPP.ini" << endl
					 << "            file should be created with a section with special options for each filter" << endl
					 << "            (see TOPP ini for an example file). For the list of options see the " << endl
					 << "            documentation of the filters." << endl
					 << endl
					 << "INI File example section:" << endl
					 << "  <ITEM name=\"in\" value=\"input.mzData\" type=\"string\"/>" << endl
					 << "  <ITEM name=\"out\" value=\"output.mzData\" type=\"string\"/>" << endl
					 << "  <ITEM name=\"rt\" value=\":100\" type=\"string\"/>" << endl;
		}
	
		void setOptionsAndFlags_()
		{
			options_["-out"] = "out";
			options_["-in"] = "in";
			options_["-rt"] = "rt";
			// filters to apply 
			options_["-filters"] = "filters";
		}
	
		ExitCodes main_(int , char**)
		{
			//-------------------------------------------------------------
			// parameter handling
			//-------------------------------------------------------------
	
			//input file names and types
			String in = getParamAsString_("in");			
			writeDebug_(String("Input file: ") + in, 1);
	
			//output file names and types
			String out = getParamAsString_("out");
			writeDebug_(String("Output file: ") + out, 1);

			//ranges
			String mz, rt, it, tmp;
			double rt_l, rt_u;
			
			//initialize ranges
			rt_l = -1 * numeric_limits<double>::max();
			rt_u = numeric_limits<double>::max();
			
			//determine rt bounds
			rt = getParamAsString_("rt",":");
			writeDebug_(String("rt bounds: ") + rt,2);	
			
			//convert bounds to numbers
			try
			{
				//rt
				tmp = rt.prefix(':');
				if (tmp!="")
				{
					rt_l = tmp.toDouble();
				}
				tmp = rt.suffix(':');
				if (tmp!="")
				{
					rt_u = tmp.toDouble();
				}
				writeDebug_("rt lower/upper bound: " + String(rt_l) + " / " + String(rt_u),1);	
			}
			catch(Exception::ConversionError& e)
			{
				writeLog_(String("Invalid boundary '") + tmp + "' given. Aborting!");
				printUsage_();
				return ILLEGAL_PARAMETERS;			
			}
		
			// get the filternames
			vector<String> filter_names;
			String filter_command = getParamAsString_("filters");
			filter_command.split(',', filter_names);
			if (filter_names.size() == 0)
			{
				filter_names.push_back(filter_command);
			}

			ClusterFactory* cluster_factory = ClusterFactory::instance();

			// get the FactoryProduct pointers from the names
			vector<FactoryProduct*> functors;
			for (vector<String>::const_iterator it = filter_names.begin(); it != filter_names.end(); ++it)
			{
				try 
				{
					functors.push_back(cluster_factory->create(*it));
				}
				catch (Exception::Base& e)
				{
					writeLog_("Unkown filter: " + *it);
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
			}

      //-------------------------------------------------------------
      // loading input
      //-------------------------------------------------------------

      MSExperiment<> exp;
      MzDataFile f;
      f.load(in,exp);

      //-------------------------------------------------------------
      // calculations
      //-------------------------------------------------------------

			RTRange<MSExperiment< >::SpectrumType> rt_predicate(rt_l, rt_u, false);

			// for every filter
			for (vector<FactoryProduct*>::const_iterator it = functors.begin(); it != functors.end(); ++it)
			{
				String ini_location = String(tool_name_) + ":" + String(instance_number_) + ":filters:";
				Param filter_param = getParamCopy_(ini_location + (*it)->getName() + ":", true);
				writeDebug_("Used filter parameters", filter_param, 3);
				
				String filter_name = (*it)->getName();
				if (filter_name == "NLargest")
				{
					NLargest filter;
					filter.getParam().insert("", filter_param);

					// apply to every spectrum in the mzdata file
					for (MSExperiment< >::iterator sit = exp.begin(); sit != exp.end(); ++sit)
					{
						if (rt_predicate(*sit))
						{
							filter.apply(*sit);
						}
					}
					continue;
				}

				if (filter_name == "Normalizer")
				{
					Normalizer filter;
					filter.getParam().insert("", filter_param);

          // apply to every spectrum in the mzdata file
          for (MSExperiment< >::iterator sit = exp.begin(); sit != exp.end(); ++sit)
          {
            if (rt_predicate(*sit))
            {
              filter.apply(*sit);
            }
          }
					continue;
				}

				if (filter_name == "BernNorm")
				{
					BernNorm filter;
					filter.getParam().insert("", filter_param);

          // apply to every spectrum in the mzdata file
          for (MSExperiment< >::iterator sit = exp.begin(); sit != exp.end(); ++sit)
          {
            if (rt_predicate(*sit))
            {
              filter.apply(*sit);
            }
          }
					continue;
				}

				if (filter_name == "ParentPeakMower")
				{
					ParentPeakMower filter;
					filter.getParam().insert("", filter_param);

          // apply to every spectrum in the mzdata file
          for (MSExperiment< >::iterator sit = exp.begin(); sit != exp.end(); ++sit)
          {
            if (rt_predicate(*sit))
            {
              filter.apply(*sit);
            }
          }
					continue;
				}

				if (filter_name == "Scaler")
				{
					Scaler filter;
					filter.getParam().insert("", filter_param);

          // apply to every spectrum in the mzdata file
          for (MSExperiment< >::iterator sit = exp.begin(); sit != exp.end(); ++sit)
          {
            if (rt_predicate(*sit))
            {
              filter.apply(*sit);
            }
          }
					continue;
				}
				
				if (filter_name == "SqrtMower")
				{
					SqrtMower filter;
					filter.getParam().insert("", filter_param);

          // apply to every spectrum in the mzdata file
          for (MSExperiment< >::iterator sit = exp.begin(); sit != exp.end(); ++sit)
          {
            if (rt_predicate(*sit))
            {
              filter.apply(*sit);
            }
          }
					continue;
				}

				if (filter_name == "ThresholdMower")
				{
					ThresholdMower filter;
					filter.getParam().insert("", filter_param);

          // apply to every spectrum in the mzdata file
          for (MSExperiment< >::iterator sit = exp.begin(); sit != exp.end(); ++sit)
          {
            if (rt_predicate(*sit))
            {
              filter.apply(*sit);
            }
          }
					continue;
				}

				if (filter_name == "WindowMower")
				{
					WindowMower filter;
					filter.getParam().insert("", filter_param);

          // apply to every spectrum in the mzdata file
          for (MSExperiment< >::iterator sit = exp.begin(); sit != exp.end(); ++sit)
          {
            if (rt_predicate(*sit))
            {
              filter.apply(*sit);
            }
          }
					continue;
				}

			}
		
			//-------------------------------------------------------------
			// writing output
			//-------------------------------------------------------------
			
			f.store(out, exp);
			
			return OK;
		}
};

/// @endcond


int main( int argc, char ** argv )
{
	TOPPSpectraFilter tool;
	return tool.main(argc,argv);
}

