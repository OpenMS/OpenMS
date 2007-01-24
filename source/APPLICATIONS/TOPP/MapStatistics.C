// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//				   OpenMS Mass Spectrometry Framework
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
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/DFeatureMapFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/DFeatureMap.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#include <map>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page MapStatistics MapStatistics
	
	@brief Computes a five-number summary of
	intensities in raw data, picked peak or feature map.
		
	This TOPP module computes a five-number summary
	of the feature intensities and qualities in a map.
	
	The five-number summary consists of median, upper
	and lower quartile, minimum and maximum. These values
	are computed for qualities and intensities. They 
	give a measure of spread and location and are stored
	in a XML format for further processing.
	
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES 

class MapStatistics
  : public TOPPBase
{
public:
	MapStatistics()
	  : TOPPBase("MapStatistics","Computes a five-number summary of peak intensities in a LC-MS map")
	{
	}

protected:

	void registerOptionsAndFlags_()
	{
		registerStringOption_("in","<file>","","input file (feature or raw data map)");
		registerStringOption_("in_type","<type>","","input file type (default: determined from input file extension)\n"
			                                          "Valid types are: 'mzData', 'mzXML', 'DTA2D', 'ANDIMS' (cdf) , 'FeatureFile'", false);
		registerStringOption_("out","<file>","","output file in XML format");
		addEmptyLine_();
		addText_("This TOPP application can be applied to raw, picked (centroided) data and feature maps.");	
	}
	
	ExitCodes main_(int , char**)
	{

		//-------------------------------------------------------------
		// parameter handling
		//-------------------------------------------------------------
		String in = getStringOption_("in");
		String out = getStringOption_("out");
		
		FileHandler fh;
		FileHandler::Type in_type = fh.nameToType(getStringOption_("in_type"));
		
		if (in_type==FileHandler::UNKNOWN)
		{
			in_type = fh.getTypeByFileName(in);
			writeDebug_(String("Input file type (from file extention): ") + fh.typeToName(in_type), 1);
		}

		if (in_type==FileHandler::UNKNOWN)
		{
			in_type = fh.getTypeByContent(in);
			writeDebug_(String("Input file type (from content): ") + fh.typeToName(in_type), 2);
		}

		if (in_type==FileHandler::UNKNOWN)
		{
			writeLog_("Error: Could not determine input file type!");
			return PARSE_ERROR;
		}
		
		//-------------------------------------------------------------
		// calculations
		//-------------------------------------------------------------

		if (in_type == FileHandler::FEATURE)
		{
			DFeatureMap<2> map;
			DFeatureMapFile().load(in,map);

			unsigned int size = map.size();

			typedef DFeatureMap<2>::FeatureType::IntensityType IntensityType;
			typedef DFeatureMap<2>::FeatureType::QualityType QualityType;

			IntensityType * intensities = new IntensityType[ size ];
			QualityType * qualities	 = new QualityType[ size ];

			for (unsigned int i = 0; i < size; 	++i)
			{
				intensities[i] = map.at(i).getIntensity();
				qualities[i]   = map.at(i).getOverallQuality();
			}

			gsl_sort(intensities, 1, size);
			gsl_sort(qualities, 1, size);

			double mean_int, var_int, max_int, min_int;
			mean_int = gsl_stats_mean(intensities,1,size);
			var_int  = gsl_stats_variance(intensities,1,size);
			max_int  = gsl_stats_max(intensities,1,size);
			min_int  = gsl_stats_min(intensities,1,size);

			double mean_q, var_q, max_q, min_q;
			mean_q = gsl_stats_mean(qualities,1,size);
			var_q  = gsl_stats_variance(qualities,1,size);
			max_q  = gsl_stats_max(qualities,1,size);
			min_q  = gsl_stats_min(qualities,1,size);

			double median_int, upperq_int, lowerq_int;
			median_int = gsl_stats_median_from_sorted_data(intensities,1,size);
			upperq_int = gsl_stats_quantile_from_sorted_data(intensities,1,size,0.75);
			lowerq_int = gsl_stats_quantile_from_sorted_data (intensities,1,size,0.25);

			double median_q, upperq_q, lowerq_q;
			median_q = gsl_stats_median_from_sorted_data(qualities,1,size);
			upperq_q = gsl_stats_quantile_from_sorted_data(qualities,1,size,0.75);
			lowerq_q = gsl_stats_quantile_from_sorted_data (qualities,1,size,0.25);

			delete [] intensities;
			delete [] qualities;

			ofstream outstream(out.c_str());
			outstream << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>" << endl;
			outstream << "<mapstatistics>" << endl;

			outstream << "\t<intensities>" << endl;
			outstream << "\t\t<mean>" << mean_int << "</mean>" << endl;
			outstream << "\t\t<median>" << median_int << "</median>" << endl;
			outstream << "\t\t<variance>" << var_int << "</variance>" << endl;
			outstream << "\t\t<min>" << min_int << "</min>" << endl;
			outstream << "\t\t<max>" << max_int << "</max>" << endl;
			outstream << "\t\t<lower_quartile>" << lowerq_int << "</lower_quartile>" << endl;
			outstream << "\t\t<upper_quartile>" << upperq_int << "</upper_quartile>" << endl;
			outstream << "\t</intensities>" << endl;

			outstream << "\t<qualities>" << endl;
			outstream << "\t\t<mean>" << mean_q << "</mean>" << endl;
			outstream << "\t\t<median>" << median_q << "</median>" << endl;
			outstream << "\t\t<variance>" << var_q << "</variance>" << endl;
			outstream << "\t\t<min>" << min_q << "</min>" << endl;
			outstream << "\t\t<max>" << max_q << "</max>" << endl;
			outstream << "\t\t<lower_quartile>" << lowerq_q << "</lower_quartile>" << endl;
			outstream << "\t\t<upper_quartile>" << upperq_q << "</upper_quartile>" << endl;
			outstream << "\t</qualities>" << endl;

			outstream << "</mapstatistics>" << endl;

			outstream.close();
		}
		else
		{
			RawMap exp;

			if (! fh.loadExperiment(in,exp,in_type) )
			{
				writeLog_("Unsupported or corrupt input file. Aborting!");
				printUsage_();
				return ILLEGAL_PARAMETERS;			
			}
			
			//copy intensities of  MS-level 1 peaks
			exp.updateRanges(1);
			unsigned int size = exp.getSize();
			DPeak<1>::IntensityType * intensities = new  DPeak<1>::IntensityType[ size ];
			UnsignedInt i = 0;
      for (RawMap::const_iterator spec = exp.begin(); spec != exp.end(); ++spec)
      {
	      if (spec->getMSLevel()!=1)
	      {
	          continue;
	      }
	      for (RawMap::SpectrumType::const_iterator it = spec->begin(); it!=spec->end(); ++it)
	      {
					intensities[i++] = it->getIntensity();
	      }
      }

			gsl_sort(intensities, 1, size);

			double mean_int, var_int, max_int, min_int;
			mean_int = gsl_stats_mean(intensities,1,size);
			var_int  = gsl_stats_variance(intensities,1,size);
			max_int  = gsl_stats_max(intensities,1,size);
			min_int  = gsl_stats_min(intensities,1,size);

			double median_int, upperq_int, lowerq_int;
			median_int = gsl_stats_median_from_sorted_data(intensities,1,size);
			upperq_int = gsl_stats_quantile_from_sorted_data(intensities,1,size,0.75);
			lowerq_int = gsl_stats_quantile_from_sorted_data (intensities,1,size,0.25);

			delete [] intensities;

			ofstream outstream(out.c_str());
			outstream << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>" << endl;
			outstream << "<mapstatistics>" << endl;

			outstream << "\t<intensities>" << endl;
			outstream << "\t\t<mean>" << mean_int << "</mean>" << endl;
			outstream << "\t\t<median>" << median_int << "</median>" << endl;
			outstream << "\t\t<variance>" << var_int << "</variance>" << endl;
			outstream << "\t\t<min>" << min_int << "</min>" << endl;
			outstream << "\t\t<max>" << max_int << "</max>" << endl;
			outstream << "\t\t<lower_quartile>" << lowerq_int << "</lower_quartile>" << endl;
			outstream << "\t\t<upper_quartile>" << upperq_int << "</upper_quartile>" << endl;
			outstream << "\t</intensities>" << endl;

			outstream << "</mapstatistics>" << endl;

			outstream.close();

		}

		return EXECUTION_OK;

	}
};


/// @endcond


int main( int argc, char ** argv )
{
  MapStatistics tool;
  return tool.main(argc,argv);
}

