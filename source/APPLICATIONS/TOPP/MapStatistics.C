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
// $Maintainer: Oliver Kohlbacher $
// $Authors: Oliver Kohlbacher $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>
#include <OpenMS/DATASTRUCTURES/Map.h>

#include <QtCore/QString>

#include <gsl/gsl_statistics.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_MapStatistics MapStatistics
	@brief Extract extended statistics on the features of a map for quality control.

<CENTER>
	<table>
		<tr>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
			<td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ FileInfo \f$ \longrightarrow \f$</td>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> FeatureFinder, FeatureMatcher</td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> - </td>
		</tr>
	</table>
</CENTER>

	This tool computes some basic statistics on the features of a map
	that are frequently used for quality control. In contrast to FileInfo
	
	Information displayed includes:
	- show information about the data range of a file (m/z, RT, intensity)
	- show a statistical summary for intensities, qualities, feature widths
	- break down the statistics for fractions of the map
  - total ion current included in the features as a function of RT

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_FileInfo.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

namespace OpenMS
{
	/// A little helper class to gather (and dump) some statistics from a vector<double>.  Uses statistical functions implemented in GSL.
	struct SomeStatistics
	{
		/**@brief Initialize SomeStatistics from data.

		@note: GSL statistics uses double and so we write double not DoubleReal here and where we use this.
		*/
		SomeStatistics& operator()(vector<double>& data)
		{
			// Sanity check: avoid core dump if no data points present.
      if ( !data.empty() )
			{
				sort(data.begin(), data.end());
				mean = gsl_stats_mean(&data.front(), 1, data.size());
				variance = gsl_stats_variance_m(&data.front(), 1, data.size(), mean);
				min = data.front();
				lowerq = gsl_stats_quantile_from_sorted_data (&data.front(), 1, data.size(), 0.25);
				median = gsl_stats_median_from_sorted_data(&data.front(), 1, data.size());
				upperq = gsl_stats_quantile_from_sorted_data(&data.front(), 1, data.size(), 0.75);
				max = data.back();
			}
			else
			{
				mean = variance = min = lowerq = median = upperq = max = 0.0;
			}
			return *this;
		}
		double mean, variance, min, lowerq, median, upperq, max;
	};

	/// Copy the statistics into a vector
	static vector<double>& operator << (vector<double>& result, const SomeStatistics& stats)
	{
		result.push_back(stats.mean);
		result.push_back(sqrt(stats.variance));
		result.push_back(stats.min);
		result.push_back(stats.max);
		result.push_back(stats.median);
		result.push_back(stats.lowerq);
		result.push_back(stats.upperq);
		return result;
	}

	/// Write SomeStatistics to a stream.
	static ostream& operator << (ostream& os, const SomeStatistics& rhs)
	{
		return os <<
			"  mean: " << rhs.mean << endl <<
			"  stddev: " << sqrt(rhs.variance) << endl <<
			"  median: " << rhs.median << endl <<
			"  min: " << rhs.min << endl <<
			"  max: " << rhs.max << endl;
	}
}

class TOPPMapStatistics
	: public TOPPBase
{
	public:
	TOPPMapStatistics()
		: TOPPBase("MapStatistics", "Extract extended statistics on the features of a map for quality control.")
	{
	}

	vector<double> sliceStatistics(const FeatureMap<>& map, Size begin, Size end) const
	{
		// If we are asked to produce stats for an empty set, return an empty vector.
		if (end <= begin || end > map.size()) 
		{
			return vector<double>(43);
		}
		
		Size size = end - begin;
		vector<double> intensities(size);
		vector<double> peak_widths(size);
		vector<double> mz(size);
		vector<double> overall_qualities(size);
		vector<double> mz_qualities(size);
		vector<double> rt_qualities(size);
		DoubleReal tic = 0.0;

		for (Size i = begin; i < end; ++i)
		{
			intensities[i - begin]	= map[i].getIntensity();
			mz[i - begin]						= map[i].getMZ();
			peak_widths[i - begin]	= map[i].getWidth();
			rt_qualities[i - begin]	= map[i].getQuality(Feature::RT);
			mz_qualities[i - begin] = map[i].getQuality(Feature::MZ);
			overall_qualities[i - begin] = map[i].getOverallQuality();
			tic += map[i].getIntensity();
		}

		vector<double> results;
		SomeStatistics some_statistics;
		results.reserve(43); // 6 7-number stats + tic
		results.push_back(tic);
		results << some_statistics(intensities);
		results << some_statistics(mz);
		results << some_statistics(peak_widths);
		results << some_statistics(overall_qualities);
		results << some_statistics(rt_qualities);
		results << some_statistics(mz_qualities);

		return results;
	}
	
	protected:

	virtual void registerOptionsAndFlags_()
	{
		registerInputFile_("in","<file>","","Input file");
		setValidFormats_("in", StringList::create("featureXML,consensusXML"));
		registerStringOption_("in_type", "<type>", "", "Input file type -- default: determined from file extension or content", false);
		setValidStrings_("in_type",StringList::create("featureXML,consensusXML"));
		registerOutputFile_("out","<file>","","Optional output file. If '-' or left out, the output is written to the command line.", false);

		registerIntOption_("n", "<n>", 4, // 4 slices is the default 
											 "Report separate statistics for each of n RT slices of the map.",
											 false, false);
		setMinInt_("n", 1);
		setMaxInt_("n", 100);

		registerFlag_("m", "Show meta information about the whole experiment");
		registerFlag_("p", "Shows data processing information");
		registerFlag_("s", "Computes a summary statistics of intensities, qualities, and widths");
	}

	ExitCodes outputTo(ostream& os)
	{
		//-------------------------------------------------------------
		// Parameter handling
		//-------------------------------------------------------------

		// File names
		String in = getStringOption_("in");

		// File type
		FileHandler fh;
		FileTypes::Type in_type = fh.nameToType(getStringOption_("in_type"));

		if (in_type == FileTypes::UNKNOWN)
		{
			in_type = fh.getType(in);
			writeDebug_(String("Input file type: ") + fh.typeToName(in_type), 2);
		}

		if (in_type == FileTypes::UNKNOWN)
		{
			writeLog_("Error: Could not determine input file type!");
			return PARSE_ERROR;
		}

		MSExperiment<Peak1D> exp;
		FeatureMap<> feat;
		ConsensusMap cons;

		if (in_type == FileTypes::FEATUREXML) //features
		{
			FeatureXMLFile().load(in, feat);
			feat.updateRanges();
		}
		else if (in_type == FileTypes::CONSENSUSXML) //consensus features
		{
			ConsensusXMLFile().load(in, cons);
			cons.updateRanges();
		}

		//-------------------------------------------------------------
		// meta information
		//-------------------------------------------------------------
		if (getFlag_("m"))
		{
			os << endl
				 << "-- General information --" << endl
				 << endl
				 << "file name: " << in << endl
				 << "file type: " <<  fh.typeToName(in_type) << endl;

			//basic info
			os << endl
				 << "-- Meta information --" << endl
				 << endl;

			if (in_type == FileTypes::FEATUREXML) //features
			{
				os << "Document id       : " << feat.getIdentifier() << endl << endl;
			}
			else if (in_type == FileTypes::CONSENSUSXML) //consensus features
			{
				os << "Document id       : " << cons.getIdentifier() << endl << endl;
			}
		}

		//-------------------------------------------------------------
		// data processing
		//-------------------------------------------------------------
		if (getFlag_("p"))
		{
			//basic info
			os << endl
				 << "-- Data processing information --" << endl
				 << endl;

			//get data processing info
			vector<DataProcessing> dp;
			if (in_type == FileTypes::FEATUREXML) //features
			{
				dp = feat.getDataProcessing();
			}
			else if (in_type == FileTypes::CONSENSUSXML) //consensus features
			{
				dp = cons.getDataProcessing();
			}
		}

		//-------------------------------------------------------------
		// statistics
		//-------------------------------------------------------------
		if (getFlag_("s"))		{
			//-------------------------------------------------------------
			// Content statistics
			//-------------------------------------------------------------
			Map<String,int> meta_names;
			if (in_type == FileTypes::FEATUREXML) //features
			{
				os << "Number of features: " << feat.size() << endl
					 << endl
					 << "Ranges:" << endl
					 << "  retention time:  " << String::number(feat.getMin()[Peak2D::RT],2) << " : " << String::number(feat.getMax()[Peak2D::RT],2) << endl
					 << "  mass-to-charge:  " << String::number(feat.getMin()[Peak2D::MZ],2) << " : " << String::number(feat.getMax()[Peak2D::MZ],2) << endl
					 << "  intensity:       " << String::number(feat.getMinInt(),2) << " : " << String::number(feat.getMaxInt(),2) << endl
					 << endl;

				// Charge distribution
				Map<UInt, UInt> charges;
				for (Size i = 0; i < feat.size(); ++i)
				{
					charges[feat[i].getCharge()]++;
				}

				os << "Charge distribution" << endl;
				for (Map<UInt, UInt>::const_iterator it = charges.begin(); 
						 it != charges.end(); ++it)
				{
					os << "charge " << it->first << ": " << it->second << endl;
				}
			}
			else if (in_type == FileTypes::CONSENSUSXML) //consensus features
			{
				map<Size, UInt> num_consfeat_of_size;
				for (ConsensusMap::const_iterator cmit = cons.begin(); 
						 cmit != cons.end(); ++cmit )
				{
					++num_consfeat_of_size[cmit->size()];
				}

				os << endl << "Number of consensus features:" << endl;
				for (map<Size, UInt>::reverse_iterator i = num_consfeat_of_size.rbegin(); i != num_consfeat_of_size.rend(); ++i )
				{
					os << "  of size " << setw(2) << i->first << ": " << setw(6) << i->second << endl;
				}
				os << "  total:      " << setw(6) << cons.size() << endl << endl;

				os << "Ranges:" << endl
					 << "  retention time:  " << String::number(cons.getMin()[Peak2D::RT],2) << " : " << String::number(cons.getMax()[Peak2D::RT],2) << endl
					 << "  mass-to-charge:  " << String::number(cons.getMin()[Peak2D::MZ],2) << " : " << String::number(cons.getMax()[Peak2D::MZ],2) << endl
					 << "  intensity:       " << String::number(cons.getMinInt(),2) << " : " << String::number(cons.getMaxInt(),2) << endl;

				// file descriptions
				const ConsensusMap::FileDescriptions& descs = cons.getFileDescriptions();
        if ( !descs.empty() )
				{
					os << endl <<
						"File descriptions:" << endl;
					for (ConsensusMap::FileDescriptions::const_iterator it=descs.begin(); it != descs.end(); ++it)
					{
						os << " - " << it->second.filename << endl
								 << "   identifier: " << it->first << endl
								 << "   label     : " << it->second.label << endl
								 << "   size      : " << it->second.size << endl;
					}
				}
			}
		
			os << endl
					 << "-- Summary Statistics --" << endl
					 << endl;

		}
		OpenMS::SomeStatistics some_statistics;

		if (in_type == FileTypes::FEATUREXML) //features
		{
			feat.sortByRT();

			vector<double> slice_stats;
			Size n = getIntOption_("n");

			Size begin = 0;
			Size end = 0;
			os << "#slice\tRT_begin\tRT_end\tnumber_of_features\ttic\t" 
				 << "int_mean\tint_stddev\tint_min\tint_max\tint_median\tint_lowerq\tint_upperq\t"
				 << "mz_mean\tmz_stddev\tmz_min\tmz_max\tmz_median\tmz_lowerq\tmz_upperq\t"
				 << "width_mean\twidth_stddev\twidth_min\twidth_max\twidth_median\twidth_lowerq\twidth_upperq\t"
				 << "qual_mean\tqual_stddev\tqual_min\tqual_max\tqual_median\tqual_lowerq\tqual_upperq\t"
				 << "rt_qual_mean\trt_qual_stddev\trt_qual_min\trt_qual_max\trt_qual_median\trt_qual_lowerq\trt_qual_upperq\t"
				 << "mz_qual_mean\tmz_qual_stddev\tmz_qual_min\tmz_qual_max\tmz_qual_median\tmz_qual_lowerq\tmz_qual_upperq"
				 << endl;

			double rt_begin = 0.0;
			for (Size slice = 0; slice < n; ++slice)
			{
				// Determine slice boundaries.
				DoubleReal rt_end = feat.rbegin()->getRT() / (double)n * (slice + 1);
				for (end = begin; end < feat.size() && feat[end].getRT() < rt_end; ++end) {}

				// Compute statistics on all features in this slice.
				slice_stats = sliceStatistics(feat, begin, end);

				// Write the beginning and end of the slices to the output as well as the slice index.
				os << slice << "\t" << rt_begin << "\t" << rt_end << "\t" << end - begin << "\t";

				// Write the statistics as a line of an csv file
				copy(slice_stats.begin(), slice_stats.end(), ostream_iterator<double>(os, "\t"));
				os << endl;

				begin = end;
				rt_begin = rt_end;
			}
		}
		else if (in_type == FileTypes::CONSENSUSXML) //consensus features
		{
			Size size = cons.size();

			vector<double> intensities;
			intensities.reserve(size);
			vector<double> qualities(size);
			qualities.reserve(size);
			vector<double> widths(size);
			widths.reserve(size);

			vector<double> rt_delta_by_elems;
			vector<double> rt_aad_by_elems;
			vector<double> rt_aad_by_cfs;
			rt_aad_by_cfs.reserve(size);

			vector<double> mz_delta_by_elems;
			vector<double> mz_aad_by_elems;
			vector<double> mz_aad_by_cfs;
			mz_aad_by_cfs.reserve(size);

			vector<double> it_delta_by_elems;
			vector<double> it_aad_by_elems;
			vector<double> it_aad_by_cfs;
			it_aad_by_cfs.reserve(size);

			for (ConsensusMap::const_iterator cm_iter = cons.begin();
					 cm_iter != cons.end(); ++cm_iter)
			{
				double rt_aad = 0;
				double mz_aad = 0;
				double it_aad = 0;
				intensities.push_back(cm_iter->getIntensity());
				qualities.push_back(cm_iter->getQuality());
				widths.push_back(cm_iter->getWidth());
				for (ConsensusFeature::HandleSetType::const_iterator hs_iter = cm_iter->begin();
							hs_iter != cm_iter->end(); ++hs_iter)
				{
					double rt_diff = hs_iter->getRT() - cm_iter->getRT();
					rt_delta_by_elems.push_back(rt_diff);
					if (rt_diff < 0)
					{	
						rt_diff = -rt_diff;
					}
					rt_aad_by_elems.push_back(rt_diff);
					rt_aad += rt_diff;
					double mz_diff = hs_iter->getMZ() - cm_iter->getMZ();
					mz_delta_by_elems.push_back(mz_diff);
					if (mz_diff < 0)
					{
						mz_diff = -mz_diff;
					}
					mz_aad_by_elems.push_back(mz_diff);
					mz_aad += mz_diff;
					double it_ratio = hs_iter->getIntensity() / ( cm_iter->getIntensity() ? cm_iter->getIntensity() : 1. );
					it_delta_by_elems.push_back(it_ratio);
					if (it_ratio < 1.)
					{
						it_ratio = 1./it_ratio;
					}
					it_aad_by_elems.push_back(it_ratio);
					it_aad += it_ratio;
				}
				if (!cm_iter->empty())
				{
					rt_aad /= cm_iter->size();
					mz_aad /= cm_iter->size();
					it_aad /= cm_iter->size();
				} // otherwise rt_aad etc. are 0 anyway
				rt_aad_by_cfs.push_back(rt_aad);
				mz_aad_by_cfs.push_back(mz_aad);
				it_aad_by_cfs.push_back(it_aad);
			}

			os.precision(writtenDigits(ConsensusFeature::IntensityType()));
			os << "Intensities of consensus features:" << endl << some_statistics(intensities) << endl;

			os.precision(writtenDigits(ConsensusFeature::QualityType()));
			os << "Qualities of consensus features:" << endl << some_statistics(qualities) << endl;

			os.precision(writtenDigits(ConsensusFeature::CoordinateType()));
			os << "Retention time differences ( element-center, weight 1 per element):" << endl << some_statistics(rt_delta_by_elems) << endl;
			os << "Absolute retention time differences ( |element-center|, weight 1 per element):" << endl << some_statistics(rt_aad_by_elems) << endl;
			os << "Average absolute differences of retention time within consensus features ( |element-center|, weight 1 per consensus features):" << endl << some_statistics(rt_aad_by_cfs) << endl;

			os.precision(writtenDigits(ConsensusFeature::CoordinateType()));
			os << "Mass-to-charge differences ( element-center, weight 1 per element):" << endl << some_statistics(mz_delta_by_elems) << endl;
			os << "Absolute differences of mass-to-charge ( |element-center|, weight 1 per element):" << endl << some_statistics(mz_aad_by_elems) << endl;
			os << "Average absolute differences of mass-to-charge within consensus features ( |element-center|, weight 1 per consensus features):" << endl << some_statistics(mz_aad_by_cfs) << endl;

			os.precision(writtenDigits(ConsensusFeature::IntensityType()));
			os << "Intensity ratios ( element/center, weight 1 per element):" << endl << some_statistics(it_delta_by_elems) << endl;
			os << "Relative intensity error ( max{(element/center),(center/element)}, weight 1 per element):" << endl << some_statistics(it_aad_by_elems) << endl;
			os << "Average relative intensity error within consensus features ( max{(element/center),(center/element)}, weight 1 per consensus features):" << endl << some_statistics(it_aad_by_cfs) << endl;
		}

		return EXECUTION_OK;
	}

	ExitCodes main_(int, const char**)
	{
		String out = getStringOption_("out");

		//output to command line
		if (out == "" || out == "-")
		{
			return outputTo(cout);
		}
		//output to file
		else
		{
			ofstream os(out.c_str());
			return outputTo(os);
		}
	}
};

int main(int argc, const char** argv)
{
	TOPPMapStatistics tool;
	return tool.main(argc, argv);
}

/// @endcond
