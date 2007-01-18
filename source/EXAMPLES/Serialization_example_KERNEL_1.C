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
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------

#include <fstream>
#include <string>

// include headers that implement an archive in simple text format
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

// include headers that implement an archive in simple xml format
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>

// boost.ref library
#include <boost/ref.hpp>

#define OPENMS_HAVE_SERIALIZATION 1
#include <OpenMS/FORMAT/Serialization.h>

#include <OpenMS/KERNEL/StandardTypes.h>

/////////////////////////////////////////////////////////////
// namespace OpenMS
// {
// 	/**
// 		@defgroup default_types Default MS data types
		
// 		@brief OpenMS default data types
		
//     The following classes and type definitions provide some
// 		default data types to handle one- and two-dimensional
// 		raw data, peak (stick) data, and feature data (maps).
// 		Depending on your application, the use of other internal
// 		data structures might be preferable, although we recommend
// 		to use these data types if you have no urgent need to use					
// 		something else.
		
// 		@ingroup Kernel
// 	*/

// 	//@{	
// 	/**
// 		 @brief One-dimensional raw data point
// 	 */
// 	typedef DRawDataPoint<1,KernelTraits> RawDataPoint;

// 	/**
// 		 @brief Two-dimensional raw data point
// 	*/
// 	typedef DRawDataPoint<2,KernelTraits> RawDataPoint2D;

// 	/**
// 		 @brief Spectrum consisting of raw data points, with meta information.

// 		 Meta information includes retention time and MS-level.
// 	*/
// 	typedef MSSpectrum<RawDataPoint> RawSpectrum;
// 	/**
// 		 @brief Two-dimensional map of raw data points, with meta information about experimental settings.
// 	*/
// 	typedef MSExperiment<RawDataPoint> RawMap;

// 	/**
// 		 @brief One-dimensional peak
// 	*/
// 	typedef DPeak<1, KernelTraits> Peak;
// 	/**
// 		 @brief Two-dimensional peak
// 	*/
// 	typedef DPeak<2, KernelTraits> Peak2D;
// 	/**
// 		 @brief Spectrum consisting of peaks with meta information.
// 	*/
// 	typedef MSSpectrum<Peak> PeakSpectrum;
// 	/**
// 		 @brief  Two-dimensional map of peaks, with meta information about experimental settings.
// 	*/
// 	typedef MSExperiment<Peak> PeakMap;

// 	/**
// 		 @brief Two-dimensional feature.
// 	*/
// 	typedef DFeature<2, KernelTraits> Feature;
// 	/**
// 		 @brief  Two-dimensional map of features, with meta information about experimental settings.
// 	*/
// 	typedef DFeatureMap<2,KernelTraits> FeatureMap;
// 	//@}

// }
/////////////////////////////////////////////////////////////

int main() {

	using namespace OpenMS;

  // create and open a character archive for output
  std::ofstream text_ofs("archived.txt");

  // create and open a character archive for output
  std::ofstream xml_ofs("archived.xml");

  // create class instances
	RawDataPoint raw_data_point;
	const RawDataPoint2D raw_data_point_2D;
	const Peak peak;
  Peak2D peak_2D;
  const Peak2D & peak_2D_cref = peak_2D;
	const DataValue data_value_sho(short(66));
	const DataValue data_value_int(int(666));
	const DataValue data_value_lon(long(6666));
	const DataValue data_value_flo(float(666.666));
	const DataValue data_value_dou(double(66666.66666));
	const DataValue data_value_str("sechshundertsechsundsechzig");
	std::string blabla("blablablabla");
	// const std::string blabla("blablablabla");

  // save data to archive
  {
    boost::archive::text_oarchive text_oa(text_ofs);
    // write class instance to archive
    text_oa
			<< makeConstReference(raw_data_point)
			<< raw_data_point_2D
			<< peak
			<< peak_2D_cref
			;
    // archive and stream closed when destructors are called
  }
  {
    boost::archive::xml_oarchive xml_oa(xml_ofs);
    // write class instance to archive
    xml_oa
			<< boost::serialization::make_nvp("RawDataPoint",raw_data_point)
			<< boost::serialization::make_nvp("RawDataPoint2D",raw_data_point_2D)
			<< boost::serialization::make_nvp("Peak",peak)
			<< boost::serialization::make_nvp("Peak2D",peak_2D)
			<< boost::serialization::make_nvp("string",blabla )
			<< boost::serialization::make_nvp("DataValue",data_value_sho )
			<< boost::serialization::make_nvp("DataValue",data_value_int )
			<< boost::serialization::make_nvp("DataValue",data_value_lon )
			<< boost::serialization::make_nvp("DataValue",data_value_flo )
			<< boost::serialization::make_nvp("DataValue",data_value_dou )
			<< boost::serialization::make_nvp("DataValue",data_value_str )
			;
		// archive and stream closed when destructors are called
	}

  // ... some time later restore the class instance to its orginal state
	RawDataPoint restored_text_raw_data_point;
	RawDataPoint2D restored_text_raw_data_point_2D;
	Peak restored_text_peak;
  Peak2D restored_text_peak_2D;
  {
    // create and open an archive for input
    std::ifstream text_ifs("archived.txt", std::ios::binary);
    boost::archive::text_iarchive text_ia(text_ifs);
    // read class state from archive
    text_ia
			>> restored_text_raw_data_point
			>> restored_text_raw_data_point_2D
			>> restored_text_peak
			>> restored_text_peak_2D
			;
    // archive and stream closed when destructors are called was the same
  }
  
	RawDataPoint restored_xml_raw_data_point;
	RawDataPoint2D restored_xml_raw_data_point_2D;
	Peak restored_xml_peak;
  Peak2D restored_xml_peak_2D;
	std::string restored_xml_blabla("blablablabla");
	DataValue restored_xml_data_value_sho(short(66));
	DataValue restored_xml_data_value_int(int(666));
	DataValue restored_xml_data_value_lon(long(6666));
	DataValue restored_xml_data_value_flo(float(666.666));
	DataValue restored_xml_data_value_dou(double(66666.66666));
	DataValue restored_xml_data_value_str("sechshundertsechsundsechzig");
  {
    // create and open an archive for input
    std::ifstream xml_ifs("archived.xml", std::ios::binary);
    boost::archive::xml_iarchive xml_ia(xml_ifs);
    // read class state from archive
    xml_ia
			>> boost::serialization::make_nvp("RawDataPoint",restored_xml_raw_data_point)
			>> boost::serialization::make_nvp("RawDataPoint2D",restored_xml_raw_data_point_2D)
			>> boost::serialization::make_nvp("Peak",restored_xml_peak)
			>> boost::serialization::make_nvp("Peak2D",restored_xml_peak_2D)
			>> boost::serialization::make_nvp("string",restored_xml_blabla )
			>> boost::serialization::make_nvp("DataValue",restored_xml_data_value_sho )
			>> boost::serialization::make_nvp("DataValue",restored_xml_data_value_int )
			>> boost::serialization::make_nvp("DataValue",restored_xml_data_value_lon )
			>> boost::serialization::make_nvp("DataValue",restored_xml_data_value_flo )
			>> boost::serialization::make_nvp("DataValue",restored_xml_data_value_dou )
			>> boost::serialization::make_nvp("DataValue",restored_xml_data_value_str )
			;
    // archive and stream closed when destructors are called was the same
  }

	// object tracking would not allow us to serialize them if they were mutable

	RawDataPoint const restored_text_raw_data_point_copy(restored_text_raw_data_point);
	RawDataPoint2D const restored_text_raw_data_point_2D_copy(restored_text_raw_data_point_2D);
	Peak const restored_text_peak_copy(restored_text_peak);
  Peak2D const restored_text_peak_2D_copy(restored_text_peak_2D);

	// RawDataPoint const& restored_xml_raw_data_point_copy(restored_xml_raw_data_point);
	RawDataPoint2D const& restored_xml_raw_data_point_2D_copy(restored_xml_raw_data_point_2D);
	Peak const&  restored_xml_peak_copy(restored_xml_peak);
  Peak2D const& restored_xml_peak_2D_copy(restored_xml_peak_2D);
	std::string const& restored_xml_blabla_copy(restored_xml_blabla);
	DataValue const& restored_xml_data_value_sho_copy(restored_xml_data_value_sho);
	DataValue const& restored_xml_data_value_int_copy(restored_xml_data_value_int);
	DataValue const& restored_xml_data_value_lon_copy(restored_xml_data_value_lon);
	DataValue const& restored_xml_data_value_flo_copy(restored_xml_data_value_flo);
	DataValue const& restored_xml_data_value_dou_copy(restored_xml_data_value_dou);
	DataValue const& restored_xml_data_value_str_copy(restored_xml_data_value_str);

  // create and open a character archive for output
  std::ofstream text_ofs2("archived2.txt");

  // create and open a character archive for output
  std::ofstream xml_ofs2("archived2.xml");

  // save data to archive
  {
    boost::archive::text_oarchive text_oa2(text_ofs2);
    // write class instance to archive
    text_oa2
			<< restored_text_raw_data_point_copy
			<< restored_text_raw_data_point_2D_copy
			<< restored_text_peak_copy
			<< restored_text_peak_2D_copy
			;
    // archive and stream closed when destructors are called
  }
  {
    boost::archive::xml_oarchive xml_oa2(xml_ofs2);
    // write class instance to archive
    xml_oa2 
			// << boost::serialization::make_nvp("RawDataPoint",restored_xml_raw_data_point_copy)
			<< boost::serialization::make_nvp("RawDataPoint",*const_cast< RawDataPoint const * const>(&restored_xml_raw_data_point))
			<< boost::serialization::make_nvp("RawDataPoint2D",restored_xml_raw_data_point_2D_copy)
			<< boost::serialization::make_nvp("Peak",restored_xml_peak_copy)
			<< boost::serialization::make_nvp("Peak2D",restored_xml_peak_2D_copy)
			<< boost::serialization::make_nvp("string",restored_xml_blabla_copy )
			<< boost::serialization::make_nvp("DataValue",restored_xml_data_value_sho_copy )
			<< boost::serialization::make_nvp("DataValue",restored_xml_data_value_int_copy )
			<< boost::serialization::make_nvp("DataValue",restored_xml_data_value_lon_copy )
			<< boost::serialization::make_nvp("DataValue",restored_xml_data_value_flo_copy )
			<< boost::serialization::make_nvp("DataValue",restored_xml_data_value_dou_copy )
			<< boost::serialization::make_nvp("DataValue",restored_xml_data_value_str_copy )
			;	
    // archive and stream closed when destructors are called
  }


  return 0;
}

