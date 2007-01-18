// // // -*- Mode: C++; tab-width: 2; -*-
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

#define OPENMS_HAVE_SERIALIZATION 1
#include <OpenMS/FORMAT/Serialization.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/DPickedPeak.h>

////////////////////////////////////////////////////////////

// In this serialization example we practice some harder exercises.  We
// serialize a DataValue, whose binary content can have many interpretations.
// Even worse, we serialize some features through pointers into the
// DPeakArray.  One is serialzed before, the other afterwards.  Have a look at
// archive2.xml to see how funnily shows up in xml ;-)
// 
// A few things are serialized to text form as well.  But then I just got bored
// of writing all this twice.  Writing xml is the harder case, anyway.
//
// IMPORTANT: We need to 'register' class feature.  Otherwise serialization of
// a polymorphic DPeakArray will fail unless a Feature has been serialized
// through a pointer elsewhere before (hmhm?).  See the documentation of boost
// serialization about that.  Such a registration is necessary because a base
// class cannot figure out what classes are derived from it.

BOOST_CLASS_EXPORT_GUID(OpenMS::Feature, "Feature")

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
  const int abracadabra_id =
    peak_2D.metaRegistry().registerName("abracadabra","this does the necessary magic","seckonds");
  peak_2D.setMetaValue("abracadabra",666);
  const Peak2D & peak_2D_cref = peak_2D;
	const DataValue data_value_sho(short(66));
	const DataValue data_value_int(int(666));
	const DataValue data_value_lon(long(6666));
	const DataValue data_value_flo(float(666.666));
	const DataValue data_value_dou(double(66666.66666));
	const DataValue data_value_str("sechshundertsechsundsechzig");
	std::string blabla("blablablabla");
	const RawSpectrum raw_spectrum;

	DPeakArray<2> dpeak_array;
	
	Feature feature;
 	feature.setPos(0,178);
 	feature.setPos(1,39);
	feature.setIntensity(353535);
	feature.setCharge(2);
	feature.setOverallQuality(38);
 	dpeak_array.push_back(feature);
	
	dpeak_array.push_back(Peak2D());
	dpeak_array.back().setPos(0,100);
	dpeak_array.back().setPos(1,1000);
 	
	feature.setPos(0,8);
 	feature.setPos(1,9);
	feature.setIntensity(6635);
	feature.setCharge(1);
	feature.setOverallQuality(399);
 	dpeak_array.push_back(feature);
		
  const Feature * feature0 = dynamic_cast<const Feature * const>(&dpeak_array[0]);
	const Feature * feature2 = dynamic_cast<const Feature * const>(&dpeak_array[2]);
	
	DPeakArray<2> dpeak_list(dpeak_array.begin(),dpeak_array.end());
	DPickedPeak<2> dpicked_peak;
	
  // save data to archive
	try
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
	catch ( const std::exception & e )
	{
		std::cerr << "Exception caught (line:" << __LINE__ << "): \"" << e.what() << "\"\n";
	}
  try
	{
    boost::archive::xml_oarchive xml_oa(xml_ofs);
   // xml_oa.register_type<Feature>();
		
    // write class instance to archive
    xml_oa
			<< boost::serialization::make_nvp("RawDataPoint",raw_data_point)
			<< boost::serialization::make_nvp("RawDataPoint2D",raw_data_point_2D)
			<< boost::serialization::make_nvp("Peak",peak)
			<< boost::serialization::make_nvp("Peak2D",peak_2D)
	    << boost::serialization::make_nvp("abracadabra_id",abracadabra_id)
      << boost::serialization::make_nvp("string",blabla )
			<< boost::serialization::make_nvp("DataValue",data_value_sho )
			<< boost::serialization::make_nvp("DataValue",data_value_int )
			<< boost::serialization::make_nvp("DataValue",data_value_lon )
			<< boost::serialization::make_nvp("DataValue",data_value_flo )
			<< boost::serialization::make_nvp("DataValue",data_value_dou )
			<< boost::serialization::make_nvp("DataValue",data_value_str )
			<< boost::serialization::make_nvp("RawSpectrum",raw_spectrum )
			<< boost::serialization::make_nvp("Feature",makeConstReference(feature0))
			<< boost::serialization::make_nvp("DPeakArray",makeConstReference(dpeak_array))
			<< boost::serialization::make_nvp("Feature",makeConstReference(feature2))
			<< boost::serialization::make_nvp("DPeakList",makeConstReference(dpeak_list))
			<< boost::serialization::make_nvp("DPickedPeak",makeConstReference(dpicked_peak))
			;
		// archive and stream closed when destructors are called
	}
	catch ( const std::exception & e )
	{
		std::cerr << "Exception caught (line:" << __LINE__ << "): \"" << e.what() << "\"\n";
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
  int restored_xml_abracadabra_id;
	std::string restored_xml_blabla("blablablabla");
	DataValue restored_xml_data_value_sho(short(66));
	DataValue restored_xml_data_value_int(int(666));
	DataValue restored_xml_data_value_lon(long(6666));
	DataValue restored_xml_data_value_flo(float(666.666));
	DataValue restored_xml_data_value_dou(double(66666.66666));
	DataValue restored_xml_data_value_str("sechshundertsechsundsechzig");
	RawSpectrum restored_xml_raw_spectrum;
	Feature * restored_xml_feature0;
	Feature * restored_xml_feature2;
	DPeakArray<2> restored_xml_dpeak_array;
	DPeakArray<2> restored_xml_dpeak_list;
	DPickedPeak<2> restored_xml_dpicked_peak;
	
	try
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
	    >> boost::serialization::make_nvp("abracadabra_id",restored_xml_abracadabra_id)
			>> boost::serialization::make_nvp("string",restored_xml_blabla )
			>> boost::serialization::make_nvp("DataValue",restored_xml_data_value_sho )
			>> boost::serialization::make_nvp("DataValue",restored_xml_data_value_int )
			>> boost::serialization::make_nvp("DataValue",restored_xml_data_value_lon )
			>> boost::serialization::make_nvp("DataValue",restored_xml_data_value_flo )
			>> boost::serialization::make_nvp("DataValue",restored_xml_data_value_dou )
			>> boost::serialization::make_nvp("DataValue",restored_xml_data_value_str )
			>> boost::serialization::make_nvp("RawSpectrum",restored_xml_raw_spectrum )
			>> boost::serialization::make_nvp("Feature",restored_xml_feature0 )
			>> boost::serialization::make_nvp("DPeakArray",restored_xml_dpeak_array )
			>> boost::serialization::make_nvp("Feature",restored_xml_feature2 )
			>> boost::serialization::make_nvp("DPeakList",restored_xml_dpeak_list )
			>> boost::serialization::make_nvp("DPickedPeak",restored_xml_dpicked_peak )
			;
		// archive and stream closed when destructors are called
	}
	catch ( const std::exception & e )
	{
		std::cerr << "Exception caught (line:" << __LINE__ << "): \"" << e.what() << "\"\n";
	}

	// create and open a character archive for output
  std::ofstream text_ofs2("archived2.txt");

  // create and open an xml archive for output
  std::ofstream xml_ofs2("archived2.xml");

  // save data to archive
	try
  {
    boost::archive::text_oarchive text_oa2(text_ofs2);
    // write class instance to archive
    text_oa2
			<< makeConstReference(restored_text_raw_data_point)
			<< makeConstReference(restored_text_raw_data_point_2D)
			<< makeConstReference(restored_text_peak)
			<< makeConstReference(restored_text_peak_2D)
			;
    // archive and stream closed when destructors are called
  }
	catch ( const std::exception & e )
	{
		std::cerr << "Exception caught (line:" << __LINE__ << "): \"" << e.what() << "\"\n";
	}
	try
  {
    boost::archive::xml_oarchive xml_oa2(xml_ofs2);
    // write class instance to archive
    xml_oa2 
			<< boost::serialization::make_nvp("RawDataPoint",makeConstReference(restored_xml_raw_data_point))
			<< boost::serialization::make_nvp("RawDataPoint2D",makeConstReference(restored_xml_raw_data_point_2D))
			<< boost::serialization::make_nvp("Peak",makeConstReference(restored_xml_peak))
			<< boost::serialization::make_nvp("Peak2D",makeConstReference(restored_xml_peak_2D))
	    << boost::serialization::make_nvp("abracadabra_id",makeConstReference(restored_xml_abracadabra_id))
			<< boost::serialization::make_nvp("string",makeConstReference(restored_xml_blabla))
			<< boost::serialization::make_nvp("DataValue",makeConstReference(restored_xml_data_value_sho))
			<< boost::serialization::make_nvp("DataValue",makeConstReference(restored_xml_data_value_int))
			<< boost::serialization::make_nvp("DataValue",makeConstReference(restored_xml_data_value_lon))
			<< boost::serialization::make_nvp("DataValue",makeConstReference(restored_xml_data_value_flo))
			<< boost::serialization::make_nvp("DataValue",makeConstReference(restored_xml_data_value_dou))
			<< boost::serialization::make_nvp("DataValue",makeConstReference(restored_xml_data_value_str))
			<< boost::serialization::make_nvp("RawSpectrum",makeConstReference(restored_xml_raw_spectrum))
			<< boost::serialization::make_nvp("Feature",makeConstReference(restored_xml_feature0))
			<< boost::serialization::make_nvp("DPeakArray",makeConstReference(restored_xml_dpeak_array))
			<< boost::serialization::make_nvp("Feature",makeConstReference(restored_xml_feature2))
			<< boost::serialization::make_nvp("DPeakList",makeConstReference(restored_xml_dpeak_list))
			<< boost::serialization::make_nvp("DPickedPeak",makeConstReference(restored_xml_dpicked_peak))
				;	
    // archive and stream closed when destructors are called
	}
	catch ( const std::exception & e )
	{
		std::cerr << "Exception caught (line:" << __LINE__ << "): \"" << e.what() << "\"\n";
	}

  return 0;
}

