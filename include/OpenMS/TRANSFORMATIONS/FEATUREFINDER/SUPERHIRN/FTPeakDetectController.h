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
// $Maintainer: Florian Zeller $
// $Authors: Lukas Mueller, Markus Mueller $
// --------------------------------------------------------------------------
//
///////////////////////////////////////////////////////////////////////////
//
//  written by Lukas N Mueller, 30.3.05
//  Lukas.Mueller@imsb.biol.ethz.ch
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
//
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
//

#ifndef USE_FT_PEAK_DETECT_CONTROLLER_H
#define USE_FT_PEAK_DETECT_CONTROLLER_H

namespace OpenMS
{

class OPENMS_DLLAPI FTPeakDetectController
{

	////////////////////////////////////////////////
	// declaration of the private members:

	private:

	////////////////////////////////////////////////
	// declaration of the public members:

	// LCMS runs
	//LCMS* THIS_LCMS;
	LCMS* lcms_;
//	std::vector<SHFeature> fakeFeatureList_;
	std::vector<LCMS> lcmsRuns_;

	// paths:
	std::string targetMzXML;
	std::string SOURCE_DIR;
	std::string OUTPUT_DIR;

	public:

	typedef std::pair<double, RawData*> Map;
	typedef std::vector<Map> Vec;

//	static bool CREATE_FEATURE_ELUTION_PROFILES;
//	static bool LCelutionPeakDebugging;
//	static double LCelutionPeakMassMin;
//	static double LCelutionPeakMassMax;

//	static MS2Feature* SearchedM2Feature;

//	static bool FEATURE_FAKE_INSERTION_BASED_ON_MS2_FEATURE;

	// class destructor
	~FTPeakDetectController();

	// class constructor
	FTPeakDetectController();
	// class copy constructor
	FTPeakDetectController(const FTPeakDetectController&);

	/////////////////////////////////////////////////////////
	// function for batch processing of mzXML data
	// parses LC-MS from runs from a directory or file of raw mzXML data:
	void parseMzXMLData();

	//////////////////////////////////////////////////
	// mzXML parsing functions for a single MzXML file:
	// start the scan parsing of a mzXML file:
	void startScanParsing(Vec datavec);

	// **** for the MS1 level post processing:
	// process MS1 level data
	void process_MS1_level_data_structure( ProcessData*);
	// adds an elution peak to the LC/MS run:
	void add_raw_peak_to_LC_MS_run( LCElutionPeak* );
	// function to add the elution profile to the feature:
	void addLCelutionProfile( SHFeature* , LCElutionPeak* );

	/////////////////////////////////////////////////////////////
	// reads already paths of existing LC-MS runs in xml format into the
	// memory
	// for now, open file system for every check, but otherwise could eb done
	// in the constructor
	bool checkIfFeatureExtractionExists( std::string );

	// **** for the MS2 level post processing:
	// process MS2 level data
	void process_MS2_level_data_structure( ProcessData*);
	// processes the tracted signals on teh MS2 level
	void extract_MS2_elution_features();
	// combine the MS2 feature trace data to the MS1 features:
	void associateMS2FeatureToMS1Feature( MS2Feature*);
	// add an observed MS2 feature to the MS1 feature
	// if an observation is already there, then
	// construct a merged MS2 feature
	void addMS2FeatureToMS1Feature( MS2Feature*, SHFeature* );

	// construct here fake ms1 features based on a observed MS2 feature
	// which however could not be matched to a exiting ms1 feature
	void constructMS1FeatureFromMS2Feature( MS2Feature* );

	/////////////////////////////////////////////////////////
	// write a parsed LC/MS into directory:
	void write_out_parsed_LC_MS(LCMS*);
	// add fake MS/MS information for the MS1 feature:
	void addFakeMSMSToFeature( SHFeature* );


	//////////////////////////////////////////////////
	// overload operators:
	FTPeakDetectController& operator=(const FTPeakDetectController&);
	FTPeakDetectController& operator<=(const FTPeakDetectController&);
	FTPeakDetectController& operator>=(const FTPeakDetectController&);
	FTPeakDetectController& operator<(const FTPeakDetectController&);
	FTPeakDetectController& operator>(const FTPeakDetectController&);

	///////////////////////////////
	// start here all the get / set
	// function to access the
	// variables of the class

	// target file:
	void set_target_file(std::string IN);
	std::string get_target_file();

	// get the vector of LC/MS runs:
//	std::vector<LCMS> getParsedData();
//	bool getParsedDataEmpty();
//	std::vector<LCMS>::iterator get_parsed_DATA_START();
//	std::vector<LCMS>::iterator get_parsed_DATA_END();
	LCMS* getLCMS();
};

	inline void FTPeakDetectController::set_target_file(std::string IN)
	{
		targetMzXML = IN;
	}

	inline std::string FTPeakDetectController::get_target_file()
	{
		return targetMzXML;
	}

	/*
	// get the vector of LC/MS runs:
	inline std::vector<LCMS> FTPeakDetectController::getParsedData()
	{
		return lcmsRuns_;
	}

	inline bool FTPeakDetectController::getParsedDataEmpty()
	{
		return LC_MS_RUNS.empty();
	}

	inline std::vector<LCMS>::iterator FTPeakDetectController::get_parsed_DATA_START()
	{
		return LC_MS_RUNS.begin();
	}

	inline std::vector<LCMS>::iterator FTPeakDetectController::get_parsed_DATA_END()
	{
		return LC_MS_RUNS.end();
	}
*/

	inline LCMS* FTPeakDetectController::getLCMS()
	{
		return lcms_;
	}
} // ns

#endif

