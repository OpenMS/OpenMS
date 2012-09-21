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
//  PEAK DETECTION OF FOURIER TRANSFORME MS INSTRUMENT DATA
//
//  by Lukas Mueller, Lukas.Mueller@imsb.biol.ethz.ch
//  October 2005
//  
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
//
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
// 
//

#include <fstream>
#include <iostream>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/RawData.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MSPeak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/CentroidPeak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/CentroidData.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/SuperHirnParameters.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/LCElutionPeak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/BackgroundIntensityBin.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/PeptideIsotopeDistribution.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/ExternalIsotopicDistribution.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/BackgroundControl.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/LCMSCData.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/ProcessData.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/FTPeakDetecMzXmlReader.h>

namespace OpenMS
{

// debugging classes:
//int	FTPeakDetecMzXmlReader::sfReportMonoPeaks = 0; // 1 if info about monoisotopic peaks should be written to mono_peaks.txt
//std::string FTPeakDetecMzXmlReader::sfDebugDirectory; // Directory where peak detection debug files are written
//int	FTPeakDetecMzXmlReader::sfReportScanNumber = -1; // if sfReportMonoPeaks is set to 1, details about this spectrum will be written to debug files 

//std::vector<double> FTPeakDetecMzXmlReader::FRAGMENT_MASS_SCAN_LEVELS;
//std::vector<double> FTPeakDetecMzXmlReader::PEAK_EXTRACTION_SCAN_LEVELS;
//bool FTPeakDetecMzXmlReader::MS2_PEAK_PROCESSING = false;
//int FTPeakDetecMzXmlReader::MS1_base_inter_scan_distance;
//int FTPeakDetecMzXmlReader::MS2_base_inter_scan_distance;

////////////////////////////////////////////////
// constructor for the object ana_summarizer:
	FTPeakDetecMzXmlReader::FTPeakDetecMzXmlReader()
	{

		dataProcessor_ = NULL;
		//MS2_LC_MS_DATA_PROCESSOR = NULL;

		// initialize the variables
		//  offset_ = 0;
		//  scan_index = NULL;
		//  total_scan = 0;
		//  nbMS2Scans = 0;

		// FLOFLO
		dataProcessor_ = new ProcessData();

	}

//////////////////////////////////////////////////
// class desctructor
	FTPeakDetecMzXmlReader::~FTPeakDetecMzXmlReader()
	{
		if (dataProcessor_ != NULL)
		{
			delete dataProcessor_;
			dataProcessor_ = NULL;
		}
	}


//////////////////////////////////////////////////
// set indexes of the current mzXML file;

//////////////////////////////////////////////////
// reads the ms data from a mzXML file opened by teh handler
	void FTPeakDetecMzXmlReader::readData(Vec datavec)
	{

		unsigned int i;
		double minrt = datavec[0].begin()->first;
		double maxrt = datavec[datavec.size() - 1].begin()->first;
		ExternalIsotopicDistribution::initRetentionTimeSegments(minrt, maxrt);

		std::map<double, RawData*>::const_iterator it;
		for (i = 0; i < datavec.size(); i++)
		{
			Map::iterator it = datavec.at(i).begin();
			RawData* pRawData = it->second;
			std::vector<double> masses, intens;
			pRawData->get(masses, intens);
			getMsScan(i, it->first, it->second);
		}
	}

//////////////////////////////////////////////////
// get a MS scan at a given scan number within
// a mass range
	void FTPeakDetecMzXmlReader::getMsScan(off_t IN, double TR, RawData* data)
	{

		if ((TR >= SuperHirnParameters::instance()->getMinTR()) && (TR <= SuperHirnParameters::instance()->getMaxTR()))
		{  //  &&( scan_header.peaksCount > 1 )
			// build up an index scan vs retention time:
			SuperHirnParameters::instance()->getScanTRIndex()->insert(std::pair<int, float>(IN, (float) TR));
			//      insert_into_scan_TR_index(IN, (float)TR);

			// set the maximal inter-monoisotopic distance
			// for the same LC-elution peak
			// set the inter-monoistopic distance:
			//int max_scan = setInterMonoIsotopicLCDistance( IN , 1, FTPeakDetecMzXmlReader::MS1_base_inter_scan_distance);

			// superhirn bug: in fact, this was alway 0
			int max_scan = 0;
			dataProcessor_->setMaxScanDistance(max_scan);

			// process the data:
			processMS1InputData(IN, (float) TR, data);
		}

	}

///////////////////////////////////////////////////////////////////////////////////////
// process the MS1 level input data:
// - construct a RawData object with input peaks
// - centoid them
// - add them to Process_Data Structure:
	void FTPeakDetecMzXmlReader::processMS1InputData(int SCAN, float TR, RawData* data)
	{

		// centroid it:
		CentroidData cd(SuperHirnParameters::instance()->getCentroidWindowWidth(), *data, TR,
				SuperHirnParameters::instance()->centroidDataModus());

		//////////////////////////////////
		//  store it:
		dataProcessor_->add_scan_raw_data(SCAN, TR, &cd);

	}

}
