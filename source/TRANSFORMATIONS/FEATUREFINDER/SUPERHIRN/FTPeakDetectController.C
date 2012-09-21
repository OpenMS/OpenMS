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

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <stdio.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/RawData.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MSPeak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/CentroidPeak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/CentroidData.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/SuperHirnParameters.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/LCElutionPeak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/BackgroundIntensityBin.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/BackgroundControl.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/LCMSCData.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/ProcessData.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2Info.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2Feature.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/FeatureLCProfile.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/SHFeature.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/LCMS.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS1FeatureMerger.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/FTPeakDetectController.h>

namespace OpenMS
{

	using namespace std;

//	bool FTPeakDetectController::CREATE_FEATURE_ELUTION_PROFILES = false;
//	bool FTPeakDetectController::LCelutionPeakDebugging = false;
//	double FTPeakDetectController::LCelutionPeakMassMin = -1;
//	double FTPeakDetectController::LCelutionPeakMassMax = -2;

//	MS2Feature* FTPeakDetectController::SearchedM2Feature;

// if this option is on, then construct fake features for available MS2 features
//	bool FTPeakDetectController::FEATURE_FAKE_INSERTION_BASED_ON_MS2_FEATURE = true;


////////////////////////////////////////////////
// constructor for the object FTPeakDetectController:
	FTPeakDetectController::FTPeakDetectController()
	{
		lcms_ = NULL;
	}

//////////////////////////////////////////////////
// class desctructor of FTPeakDetectController
	FTPeakDetectController::~FTPeakDetectController()
	{

		lcmsRuns_.clear();
		if (lcms_ != NULL)
		{
			delete lcms_;
			lcms_ = NULL;
		}
	}

//////////////////////////////////////////////////
// start the scan parsing of a mzXML file
	void FTPeakDetectController::startScanParsing(Vec datavec)
	{

		// set the titel of the current LC_MS run:
		string name = "tmplcms";

		// create a new LC/MS:
		lcms_ = new LCMS(name);
		lcms_->set_spectrum_ID((int) this->lcmsRuns_.size());

		ProcessData *dataProcessor = new ProcessData();
		unsigned int i;

		for (i = 0; i < datavec.size(); i++)
		{
			Map it = datavec.at(i);

			dataProcessor->setMaxScanDistance(0);
			if ((it.first >= SuperHirnParameters::instance()->getMinTR()) &&
					(it.first <= SuperHirnParameters::instance()->getMaxTR()))
			{
				SuperHirnParameters::instance()->getScanTRIndex()->insert(std::pair<int, float>(i, (float) it.first));

				// centroid it:
				CentroidData cd(SuperHirnParameters::instance()->getCentroidWindowWidth(), *it.second, it.first,
						SuperHirnParameters::instance()->centroidDataModus());

				//  store it:
				dataProcessor->add_scan_raw_data(i, it.first, &cd);

			}
		}

		//////////////////////////////////////
		// post processing of mzXML data of a file:
		// !!! MS1 LEVEL !!!
		process_MS1_level_data_structure(dataProcessor);
		// !!! MS2 LEVEL !!!
		//process_MS2_level_data_structure( reader );

		lcms_->order_by_mass();

		if (SuperHirnParameters::instance()->ms1FeatureClustering())
		{
			MS1FeatureMerger* merg = new MS1FeatureMerger(lcms_);
			merg->startFeatureMerging();
			delete merg;
		}

		lcms_->show_info();

		/* Debug-output, commented out for brutus */
		/*
		 string SEP = "";
		 FILE *file;
		 file = fopen("ffsh-features.txt","w+");
		 fprintf(file,"%s", "Features\n");

		 vector<feature>::iterator p = lcms_->get_feature_list_begin();
		 while(p != lcms_->get_feature_list_end()){
		 fprintf(file, "MS1 Feature#:%d,%s", (*p).get_feature_ID(),SEP.c_str());
		 fprintf(file, "m/z:%0.5f%s",(*p).get_MZ(),SEP.c_str());
		 fprintf(file, "[+%d],%s",(*p).get_charge_state(),SEP.c_str());
		 fprintf(file, "Area:%0.2f%s",(*p).get_peak_area(),SEP.c_str());
		 fprintf(file, ",apex:%0.2f[%0.2f:%0.2f][%d:%d:%d],s/n:%0.2f,%0.2f%s",(*p).get_retention_time(),(*p).get_retention_time_START(),(*p).get_retention_time_END(),(*p).get_scan_start(),(*p).get_scan_number(),(*p).get_scan_end(), (*p).getSignalToNoise(), (*p).get_peak_score(),SEP.c_str());
		 fprintf(file, ",matches:%d%s",(*p).get_replicate_match_nb(),SEP.c_str());
		 fprintf(file, ",LCMS-ID: %d",(*p).get_spectrum_ID());
		 fprintf(file,"%s", "\n");
		 p++;
		 }
		 fclose(file);
		 */

		lcmsRuns_.push_back(*lcms_);
		delete dataProcessor;

	}

	typedef multimap<double, MZ_series> MAIN_DATA_STRUCTURE;
	typedef MAIN_DATA_STRUCTURE::iterator MAIN_ITERATOR;

////////////////////////////////////////////////////
// process MS1 level data
	void FTPeakDetectController::process_MS1_level_data_structure(ProcessData* rawData)
	{

		// FLOFLO
		//int mzListSize = rawData->getNbMSTraces();
		//std::cout << "mzListSize: " << mzListSize << "\n";

		// extract LC elution features on the MS1 level:
		rawData->extract_elution_peaks();

		// get the new structure with the LC features:
		LCMSCData* currentData = rawData->getProcessedData();

		// iterator over the extracted features, convert
		vector<LCElutionPeak*> PEAKS = currentData->get_ALL_peak();
		// show program status:
		printf("\t* Processing of %d MS1 level features...\n", (int) PEAKS.size());

		vector<LCElutionPeak*>::iterator P = PEAKS.begin();
		while (P != PEAKS.end())
		{
			// add the LC peak for conversion to a feature structure
			// and to insert into the current LC-MS run
			add_raw_peak_to_LC_MS_run((*P));
			P++;
		}

		// order the run again:
		lcms_->order_by_mass();

		// clear:
		rawData = NULL;
		currentData = NULL;
	}

////////////////////////////////////////////////////
// add an observed MS2 feature to the MS1 feature
// if an observation is already there, then
// construct a merged MS2 feature
	void FTPeakDetectController::addMS2FeatureToMS1Feature(MS2Feature* ms2, SHFeature* ms1)
	{

		if (ms1->getMS2Feature() == NULL)
		{
			ms1->addMS2Feature(ms2);
		}
		else
		{

			MS2Feature* previousMS2 = ms1->getMS2Feature();
			previousMS2->addMS2ConsensusSpectrum(ms2);

			// in the case where a MS1 feature was generated from a MS1,
			// then adjust the retention time levels
			if (ms1->get_peak_area() == -1)
			{

				// start:
				if (ms2->getStartTR() < ms1->get_retention_time_START())
				{
					ms1->set_retention_time_START(ms2->getStartTR());
				}
				// end:
				if (ms2->getEndTR() > ms1->get_retention_time_END())
				{
					ms1->set_retention_time_END(ms2->getEndTR());
				}

			}
		}
	}

////////////////////////////////////////////////////
// construct here fake ms1 features based on a observed MS2 feature
// which however could not be matched to a exiting ms1 feature
	void FTPeakDetectController::constructMS1FeatureFromMS2Feature(MS2Feature* in)
	{

		SHFeature* fakeMS1 = new SHFeature(in);
		lcms_->add_feature(fakeMS1);
		delete fakeMS1;
		fakeMS1 = NULL;
		in = NULL;
	}

////////////////////////////////////////////////////
// adds an elution peak to the LC/MS run:
	void FTPeakDetectController::add_raw_peak_to_LC_MS_run(LCElutionPeak* PEAK)
	{

		////////////////////////////////
		// get parameters of the peak:
		// apex:
		int apex_scan = PEAK->get_scan_apex();
		double apex_MZ = PEAK->get_apex_MZ();
		double apex_TR = PEAK->get_apex_retention_time();
		float apex_INTENSITY = (float) PEAK->get_apex_intensity();

		// peak shape:
		float peak_area = (float) PEAK->get_total_peak_area();
		int charge_state = PEAK->get_charge_state();
		int peak_start = PEAK->get_start_scan();
		int peak_end = PEAK->get_end_scan();

		if ((apex_TR <= SuperHirnParameters::instance()->getMaxTR())
				&& (apex_TR >= SuperHirnParameters::instance()->getMinTR()))
		{

			////////////////////////////////////
			// construct a feature:
			SHFeature* TMP = new SHFeature(apex_MZ, apex_TR, apex_scan, peak_start, peak_end, charge_state, peak_area,
					apex_INTENSITY, 0);

			// set start / end TR:
			TMP->set_retention_time_START(PEAK->get_start_retention_time());
			TMP->set_retention_time_END(PEAK->get_end_retention_time());

			// set ids:
			TMP->set_spectrum_ID(lcms_->get_spectrum_ID());
			TMP->set_feature_ID(lcms_->get_nb_features());

			// set S/N
			TMP->setSignalToNoise(PEAK->getSignalToNoise());
			// set background noise:
			TMP->setBackgroundNoiseLevel(PEAK->getSignalToNoiseBackground());

			// set feature extract information, if available:
			if (!PEAK->getElutionPeakExtraInfo().empty())
			{
				TMP->setFeatureExtraInformation(PEAK->getElutionPeakExtraInfo());
				// add fake MS/MS information for the MS1 feature:
				addFakeMSMSToFeature(TMP);
			}

			// function to add the elution profile to the feature:
			// (if turned on)
			if (SuperHirnParameters::instance()->createFeatureElutionProfiles())
			{
				addLCelutionProfile(TMP, PEAK);
			}

			// add it to the LC_MS run:
			lcms_->add_feature(TMP);
			TMP = NULL;

		}
		PEAK = NULL;
	}

////////////////////////////////////////////////////////
// add fake MS/MS information for the MS1 feature:
	void FTPeakDetectController::addFakeMSMSToFeature(SHFeature* in)
	{

		string tmp = in->getFeatureExtraInformation();
		string tag = "INFO:";
		string sep = ";";
		// parse name etc out:
		tmp = tmp.substr(tmp.find(tag) + tag.size());

		string AC = tmp.substr(0, tmp.find(sep));
		tmp = tmp.substr(tmp.find(sep) + sep.size());

		string SQ = tmp.substr(0, tmp.find(sep));
		tmp = tmp.substr(tmp.find(sep) + sep.size());

		MS2Info* info = new MS2Info(AC, SQ, in->get_charge_state(), 1.0);
		info->set_MONO_MZ(in->get_MZ());

		info->set_SCAN_START(in->get_scan_number());
		info->set_SCAN_END(in->get_scan_number());
		info->setRetentionTime(in->get_retention_time());
		info->set_PREV_AA("R/K");

		in->add_MS2_info(info);

		delete info;
	}

////////////////////////////////////////////////////////
// function to add the elution profile to the feature:
	void FTPeakDetectController::addLCelutionProfile(SHFeature* inF, LCElutionPeak* PEAK)
	{

		////////////////////////////////
		// get parameters of the peak:
		// apex:
		int apex_scan = PEAK->get_scan_apex();
		double apex_MZ = PEAK->get_apex_MZ();
		double apex_TR = PEAK->get_apex_retention_time();
		float apex_Intensity = (float) PEAK->get_apex_intensity();

		// peak shape:
		float peak_area = (float) PEAK->get_total_peak_area();
		int charge_state = PEAK->get_charge_state();

		// create the class:
		FeatureLCProfile* myProfile = new FeatureLCProfile(apex_MZ, apex_TR, apex_Intensity, apex_scan, charge_state,
				peak_area);

		// add the raw individual mono isotopic peaks peaks:
		SIGNAL_iterator P = PEAK->get_signal_list_start();
		while (P != PEAK->get_signal_list_end())
		{

			MSPeak* MSpeak = &(*P).second;
			myProfile->addMS1elutionSignal(MSpeak->get_MZ(), MSpeak->get_intensity(), MSpeak->get_scan_number(),
					MSpeak->get_charge_state(), MSpeak->get_retention_time());
			P++;

		}

		inF->setLCelutionProfile(myProfile);

		myProfile = NULL;
	}

}
