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
// $Authors: Florian Zeller $
// --------------------------------------------------------------------------

#include <list>
#include <vector>
#include <map>


#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/SuperHirnParameters.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2Feature.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MSPeak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/CentroidData.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/LCElutionPeak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/BackgroundIntensityBin.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/BackgroundControl.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/LCMSCData.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/ProcessData.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2Info.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/FeatureLCProfile.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/SHFeature.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/LCMS.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/FTPeakDetectController.h>

#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmSHCtrl.h>


namespace OpenMS
{
	typedef std::pair<double, RawData*> MyMap;
	typedef std::vector<MyMap> Vec;

	std::vector<Feature> FeatureFinderAlgorithmSHCtrl::extractPeaks(Vec datavec)
	{

		SuperHirnParameters::instance()->initIsotopeDist_ = false;  // reset this so that the IsotopeDist gets reinitalized

		FTPeakDetectController controller;
		controller.startScanParsing(datavec);

		std::vector<Feature> thefeatures;

		std::vector<SHFeature>::iterator p = controller.getLCMS()->get_feature_list_begin();
		while (p != controller.getLCMS()->get_feature_list_end())
		{

			Feature f;

			double mz = (*p).get_MZ();
			f.setMZ(mz);

			int charge = (*p).get_charge_state();
			f.setCharge(charge);

			double rt = (*p).get_retention_time();
			rt *= 60.0; // convert back
			f.setRT(rt);

			double darea = (*p).get_peak_area();
			float area = (float) darea;
			f.setIntensity(area);

// ------------------------------------------------------------------------------
// Convex hull -- needs to be calculated differently according to Markus Mueller
// ------------------------------------------------------------------------------
//      FeatureLCProfile* profile = (*p).getLCelutionProfile();
//      ConvexHull2D::PointArrayType hull_points(profile->getNbLCelutionSignals());
//      
//      // the key is SCAN
//      unsigned int j = 0;
//      std::map<int, MS1Signal>::iterator lcit;
//      for (lcit = profile->getLCelutionSignalsStart(); lcit != profile->getLCelutionSignalsEnd(); lcit++) {
//        //int scan = lcit->first;
//        MS1Signal signal = lcit->second;
//        
//        hull_points[j][0] = signal.TR * 60.0; // convert back
//        hull_points[j][1] = signal.mass;
//        j++;
//      }
//      
//      ConvexHull2D hull;
//      hull.addPoints(hull_points);
//      f.getConvexHulls().push_back(hull);
// ------------------------------------------------------------------------------

			thefeatures.push_back(f);
			p++;
		}

		return thefeatures;
	}

	void FeatureFinderAlgorithmSHCtrl::initParams(Param param)
	{

		// MS1 data centroid data
		// Key: ms1:data_is_centroided_already
		// 1 == data is centroided already
		SuperHirnParameters::instance()->centroidDataModus_ = !(param.getValue("centroiding:active").toBool());

		/*    
		 //def->search_tag("Precursor detection scan levels", &vInt);
		 // Key: ms1:precursor_detection_scan_levels

		 IntList list = (IntList)(param.getValue("ms1:precursor_detection_scan_levels"));
		 for (unsigned int i = 0; i < list.size(); i++)
		 {
		 FTPeakDetecMzXmlReader::PEAK_EXTRACTION_SCAN_LEVELS.push_back(list.at(i));
		 }

		 // belongs to SPECIFIC MS2 PEAK DETECTION PARAMETERS:
		 //def->search_tag("Fragment Mass Scan levels", &vInt);
		 FTPeakDetecMzXmlReader::FRAGMENT_MASS_SCAN_LEVELS.push_back(2);

		 //def->search_tag("MS1 max inter scan distance", &INT);
		 // Key: ms1:max_inter_scan_distance
		 FTPeakDetecMzXmlReader::MS1_base_inter_scan_distance = param.getValue("ms1:max_inter_scan_distance");
		 */

		//def->search_tag("MS1 LC retention time resolution", &DB);
		// Key: ms1:tr_resolution
		SuperHirnParameters::instance()->ms1TRResolution_ = param.getValue("ms1:tr_resolution");  //0.01;

		// Key: ms1:intensity_threshold
		// float thresh = 1000;
		// def->search_tag("FT peak detect MS1 intensity min threshold", &DB);
		//    LCMSCData::intensity_min_threshold = param.getValue("ms1:intensity_threshold");;

		// NOTE: Was hardcoded to 1000, however:
		//    - def->search_tag("FT peak detect MS2 intensity min threshold", &INTENSITY_THRESHOLD);
		//    - but the above is probably something different
		SuperHirnParameters::instance()->intensityThreshold_ = param.getValue("ms1:intensity_threshold");

		// MS1 max inter scan distance
		// Key: ms1:max_inter_scan_rt_distance
		SuperHirnParameters::instance()->maxInterScanRetentionTimeDistance_ = param.getValue(
				"ms1:max_inter_scan_rt_distance");    // 0.1;

		// def->search_tag("FT peak detect MS1 min nb peak members", &min_nb_cluster_members);
		// Key: ms1:min_nb_cluster_members
		SuperHirnParameters::instance()->minNbClusterMembers_ = param.getValue("ms1:min_nb_cluster_members"); // 4

		//def->search_tag("Detectable isotope factor",&DB);
		// Key: ms1:detectable_isotope_factor
		SuperHirnParameters::instance()->detectableIsotopeFactor_ = param.getValue("ms1:detectable_isotope_factor"); // 0.05;

		// def->search_tag("IntensityCV",&DB);
		// Key: ms1:intensity_cv
		SuperHirnParameters::instance()->intensityCV_ = param.getValue("ms1:intensity_cv"); // 0.9;

		// ----------------------------------------------------------------------
		// from void FT_PEAK_DETECT_initializer::init_all(){
		// ----------------------------------------------------------------------

		//def->search_tag("Centroid window width",&INT);
		// Key: centroiding:window_width
		SuperHirnParameters::instance()->centroidWindowWidth_ = param.getValue("centroiding:window_width"); // 5;

		//def->search_tag("Absolute isotope mass precision",&DB);
		// Key: centroiding:absolute_isotope_mass_precision
		SuperHirnParameters::instance()->massTolDa_ = param.getValue("centroiding:absolute_isotope_mass_precision"); // 0.01;

		//def->search_tag("Relative isotope mass precision",&DB);
		// Key: centroiding:relative_isotope_mass_precision
		SuperHirnParameters::instance()->massTolPpm_ = param.getValue("centroiding:relative_isotope_mass_precision"); // 10;

		//def->search_tag("Minimal peak height",&DB);
		// Key: centroiding:minimal_peak_height
		SuperHirnParameters::instance()->minIntensity_ = param.getValue("centroiding:minimal_peak_height"); // 0.0;

		//def->search_tag("Min. Centroid MS Signal Intensity",&DB);
		// Key: centroiding:min_ms_signal_intensity
		SuperHirnParameters::instance()->intensityFloor_ = param.getValue("centroiding:min_ms_signal_intensity"); // 50; //in config its 50, but in CentroidData it's 1;

		/*    
		 //def->search_tag("Report mono peaks",&INT);
		 // NO OpenMS equivalent
		 FTPeakDetecMzXmlReader::sfReportMonoPeaks = 0;

		 //def->search_tag("Report scan number",&INT);
		 // NO OpenMS equivalent
		 FTPeakDetecMzXmlReader::sfReportScanNumber = 0;
		 */

		// ----------------------------------------------------------------------
		// aus initializer
		// ----------------------------------------------------------------------
		// feature parameters:
		//def->search_tag("MS1 retention time tolerance", &TMP);
		// Key: ms1:retention_time_tolerance
		// Unit: min
		SuperHirnParameters::instance()->trTol_ = param.getValue("ms1:retention_time_tolerance"); // 0.5;

		//  def->search_tag("MS1 m/z tolerance", &TMP);
		// Key: ms1:mz_tolerance
		// Unit: ppm
		SuperHirnParameters::instance()->mzTolPpm_ = param.getValue("ms1:mz_tolerance"); // 0.0;
		//		ConsensusIsotopePattern::FT_MZ_TOLERANCE = SHFeature::PPM_MZ_TOL;

		// MS2_M2_matcher parameters:
		//def->search_tag("MS2 mass matching modus", &TMP_B);
		// Key: general:ms2_mass_matching_modus
		// Unit: define which modus used to match ms2 assignments to ms1 peaks
		//                      - theoretical mass [1]  : use theoretical mass calculated from sequence
		//                      - MS1 precursor mass [0]: use measured ms1 mass of precursor ion
		// MS2Info::THEO_MATCH_MODUS = 1;

		// MS2 matching PPM parameters:
		//def->search_tag("MS2 PPM m/z tolerance", &TMP);
		// Key: general:ms2_ppm_mz_tolerance
		// Unit: ppm
		// MS2Info::MS2_MZ_PPM_TOLERANCE = 30;

		// MS2 retention time tolerance:
		//def->search_tag("MS2 retention time tolerance", &TMP);
		// Key: general:ms2_retention_time_tolerance
		// Unit: retention time tolerance with which MS2 identifications will be associated
		//					to a defined MS1 LC elution peak [min]
		//					(if set to -1, then the MS1 retention time tolerance will be used
		//if( TMP > 0 ){
		//  MS2Info::MS2_TR_TOL = TMP;
		//}
		//else{
		//		MS2Info::MS2_TR_TOL = SHFeature::TR_TOL;
		//}

		//  def->search_tag("Peptide Prophet Threshold", &TMP);
		// Key: general:peptide_prophet_threshold
		// only MS2 related
		//double d = param.getValue("general:peptide_prophet_threshold"); // 0.9;

		//		double d = 0.9;
		// ------ peptide_DELTA_group::PEPTIDE_PROBABILITY_THRESHOLD = d;
		//		SHFeature::PEPTIDE_PROBABILITY_THRESHOLD = d;
		// PKUNSZT: commented out: the default is 0.9 already.
		// SuperHirnParameters::instance()->peptideProbabilityThreshold_ = d;
		// ------ interact_parser::PEPTIDE_PROBABILITY_THRESHOLD = d;
		//		LCMS::PEP_PROPHET_THERSHOLD = d;    // THIS IS NEVER USED

		//def->search_tag("Create monoisotopic LC profile", &TMP_B);
		// ------- LCMSDataImporter::CREATE_FEATURE_ELUTION_PROFILES = 1;
		//		FTPeakDetectController::CREATE_FEATURE_ELUTION_PROFILES = 1;
		SuperHirnParameters::instance()->createFeatureElutionProfiles_ = true;
		// ------- LC_MS_XML_writer::STORE_FEATURE_ELUTION_PROFILES = 1;

		/////////////////////////////////////////////////
		// Parameters for the peak merging:
		//def->search_tag("Activation of MS1 feature merging post processing", &TMP_B);
		// Key: ms1_feature_merger:active
//		MS1FeatureMerger::MS1_FEATURE_CLUSTERING = param.getValue("ms1_feature_merger:active").toBool(); //1;
		SuperHirnParameters::instance()->ms1FeatureClustering_ = param.getValue("ms1_feature_merger:active").toBool(); //1;

		//def->search_tag("MS1 LC retention time resolution", &TMP); // belongs to MS1 PEAK DETECTION PARAMETERS FOR THE DIFFERENT FILTER METHODS:
		// Key: ms1_feature_merger:tr_resolution
//		MS1FeatureMerger::MS1_PEAK_AREA_TR_RESOLUTION = param.getValue("ms1_feature_merger:tr_resolution"); //0.01;
		SuperHirnParameters::instance()->ms1PeakAreaTrResolution_ = param.getValue("ms1_feature_merger:tr_resolution"); //0.01;

		//def->search_tag("Initial Apex Tr tolerance", &TMP);
		// Key: ms1_feature_merger:initial_apex_tr_tolerance
//		MS1FeatureMerger::INITIAL_TR_TOLERANCE = param.getValue("ms1_feature_merger:initial_apex_tr_tolerance"); //5.0;
		SuperHirnParameters::instance()->initialTrTolerance_ = param.getValue(
				"ms1_feature_merger:initial_apex_tr_tolerance"); //5.0;

		//def->search_tag("MS1 feature Tr merging tolerance", &TMP);
		// Key: ms1_feature_merger:feature_merging_tr_tolerance
//		MS1FeatureMerger::MS1_FEATURE_MERGING_TR_TOLERANCE = param.getValue(
//				"ms1_feature_merger:feature_merging_tr_tolerance"); //1.0;
		SuperHirnParameters::instance()->ms1FeatureMergingTrTolerance_ = param.getValue(
				"ms1_feature_merger:feature_merging_tr_tolerance"); //1.0;

		//def->search_tag("Percentage of intensity variation between LC border peaks", &TMP);
		// Key: ms1_feature_merger:intensity_variation_percentage
//		MS1FeatureMerger::PERCENTAGE_INTENSITY_ELUTION_BORDER_VARIATION = param.getValue(
//				"ms1_feature_merger:intensity_variation_percentage"); //25;
		SuperHirnParameters::instance()->percentageIntensityElutionBorderVariation_ = param.getValue(
				"ms1_feature_merger:intensity_variation_percentage"); //25;

		//def->search_tag("PPM value for the m/z clustering of merging candidates", &TMP);
		// Key: ms1_feature_merger:ppm_tolerance_for_mz_clustering
		SuperHirnParameters::instance()->ppmToleranceForMZClustering_ = param.getValue(
				"ms1_feature_merger:ppm_tolerance_for_mz_clustering"); //10;

		/////////////////////////////////////////////////
		// what information is extracted from the LC/MS or mastermap:
		//  def->search_tag("start elution window", &TMP);
		// Key: ms1_feature_selection_options:start_elution_window
		// Unit: min
		SuperHirnParameters::instance()->minTR_ = param.getValue("ms1_feature_selection_options:start_elution_window"); //0;

		//  def->search_tag("end elution window", &TMP);
		// Key: ms1_feature_selection_options:end_elution_window
		// Unit: min
		SuperHirnParameters::instance()->maxTR_ = param.getValue("ms1_feature_selection_options:end_elution_window"); //180;

		//def->search_tag("MS1 feature mz range min", &TMP);
		// Key: ms1_feature_selection_options:mz_range_min
		SuperHirnParameters::instance()->minFeatureMZ_ = param.getValue("ms1_feature_selection_options:mz_range_min"); //0;

		//def->search_tag("MS1 feature mz range max", &TMP );
		// Key: ms1_feature_selection_options:mz_range_max
		SuperHirnParameters::instance()->maxFeatureMZ_ = param.getValue("ms1_feature_selection_options:mz_range_max"); //2000;

		//def->search_tag("MS1 feature CHRG range min", &TMP_I );
		// Key: ms1_feature_selection_options:chrg_range_min
		SuperHirnParameters::instance()->minFeatureChrg_ = param.getValue("ms1_feature_selection_options:chrg_range_min"); //1;
		//Deisotoper::sfMinCharge = param.getValue("ms1_feature_selection_options:chrg_range_min"); //1;

		//def->search_tag("MS1 feature CHRG range max", &TMP_I );
		// Key: ms1_feature_selection_options:chrg_range_max
		SuperHirnParameters::instance()->maxFeatureChrg_ = param.getValue("ms1_feature_selection_options:chrg_range_max"); //5;
		//Deisotoper::sfMaxCharge = param.getValue("ms1_feature_selection_options:chrg_range_max"); //5;

		/////////////////////////////////////////////////
		// what and how data is stored during superhirn processing:
		// ms2 information of a feature:
		// only the best ms2 info / feature stored:
		//def->search_tag("progress all low probability ms2 info in MS1 feature", &TMP_B);
		// SHFeature::STORE_ALL_LOW_PROBABILITY_MS2_SCANS = 0;
	}
}
