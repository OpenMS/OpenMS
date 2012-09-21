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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/IsotopicDist.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/RawData.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/ms_peak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/CentroidPeak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/CentroidData.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/LC_elution_peak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/BackgroundIntensityBin.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/BackgroundControl.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/Deisotoper.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/LCMSCData.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/simple_math2.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/Process_Data.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/ms2_info.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2_feature.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/featureLCprofile.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/feature.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/LC_MS.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/LC_MS_XML_reader.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS1_feature_merger.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/FT_PEAK_DETEC_mzXML_reader.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/FT_PeakDetectController.h>

#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmSHCtrl.h>

#include <vector>
#include <map>

namespace OpenMS
{ 
  typedef std::map<double, RawData*> MyMap;
  typedef std::vector<MyMap> Vec;
  
  std::vector<Feature> FeatureFinderAlgorithmSHCtrl::extractPeaks(Vec datavec) {
    std::cout << "Extracting the peaks\n";
    //initParams();
    IsotopicDist::init();
    
    FT_PeakDetectController controller;
    controller.start_scan_parsing_of_mzXML_file(datavec);
    
    std::vector<Feature> thefeatures;
    
    std::vector<feature>::iterator p = controller.THIS_LC_MS->get_feature_list_begin();
    while(p != controller.THIS_LC_MS->get_feature_list_end()){
      //fprintf(file, "MS1 Feature#:%d,%s", (*p).get_feature_ID(),SEP.c_str()); 
      //fprintf(file, "m/z:%0.5f%s",(*p).get_MZ(),SEP.c_str()); 
      //fprintf(file, "[+%d],%s",(*p).get_charge_state(),SEP.c_str()); 
      //fprintf(file, "Area:%0.2f%s",(*p).get_peak_area(),SEP.c_str()); 
      //fprintf(file, ",apex:%0.2f[%0.2f:%0.2f][%d:%d:%d],s/n:%0.2f,%0.2f%s",(*p).get_retention_time(),(*p).get_retention_time_START(),(*p).get_retention_time_END(),(*p).get_scan_start(),(*p).get_scan_number(),(*p).get_scan_end(), (*p).getSignalToNoise(), (*p).get_peak_score(),SEP.c_str()); 
      //fprintf(file, ",matches:%d%s",(*p).get_replicate_match_nb(),SEP.c_str()); 
      //fprintf(file, ",LCMS-ID: %d",(*p).get_spectrum_ID());
      //fprintf(file,"%s", "\n");
      Feature f;
      
      double mz = (*p).get_MZ();
      f.setMZ(mz);
      
      int charge = (*p).get_charge_state();
      f.setCharge(charge);
      
      double rt = (*p).get_retention_time();
      rt *= 60.0; // convert back
      f.setRT(rt);
      
      double darea = (*p).get_peak_area();
      float area = (float)darea;
      f.setIntensity(area);
      
      featureLCprofile* profile = (*p).getLCelutionProfile();
      
      ConvexHull2D::PointArrayType hull_points(profile->getNbLCelutionSignals());
      
      // das erste ist SCAN
      unsigned int j = 0;
      std::map<int, MS1Signal>::iterator lcit;
      for (lcit = profile->getLCelutionSignalsStart(); lcit != profile->getLCelutionSignalsEnd(); lcit++) {
        //int scan = lcit->first;
        MS1Signal signal = lcit->second;
        
        hull_points[j][0] = signal.TR * 60.0; // convert back
        hull_points[j][1] = signal.mass;
        j++;
      }
      
      ConvexHull2D hull;
      hull.addPoints(hull_points);
      f.getConvexHulls().push_back(hull);
      
      thefeatures.push_back(f);
      p++;
    }
    
    //Feature f;
    //set label
    //f.setMetaValue(3,plot_nr);
    //f.setCharge(c);
    //f.setOverallQuality(final_score);
    //f.setRT(fitter->getCenter());
    //f.setMZ(mono_mz);
    //f.setIntensity(0);
    //add convex hulls of mass traces
    //for (Size j=0; j<traces.size(); ++j)
    //{
    //  f.getConvexHulls().push_back(traces[j].getConvexhull());
    //}
    
    //delete controller.THIS_LC_MS;
    
    return thefeatures;
  }
  
  void FeatureFinderAlgorithmSHCtrl::initParams(Param param) {
    
    // MS1 data centroid data
    // Key: ms1:data_is_centroided_already
    // 1 == data is centroided already
    Process_Data::CENTROID_DATA_MODUS = !(param.getValue("centroiding:active").toBool());
    
    //def->search_tag("Precursor detection scan levels", &vInt);
    // Key: ms1:precursor_detection_scan_levels
    IntList list = (IntList)(param.getValue("ms1:precursor_detection_scan_levels"));
    for (unsigned int i = 0; i < list.size(); i++) 
    {
      FT_PEAK_DETEC_mzXML_reader::PEAK_EXTRACTION_SCAN_LEVELS.push_back(list.at(i));
    }
    
    // belongs to SPECIFIC MS2 PEAK DETECTION PARAMETERS:
    //def->search_tag("Fragment Mass Scan levels", &vInt);
    FT_PEAK_DETEC_mzXML_reader::FRAGMENT_MASS_SCAN_LEVELS.push_back(2);
    
    //def->search_tag("MS1 max inter scan distance", &INT);
    // Key: ms1:max_inter_scan_distance
    FT_PEAK_DETEC_mzXML_reader::MS1_base_inter_scan_distance = param.getValue("ms1:max_inter_scan_distance");
    
    //def->search_tag("MS1 LC retention time resolution", &DB);
    // Key: ms1:tr_resolution
    Process_Data::MS1_TR_RESOLUTION = param.getValue("ms1:tr_resolution");  //0.01;
    
    
    // Key: ms1:intensity_threshold
    // float thresh = 1000;
    // def->search_tag("FT peak detect MS1 intensity min threshold", &DB);
    LCMSCData::intensity_min_threshold = param.getValue("ms1:intensity_threshold");;
    
    // NOTE: Was hardcoded to 1000, however: 
    //    - def->search_tag("FT peak detect MS2 intensity min threshold", &INTENSITY_THRESHOLD);
    //    - but the above is probably something different
    Process_Data::INTENSITY_THRESHOLD = param.getValue("ms1:intensity_threshold");
    
    
    // MS1 max inter scan distance
    // Key: ms1:max_inter_scan_rt_distance
    Process_Data::max_inter_scan_retention_time_distance = param.getValue("ms1:max_inter_scan_rt_distance");    // 0.1;
    
    // def->search_tag("FT peak detect MS1 min nb peak members", &min_nb_cluster_members);
    // Key: ms1:min_nb_cluster_members
    Process_Data::min_nb_cluster_members = param.getValue("ms1:min_nb_cluster_members"); // 4
    
    //def->search_tag("Detectable isotope factor",&DB);
    // Key: ms1:detectable_isotope_factor
    IsotopicDist::sfDetectableIsoFact = param.getValue("ms1:detectable_isotope_factor"); // 0.05;
    
    // def->search_tag("IntensityCV",&DB);
    // Key: ms1:intensity_cv
    IsotopicDist::sfIntensityCV = param.getValue("ms1:intensity_cv"); // 0.9;
    
    
    // ----------------------------------------------------------------------
    // aus void FT_PEAK_DETECT_initializer::init_all(){
    // ----------------------------------------------------------------------
    
    //def->search_tag("Centroid window width",&INT);
    // Key: centroiding:window_width
    CentroidPeak::sfCentroidWindowWidth = param.getValue("centroiding:window_width"); // 5;
    
    //def->search_tag("Absolute isotope mass precision",&DB);
    // Key: centroiding:absolute_isotope_mass_precision
    CentroidData::sfMassTolDa = param.getValue("centroiding:absolute_isotope_mass_precision"); // 0.01;
    
    //def->search_tag("Relative isotope mass precision",&DB);
    // Key: centroiding:relative_isotope_mass_precision
    CentroidData::sfMassTolPpm = param.getValue("centroiding:relative_isotope_mass_precision"); // 10;
    
    //def->search_tag("Minimal peak height",&DB);
    // Key: centroiding:minimal_peak_height
    CentroidData::sfMinIntensity = param.getValue("centroiding:minimal_peak_height"); // 0.0;
    
    //def->search_tag("Min. Centroid MS Signal Intensity",&DB);
    // Key: centroiding:min_ms_signal_intensity
    CentroidData::sfIntensityFloor = param.getValue("centroiding:min_ms_signal_intensity"); // 50; //in config its 50, but in CentroidData it's 1;
    
    //def->search_tag("Report mono peaks",&INT);
    // NO OpenMS equivalent
    FT_PEAK_DETEC_mzXML_reader::sfReportMonoPeaks = 0;
    
    //def->search_tag("Report scan number",&INT);
    // NO OpenMS equivalent
    FT_PEAK_DETEC_mzXML_reader::sfReportScanNumber = 0;
    
    
    // ----------------------------------------------------------------------
    // aus initializer
    // ----------------------------------------------------------------------
    
    // feature parameters:
    //def->search_tag("MS1 retention time tolerance", &TMP);
    // Key: ms1:retention_time_tolerance
    // Unit: min
    feature::TR_TOL = param.getValue("ms1:retention_time_tolerance"); // 0.5;
    
    //  def->search_tag("MS1 m/z tolerance", &TMP);
    // Key: ms1:mz_tolerance
    // Unit: ppm
    feature::PPM_MZ_TOL = param.getValue("ms1:mz_tolerance"); // 0.0;
    consensIsotopePattern::FT_MZ_TOLERANCE = feature::PPM_MZ_TOL;
    
    // MS2_M2_matcher parameters:
    //def->search_tag("MS2 mass matching modus", &TMP_B);
    // Key: general:ms2_mass_matching_modus
    // Unit: define which modus used to match ms2 assignments to ms1 peaks
    //                      - theoretical mass [1]  : use theoretical mass calculated from sequence
    //                      - MS1 precursor mass [0]: use measured ms1 mass of precursor ion
    ms2_info::THEO_MATCH_MODUS = 1;
    
    
    /*
     
     //def->search_tag("Peptide Prophet Threshold", &TMP);
     // Key: general:peptide_prophet_threshold
     // Unit: Probability value at which peptide identification will be evaluated 
     //					(originally adopted from the Peptide Prophet, i.e. 1.0 refers to high probability)
     MS2_MS1_matcher::PEP_PROB_CUT_OFF = 0.9;
     
     //def->search_tag("MS2 SCAN tolerance", &TMP_I);
     // Key: general:ms2_scan_tolerance
     // Unit: SCAN tolerance with which MS2 identifications will be associated 
     //					to a defined MS1 LC elution peak
     MS2_MS1_matcher::SCAN_TOL = 200;
     
     //def->search_tag("INCLUSIONS LIST MS2 SCAN tolerance", &TMP_I);
     // Key: general:inclusions_list_ms2_scan_tolerance
     // Unit: SCAN tolerance with which MS2 info FROM INCLUSION LIST will be associated 
     //					to a defined MS1 LC elution peak []
     MS2_MS1_matcher::POST_SCAN_TOL = 100000;
     
     */
    
    
    // MS2 matching PPM parameters:
    //def->search_tag("MS2 PPM m/z tolerance", &TMP);
    // Key: general:ms2_ppm_mz_tolerance
    // Unit: ppm
    ms2_info::MS2_MZ_PPM_TOLERANCE = 30;
    
    // MS2 retention time tolerance:
    //def->search_tag("MS2 retention time tolerance", &TMP);
    // Key: general:ms2_retention_time_tolerance
    // Unit: retention time tolerance with which MS2 identifications will be associated 
    //					to a defined MS1 LC elution peak [min]
    //					(if set to -1, then the MS1 retention time tolerance will be used
    //if( TMP > 0 ){
    //  ms2_info::MS2_TR_TOL = TMP; 
    //}
    //else{
    ms2_info::MS2_TR_TOL = feature::TR_TOL;
    //}
    
    // peptide_DELTA_group:
    //def->search_tag("MS1 retention time tolerance", &TMP);
    // ------ peptide_DELTA_group::TR_TOL = 0.5;
    //def->search_tag("Delta pair TR tolerance", &TMP);  
    //peptide_DELTA_group::MOD_TR_TOL = TMP;
    // if smaller than -1, then use the general TR tolerance:
    //if( peptide_DELTA_group::MOD_TR_TOL < peptide_DELTA_group::TR_TOL ){
    // ------ peptide_DELTA_group::MOD_TR_TOL = peptide_DELTA_group::TR_TOL;
    //}
    //def->search_tag("Delta pair TR tolerance", &TMP);  
    
    
    // protein group:
    //def->search_tag("Peptide Proteotype Mode", &TMP_B);
    // ------ protein_group::PEPTIDE_PROTEOTYPE_MODE = 1;
    
    
    //  def->search_tag("Peptide Prophet Threshold", &TMP);
    // Key: general:peptide_prophet_threshold
    // only MS2 related
    //double d = param.getValue("general:peptide_prophet_threshold"); // 0.9;
    double d = 0.9;
    // ------ peptide_DELTA_group::PEPTIDE_PROBABILITY_THRESHOLD = d;
    feature::PEPTIDE_PROBABILITY_THRESHOLD = d;
    // ------ interact_parser::PEPTIDE_PROBABILITY_THRESHOLD = d;
    LC_MS::PEP_PROPHET_THERSHOLD = d;
    
    // FLO: TODO - what is this about?
    // slide window over TR
    //map<int, vector<double> > TMP2;
    //def->search_tag("Delta M/z list", &TMP2);
    //peptide_DELTA_group::LC_MS_Modification_Masses = TMP2; 
    //MS1_feature_ratiolizer::LC_MS_Modification_Masses = TMP2;
    //peptide_ratio_analyzer::LC_MS_Modification_Masses = TMP2;
    
    
    
    
    /* ------------------------- classes I do not have...
     
     
     // DELTA Pair MATCHING
     //def->search_tag("MS2 Delta clustering filter", &TMP_B);
     DELTA_grouper::MS_2_SELECTION = 0;
     
     
     // STATIC MODIFICATIONS:
     //def->search_tag("static glyco modifictation", &TMP_B);
     interact_parser::STATIC_GLYCOSYLATION_MOD = 0;
     //def->search_tag("static C-term. modification", &TMP);
     interact_parser::STATIC_C_TERM_MODIFICATION = 0;
     
     // FLO: TODO What kind of paramters are this...
     //map<double,double> M_TMP;
     //def->search_tag("INTERACT AA MOD transform table", &M_TMP);
     //interact_parser::MOD_MASS_TRANSFORM_TABLE = M_TMP;
     
     //def->search_tag("MS2 mass type check", &TMP_B);
     interact_parser::MS2_MASS_COMPARE = 0;
     
     //////////////////
     // log modus or not:
     bool LOG_MODE_INTENSITY = false;
     consens_profile_builder::LOG_INTENSITY_TRANSFORMATION = LOG_MODE_INTENSITY;
     profile_scorer::LOG_INTENSITY_TRANSFORMATION = LOG_MODE_INTENSITY;
     
     
     /////////////////////////
     // LC-MS correlation:
     
     //def->search_tag("pairwise correlation analysis", &TMP_B);
     LC_MS_correlation::VIEW = 0;
     
     // Tr tolerance used to compare
     // spec A to B after smoothing and
     // also used in the merging
     //def->search_tag("MS1 retention time tolerance", &TMP);
     LC_MS_correlation::tr_tol = 0.5;
     
     // M/z tolerance used to compare
     // spec A to B after smoothing and
     // also used in the merging
     //def->search_tag("MS1 m/z tolerance", &TMP);
     LC_MS_correlation::mz_tol = 10;
     
     // numbers of bins to divide the peak areas
     // used afterwards to compare peak areas of 2 peaks
     //def->search_tag("intensity bin size", &TMP_I);
     LC_MS_correlation::intensity_bin_size = 2000;
     
     // minimal alignment error = TR tolerance / 2:
     //  def->search_tag("MS1 retention time tolerance", &TMP);
     //  TMP /= 2.0;
     LC_MS_correlation::min_align_error = 0.25;
     
     // maximal alignment error:
     //def->search_tag("maximal smoothing error", &TMP);
     LC_MS_correlation::max_align_error = 3.0;
     
     // tolerance of bins, how far 2 bins can be apart and still
     // be a match
     //def->search_tag("intensity bin tolerance", &TMP);
     LC_MS_correlation::intensity_bin_tolerance = 2;
     
     // represents the worst score possible, this one will be used to 
     // normaize the observed scores between 0(bad) and 1(good) [ 0 ... 1]
     //def->search_tag("minimal LC/MS score", &TMP);
     LC_MS_correlation::min_LC_MS_score = 0.1;
     
     
     
     //////////////////////////////////////////////////////////
     // LC_MS alignmnt parameters:
     
     // retention time tolerance:
     //def->search_tag("retention time window", &TMP);
     LC_MS_aligner::Tr_window = 5.0;
     
     // mass tolerance:
     //def->search_tag("mass / charge window", &TMP);
     LC_MS_aligner::Mz_window = 20;  
     
     // percentage of outside delta error points
     //def->search_tag("perc. outside error delta points", &TMP);
     LC_MS_aligner::ERROR_DELTA_POINTS = 0.75;
     
     // Tr tolerance used to compare
     // spec A to B after smoothing and
     // also used in the merging
     //def->search_tag("MS1 retention time tolerance", &TMP);
     LC_MS_aligner::Tr_tolerance = 0.5;
     
     // M/z tolerance used to compare
     // spec A to B after smoothing and
     // also used in the merging
     //def->search_tag("MS1 m/z tolerance", &TMP);
     LC_MS_aligner::Mz_tolerance = 10;
     
     // defines if also identifications should be checked
     // in finding common lc/ms peaks between runs:
     //def->search_tag("MS2 info alignment weight", &TMP_B);
     LC_MS_aligner::ID_ALIGNMENT_WEIGHT = 1;
     
     // plotting:
     //def->search_tag("pairwise alignment plotting", &TMP_B);
     LC_MS_aligner::print_out = 0;
     
     
     /////////////////////////////////////////////////
     // these are the important parameters for the
     // for the lowess fitter druing the LC/MS alignment:
     
     // minimal alignment error = TR tolerance / 2:
     //def->search_tag("MS1 retention time tolerance", &TMP);
     regressor::min_error = 0.5;
     
     // number of boostrape cycles:
     //def->search_tag("maximal smoothing error", &TMP);
     regressor::max_error = 3.0;
     
     // max stripes for the plot smoother:
     //def->search_tag("max. nb. stripes", &TMP_I);
     plot_smoother::max_nb_stripes = 1;  
     
     // used to copmute the alignment error, use a tr window to
     // calculate the standard deviations of raw data to predicted
     // delta shift
     //def->search_tag("smoothing error TR window", &TMP);
     regressor::tr_error_smooth = 1.0;
     
     
     
     */
    
    
    /////////////////////////////////////////////////
    // parameter if the LC elution profile will be stored 
    // in the XML:
    //def->search_tag("Create monoisotopic LC profile", &TMP_B);
    // ------- LCMSDataImporter::CREATE_FEATURE_ELUTION_PROFILES = 1;
    FT_PeakDetectController::CREATE_FEATURE_ELUTION_PROFILES = 1;
    // ------- LC_MS_XML_writer::STORE_FEATURE_ELUTION_PROFILES = 1;
    
    /////////////////////////////////////////////////
    // what and how XML data is stored in the mastermap:
    // ms2 information of a feature:
    // only the best ms2 info / feature stored:
    
    /*def->search_tag("store only best MS2 per feature", &TMP_B);
     LC_MS_XML_writer::STORE_BEST_MS2_SCAN_PER_FEATURE = TMP_B;
     // only the best ms2 info / aligned feature stored:
     def->search_tag("store only best MS2 per ALIGNED feature", &TMP_B);
     LC_MS_XML_writer::STORE_BEST_MS2_SCAN_PER_FEATURE = TMP_B;
     // how many alternative protein names to store:
     def->search_tag("nb. max. alternative protein names", &TMP_I);
     LC_MS_XML_writer::MAXIMAL_NB_ALTERNATIVE_PROTEIN_NAMES = TMP_I;
     // if to store ms2 traces
     def->search_tag("MS2 fragment mass tracing", &TMP_B);
     LC_MS_XML_writer::STORE_MS2_FRAGMENT_TRACE_DATA = TMP_B;
     */
    
    
    // OK bis auf einen
    /////////////////////////////////////////////////
    // Parameters for the peak merging:
    //def->search_tag("Activation of MS1 feature merging post processing", &TMP_B);
    // Key: ms1_feature_merger:active
    MS1_feature_merger::MS1_FEATURE_CLUSTERING = param.getValue("ms1_feature_merger:active").toBool(); //1;
    
    //def->search_tag("MS1 LC retention time resolution", &TMP); // belongs to MS1 PEAK DETECTION PARAMETERS FOR THE DIFFERENT FILTER METHODS:
    // Key: ms1_feature_merger:tr_resolution
    MS1_feature_merger::MS1_PEAK_AREA_TR_RESOLUTION = param.getValue("ms1_feature_merger:tr_resolution"); //0.01;
    
    //def->search_tag("Initial Apex Tr tolerance", &TMP);
    // Key: ms1_feature_merger:initial_apex_tr_tolerance
    MS1_feature_merger::INITIAL_TR_TOLERANCE = param.getValue("ms1_feature_merger:initial_apex_tr_tolerance"); //5.0;
    
    //def->search_tag("MS1 feature Tr merging tolerance", &TMP);
    // Key: ms1_feature_merger:feature_merging_tr_tolerance
    MS1_feature_merger::MS1_FEATURE_MERGING_TR_TOLERANCE = param.getValue("ms1_feature_merger:feature_merging_tr_tolerance"); //1.0;
    
    //def->search_tag("Percentage of intensity variation between LC border peaks", &TMP);
    // Key: ms1_feature_merger:intensity_variation_percentage
    MS1_feature_merger::PERCENTAGE_INTENSITY_ELUTION_BORDER_VARIATION = param.getValue("ms1_feature_merger:intensity_variation_percentage"); //25;
    
    //def->search_tag("PPM value for the m/z clustering of merging candidates", &TMP);
    // Key: ms1_feature_merger:ppm_tolerance_for_mz_clustering
    MS1_feature_merger::PPM_TOLERANCE_FOR_MZ_CLUSTERING = param.getValue("ms1_feature_merger:ppm_tolerance_for_mz_clustering"); //10;
    
    
    
    // OK
    /////////////////////////////////////////////////
    // what information is extracted from the LC/MS or mastermap:
    // TR min:
    
    //  def->search_tag("start elution window", &TMP);
    // Key: ms1_feature_selection_options:start_elution_window
    // Unit: min
    LC_MS_XML_reader::TR_MIN = param.getValue("ms1_feature_selection_options:start_elution_window"); //0;
    FT_PEAK_DETEC_mzXML_reader::TR_MIN = param.getValue("ms1_feature_selection_options:start_elution_window"); //0; // aus initializer
    
    
    
    
    //  def->search_tag("end elution window", &TMP);
    // Key: ms1_feature_selection_options:end_elution_window
    // Unit: min
    LC_MS_XML_reader::TR_MAX = param.getValue("ms1_feature_selection_options:end_elution_window"); //180;
    FT_PEAK_DETEC_mzXML_reader::TR_MAX = param.getValue("ms1_feature_selection_options:end_elution_window"); //180; // aus initializer
    
    //def->search_tag("MS1 feature mz range min", &TMP);
    // Key: ms1_feature_selection_options:mz_range_min
    LC_MS_XML_reader::FEATURE_MZ_MIN = param.getValue("ms1_feature_selection_options:mz_range_min"); //0;
    
    //def->search_tag("MS1 feature mz range max", &TMP );
    // Key: ms1_feature_selection_options:mz_range_max
    LC_MS_XML_reader::FEATURE_MZ_MAX = param.getValue("ms1_feature_selection_options:mz_range_max"); //2000;
    
    //def->search_tag("MS1 feature signal to noise threshold", &TMP );
    // Key: ms1_feature_selection_options:signal_to_noise_threshold
    //LC_MS_XML_reader::SIGNAL_TO_NOISE_THERSHOLD = TMP;
    
    //  def->search_tag("MS1 feature intensity cutoff", &TMP );
    // Key: ms1_feature_selection_options:ms1_feature_intensity_cutoff
    //  LC_MS_XML_reader::PEAK_INTENSITY_THRESHOLD = TMP; // ich glaube 10000
    
    //def->search_tag("MS1 feature CHRG range min", &TMP_I );
    // Key: ms1_feature_selection_options:chrg_range_min
    LC_MS_XML_reader::FEATURE_CHRG_MIN = param.getValue("ms1_feature_selection_options:chrg_range_min"); //1;
    Deisotoper::sfMinCharge = param.getValue("ms1_feature_selection_options:chrg_range_min"); //1;
    
    //def->search_tag("MS1 feature CHRG range max", &TMP_I );
    // Key: ms1_feature_selection_options:chrg_range_max
    LC_MS_XML_reader::FEATURE_CHRG_MAX = param.getValue("ms1_feature_selection_options:chrg_range_max"); //5;
    Deisotoper::sfMaxCharge = param.getValue("ms1_feature_selection_options:chrg_range_max"); //5;
    
    //  Create monoisotopic LC profile:	to create and store the original profile of the detected
    //					monosiotopic pecursors in the XML (!!! increases the
    //					XML file size!!! (on[1]/off[0])
    //def->search_tag("Create monoisotopic LC profile", &TMP_B);
    //LC_MS_XML_reader::EXTRACT_MONO_ISOTOPE_PROFILE = TMP_B;
    
    
    
    /////////////////////////////////////////////////
    // what and how data is stored during superhirn processing:
    // ms2 information of a feature:
    // only the best ms2 info / feature stored:
    //def->search_tag("progress all low probability ms2 info in MS1 feature", &TMP_B);
    feature::STORE_ALL_LOW_PROBABILITY_MS2_SCANS = 0;
    
    
    // XML Data format to use during SuperHirn processing:
    //LC_MS_XML_reader::DATA_STORAGE_XML_FORMAT_TYPE = def->search_tag("SuperHirn Data Storage XML Output Format");
  }
}
