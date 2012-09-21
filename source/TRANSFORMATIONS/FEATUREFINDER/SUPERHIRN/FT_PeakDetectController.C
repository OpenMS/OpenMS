///////////////////////////////////////////////////////////////////////////
//
//  written by Lukas N Mueller, 30.3.05
//  Lukas.Mueller@imsb.biol.ethz.ch
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
//
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
//
// **********************************************************************//
// CLASS FT_PeakDetectController:
//
// variable description:
//
//
// function description:
//
//
// **********************************************************************//

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <stdio.h>

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
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/LC_MS.h>;
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS1_feature_merger.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/FT_PEAK_DETEC_mzXML_reader.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/FT_PeakDetectController.h>

namespace OpenMS
{

using namespace std;

bool FT_PeakDetectController::CREATE_FEATURE_ELUTION_PROFILES = false;
bool FT_PeakDetectController::LCelutionPeakDebugging = false;
double FT_PeakDetectController::LCelutionPeakMassMin = -1;
double FT_PeakDetectController::LCelutionPeakMassMax = -2;

MS2_feature* FT_PeakDetectController::SearchedM2Feature;


// if this option is on, then construct fake features for available MS2 features
bool FT_PeakDetectController::FEATURE_FAKE_INSERTION_BASED_ON_MS2_FEATURE = true;

////////////////////////////////////////////////
// constructor for the object FT_PeakDetectController:
FT_PeakDetectController::FT_PeakDetectController(){

 
  THIS_LC_MS = NULL;
  
}

//////////////////////////////////////////////////
// class desctructor of FT_PeakDetectController
FT_PeakDetectController::~FT_PeakDetectController(){
  
  LC_MS_RUNS.clear();
  if( THIS_LC_MS != NULL ){
    delete THIS_LC_MS;
    THIS_LC_MS = NULL;
  }
}

//////////////////////////////////////////////////
// class copy constructor of FT_PeakDetectController
FT_PeakDetectController::FT_PeakDetectController(const FT_PeakDetectController& tmp){}


//////////////////////////////////////////////////
// copy constructor:
FT_PeakDetectController& FT_PeakDetectController::operator=(const FT_PeakDetectController& tmp){
  return *this;
}
    


//////////////////////////////////////////////////
// start the scan parsing of a mzXML file
void FT_PeakDetectController::start_scan_parsing_of_mzXML_file(Vec datavec){
  
  // parse a mzXML file:
  printf("\n\t-- SuperHirn feature extraction of mzXML file '%s'\n", get_target_file().c_str() );
  
  FT_PEAK_DETEC_mzXML_reader* FT_READER = new FT_PEAK_DETEC_mzXML_reader();
  
  // set input file and open:
//  FT_READER->set_current_file( get_target_file() );
//  FT_READER->delete_existing_debug_file(); // old debug files are overwritten
//  FT_READER->open_mzxml_file();
  
  // get the titel of the current LC_MS run:
  string name = "tmplcms";
//  name.erase(0, name.rfind("/")+1);
//  name.erase(name.rfind(".mzXML"), name.size() - name.rfind(".mzXML"));
  
  // create a new LC/MS:
  THIS_LC_MS = new LC_MS( name );
  THIS_LC_MS->set_spectrum_ID( this->LC_MS_RUNS.size() );
  
  // start the mzXML reading process:
  FT_READER->read_mzXML_DATA(datavec);
    
  
  //////////////////////////////////////
  // post porocessing of mzXML data of a file:  
  // !!! MS1 LEVEL !!!
  process_MS1_level_data_structure( FT_READER );
  // !!! MS2 LEVEL !!!
  //process_MS2_level_data_structure( FT_READER );

  THIS_LC_MS->order_by_mass();
  
  if( MS1_feature_merger::MS1_FEATURE_CLUSTERING ){
    MS1_feature_merger* merg = new MS1_feature_merger( THIS_LC_MS );
    merg->startFeatureMerging();
    delete merg;
    merg = NULL;
  }
  
  THIS_LC_MS->show_info();
  
  
 /* Commented out for brutus */
  string SEP = "";
  FILE *file; 
  file = fopen("ffsh-features.txt","w+"); 
  fprintf(file,"%s", "Features\n");
  
  vector<feature>::iterator p = THIS_LC_MS->get_feature_list_begin();
  while(p != THIS_LC_MS->get_feature_list_end()){
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

  
  
  
  
  // store the extracted run to XML:
//  write_out_parsed_LC_MS( THIS_LC_MS );
  
  // add it to the spectrum list:
  LC_MS_RUNS.push_back( *THIS_LC_MS );
  
  // clear up:
  //delete THIS_LC_MS; 
  //THIS_LC_MS = NULL;
  
  delete FT_READER;
  FT_READER = NULL;
  
}

typedef multimap< double, MZ_series> MAIN_DATA_STRUCTURE;
typedef MAIN_DATA_STRUCTURE::iterator MAIN_ITERATOR;

////////////////////////////////////////////////////
// process MS1 level data
void FT_PeakDetectController::process_MS1_level_data_structure( FT_PEAK_DETEC_mzXML_reader* FT_READER){
  
  // get the processsed raw data structure back:
  // processing classes:  
  Process_Data* current_raw_mzXML_file =  FT_READER->get_processed_MS1_data_structure();
  
  // FLOFLO
  int mzListSize = current_raw_mzXML_file->getNbMSTraces();
  std::cout << "mzListSize: " << mzListSize << "\n";
    
  // extract LC elution features on the MS1 level:
  current_raw_mzXML_file->extract_elution_peaks();

  // get the new structure with the LC features:
  LCMSCData* current_processed_mzXML_file = current_raw_mzXML_file->get_processed_data();
  
  // iterator over the extracted features, convert
  vector<LC_elution_peak*> PEAKS = current_processed_mzXML_file->get_ALL_peak();
  // show program status:
  printf("\t\t\t* Processing of %d MS1 level features...\n",  PEAKS.size() );
  
  vector<LC_elution_peak*>::iterator P = PEAKS.begin();
  while( P != PEAKS.end()){    
    // add the LC peak for conversion to a feature structure
    // and to insert into the current LC-MS run
    add_raw_peak_to_LC_MS_run( (*P) );
    P++; 
  }
  
  // order the run again:
  THIS_LC_MS->order_by_mass();
  
  // clear:
  current_raw_mzXML_file = NULL;
  current_processed_mzXML_file = NULL;
  FT_READER = NULL;
}



/* FLO
////////////////////////////////////////////////////
// process MS2 level data
void FT_PeakDetectController::process_MS2_level_data_structure( FT_PEAK_DETEC_mzXML_reader* FT_READER){

  
  // get the processsed raw data structure back:
  MS2_Process_Data* current_MS2_raw_mzXML_file =  FT_READER->get_processed_MS2_data_structure();
    
  if( current_MS2_raw_mzXML_file->getNbMSTraces() > 0 ){
    
    // show program status:
    // printf("\t\t\t* Processing of %d MS2 Fragments Traces...\n", current_MS2_raw_mzXML_file->getNbMSTraces() );
    
    // convert extracted traces:
    current_MS2_raw_mzXML_file->constructMS2ConsensusSpectra( );
    
    // show how many consensus spectrum extracted:
    printf("\t\t\t* Processing of %d MS2 level features...\n", current_MS2_raw_mzXML_file->getNbMS2Features() );
      
    //////////////////////////////////////////
    // now get the list of MS2_features
    // and associate them to MS1 features:
    progress_bar bar(current_MS2_raw_mzXML_file->getNbMS2Features(),"");
    vector< MS2_feature >::iterator F = current_MS2_raw_mzXML_file->getMS2FeaturesStart();
    while( F != current_MS2_raw_mzXML_file->getMS2FeaturesEnd() ){
      associateMS2FeatureToMS1Feature( &(*F) );
      F++;
      bar.update_progress();
    }
    
    
    printf("\n\t\t\t\t- %d MS2 features added to MS1 features\n", current_MS2_raw_mzXML_file->getNbMS2Features() - fakeFeatureList.size());
     
   }
  
  // clear:
  current_MS2_raw_mzXML_file = NULL;
  FT_READER = NULL;
}


////////////////////////////////////////////////////
// combine the MS2 feature trace data to the MS1 features:
void FT_PeakDetectController::associateMS2FeatureToMS1Feature( MS2_feature* ms2){
  
  FT_PeakDetectController::SearchedM2Feature = ms2;
  vector<LC_MS_FEATURE>::iterator F = THIS_LC_MS->get_feature_list_end();
  if( MS2_Process_Data::POST_MSMS_SPECTRUM_CLUSTERING ){
    F = find_if( THIS_LC_MS->get_feature_list_begin(), THIS_LC_MS->get_feature_list_end(), MS2ToMS1Comparer() );
    if( F != THIS_LC_MS->get_feature_list_end()  ){
      addMS2FeatureToMS1Feature( ms2, &(*F) );
    }
  }
  
  // construct here fake ms1 features based on a observed MS2 feature
  // which however could not be matched to a exiting ms1 feature
  if( ( F == THIS_LC_MS->get_feature_list_end()) && FEATURE_FAKE_INSERTION_BASED_ON_MS2_FEATURE ){
    constructMS1FeatureFromMS2Feature( ms2 );
    
  }
}
 */

////////////////////////////////////////////////////
// add an observed MS2 feature to the MS1 feature
// if an observation is already there, then
// construct a merged MS2 feature
void FT_PeakDetectController::addMS2FeatureToMS1Feature( MS2_feature* ms2, feature* ms1 ){
  
  if( ms1->getMS2Feature() == NULL ){
    ms1->addMS2Feature( ms2 );  
  }
  else{
    
    MS2_feature* previousMS2 = ms1->getMS2Feature();
    previousMS2->addMS2ConsensusSpectrum( ms2 );
  
    // in the case where a MS1 feature was generated from a MS1, 
    // then adjust the retention time levels 
    if( ms1->get_peak_area() == -1 ){
      
      // start:
      if( ms2->getStartTR() < ms1->get_retention_time_START() ){
        ms1->set_retention_time_START( ms2->getStartTR() );
      }
      // end:
      if( ms2->getEndTR() > ms1->get_retention_time_END() ){
        ms1->set_retention_time_END( ms2->getEndTR() );
      }
      
    }
  }
}


////////////////////////////////////////////////////
// construct here fake ms1 features based on a observed MS2 feature
// which however could not be matched to a exiting ms1 feature
void FT_PeakDetectController::constructMS1FeatureFromMS2Feature( MS2_feature* in ){

  feature* fakeMS1 = new feature( in );    
  THIS_LC_MS->add_feature( fakeMS1 );
  delete fakeMS1; 
  fakeMS1 = NULL;
  in = NULL;
}
  

////////////////////////////////////////////////////
// adds an elution peak to the LC/MS run:
void FT_PeakDetectController::add_raw_peak_to_LC_MS_run( LC_elution_peak* PEAK ){
  
  
  ////////////////////////////////
  // get parameters of the peak:
  // apex:
  int apex_scan = PEAK->get_scan_apex();
  double apex_MZ = PEAK->get_apex_MZ();
  double apex_TR = PEAK->get_apex_retention_time();
  float apex_INTENSITY = PEAK->get_apex_intensity();
  
  // peak shape:
  float peak_area = PEAK->get_total_peak_area();
  int charge_state = PEAK->get_charge_state();
  int peak_start = PEAK->get_start_scan();
  int peak_end = PEAK->get_end_scan();
  
  if( (apex_TR <= FT_PEAK_DETEC_mzXML_reader::TR_MAX) && (apex_TR >= FT_PEAK_DETEC_mzXML_reader::TR_MIN) ){
    
    ////////////////////////////////
    // in case of debugging mode:
    // print out the LC peak info:

    
    ////////////////////////////////////
    // construct a feature:
    feature* TMP = new feature(apex_MZ, apex_TR, apex_scan, peak_start,peak_end, charge_state, peak_area, apex_INTENSITY, 0);
    
    // set start / end TR:
    TMP->set_retention_time_START(PEAK->get_start_retention_time());
    TMP->set_retention_time_END(PEAK->get_end_retention_time());

    // set ids:
    TMP->set_spectrum_ID( THIS_LC_MS->get_spectrum_ID() );
    TMP->set_feature_ID( THIS_LC_MS->get_nb_features() );
    
    // set S/N
    TMP->setSignalToNoise( PEAK->getSignalToNoise() );
    // set background noise:
    TMP->setBackgroundNoiseLevel( PEAK->getSignalToNoiseBackground() );
    
    // set feature extract information, if available:
    if( ! PEAK->getElutionPeakExtraInfo( ).empty()){
      TMP->setFeatureExtraInformation( PEAK->getElutionPeakExtraInfo( ) );
      // add fake MS/MS information for the MS1 feature:
      addFakeMSMSToFeature( TMP );
    }
    

    // function to add the elution profile to the feature:
    // (if turned on)
    if( CREATE_FEATURE_ELUTION_PROFILES ){
      addLCelutionProfile( TMP, PEAK );
    }

    // add it to the LC_MS run:
    THIS_LC_MS->add_feature(TMP);
    TMP = NULL;
  
  }
  PEAK = NULL;
}

////////////////////////////////////////////////////////
// add fake MS/MS information for the MS1 feature:
void FT_PeakDetectController::addFakeMSMSToFeature( feature* in ){
  
  string tmp = in->getFeatureExtraInformation();
  string tag = "INFO:";
  string sep = ";";
  // parse name etc out:
  tmp = tmp.substr( tmp.find( tag ) + tag.size() ); 
  
  string AC = tmp.substr( 0, tmp.find( sep )); 
  tmp = tmp.substr( tmp.find( sep ) + sep.size() ); 
  
  string SQ = tmp.substr( 0, tmp.find( sep )); 
  tmp = tmp.substr( tmp.find( sep ) + sep.size() ); 

  // int id = atoi( (tmp.substr( 0, tmp.find( sep )) ).c_str() ); 
  
  ms2_info* info = new ms2_info( AC, SQ, in->get_charge_state(), 1.0);
  info->set_MONO_MZ( in->get_MZ() );

  info->set_SCAN_START( in->get_scan_number() );
  info->set_SCAN_END( in->get_scan_number() );
  info->setRetentionTime( in->get_retention_time());
  info->set_PREV_AA( "R/K" );
  
  in->add_MS2_info( info );

  delete info;
  info = NULL;
}



////////////////////////////////////////////////////////
// function to add the elution profile to the feature:
void FT_PeakDetectController::addLCelutionProfile( feature* inF, LC_elution_peak* PEAK ){

  ////////////////////////////////
  // get parameters of the peak:
  // apex:
  int apex_scan = PEAK->get_scan_apex();
  double apex_MZ = PEAK->get_apex_MZ();
  double apex_TR = PEAK->get_apex_retention_time();
  float apex_Intensity = PEAK->get_apex_intensity();

  // peak shape:
  float peak_area = PEAK->get_total_peak_area();
  int charge_state = PEAK->get_charge_state();
  
  // create the class:
  featureLCprofile* myProfile = new featureLCprofile(apex_MZ, apex_TR, apex_Intensity, apex_scan, charge_state, peak_area);

  // add the raw individual mono isotopic peaks peaks:
  SIGNAL_iterator P = PEAK->get_signal_list_start();
  while( P != PEAK->get_signal_list_end() ){

    ms_peak* MSpeak = &(*P).second;
    myProfile->addMS1elutionSignal( MSpeak->get_MZ(), MSpeak->get_intensity(), MSpeak->get_scan_number(), MSpeak->get_charge_state(), MSpeak->get_retention_time() );    
    P++;
  
  }

  inF->setLCelutionProfile( myProfile );
  
  myProfile = NULL;
}

}
