///////////////////////////////////////////////////////////////////////////
//
//  PEAK DETECTION OF FOURIER TRANSFORME MS INSTRUMENT DATA
//
//  written by Markus Mueller, markus.mueller@imsb.biol.ethz.ch
//  and Lukas Mueller, Lukas.Mueller@imsb.biol.ethz.ch
//  October 2005
//
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
//
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
// 
//

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/ms_peak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/CentroidData.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/LC_elution_peak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/LCMSCData.h>

namespace OpenMS
{

// minimal MS peak intensity allowed:
float LCMSCData::intensity_min_threshold;

// Constructor & destructor
///////////////////////////////////////////////////////////////////////////////

using namespace std;

LCMSCData::LCMSCData(){
  // Init();
}
///////////////////////////////////////////////////////////////////////////////

LCMSCData::~LCMSCData(){
  DATA.clear();
}

////////////////////////////////////////////////
// constructor for the object :
LCMSCData::LCMSCData(const LCMSCData* tmp){
  DATA = tmp->DATA;
}

////////////////////////////////////////////////
// constructor for the object :
LCMSCData::LCMSCData(const LCMSCData& tmp){
  DATA = tmp.DATA;
}

////////////////////////////////////////////////
// constructor for the object :
LCMSCData& LCMSCData::operator=(const LCMSCData& tmp){
  DATA = tmp.DATA;
  return *this;
}


// Public methods
///////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////
// find data of a specific m/z:
MZ_LIST_ITERATOR LCMSCData::get_MZ_by_iterator(double MZ){
  
  MZ_LIST_ITERATOR P = DATA.find(MZ);  
  
  return P;  
}


///////////////////////////////////////////////////////////////////////////////
// add data into the structure:
void LCMSCData::add_LC_elution_peak(double MZ, LC_elution_peak* IN){

  // get the scan apex:
  int APEX = IN->get_scan_apex();
  
  // check if this mass is already stored:
  MZ_LIST_ITERATOR P = get_MZ_by_iterator(MZ);
  
  if( P == get_DATA_end()){
    
    // ok, not yet inserted this mass:
    elution_peak_list tmp;
    tmp.insert(pair< int, LC_elution_peak>(APEX, *IN)); 
    IN = NULL;
    
    // insert into m/z list:
    DATA.insert(pair< double, elution_peak_list>( MZ, tmp));
    
    
  }
  else{
    
    // ok, already existing:
    (*P).second.insert(pair< int, LC_elution_peak>(APEX, *IN)); 
    IN = NULL;
  }
}


///////////////////////////////////////////////////////////////////////////////
// get a list of m/z observed in a scan +/- scan_tolerance:
// return the area of the LC elution peaks 
vector<LC_elution_peak> LCMSCData::get_MZ_list(int SCAN){
  
  int start_scan = SCAN;
  int end_scan = SCAN;
  
  vector<LC_elution_peak> OUT;  
  
  // go through the structure and find all m/z at this scan:
  MZ_LIST_ITERATOR P = get_DATA_start();
  
  while( P != get_DATA_end()){
    
    double this_INT = 0;
    LC_elution_peak* TMP = NULL;
    
    // search around some scan region:
    for( int this_scan = start_scan; this_scan < end_scan; this_scan++){
      
      // find elution peaks by the apex:
      elution_peak_list_ITERATOR Q = (*P).second.find( this_scan );
      
      if( Q != (*P).second.end()){
        double tmp = (*Q).second.get_total_peak_area();
        if(this_INT < tmp){
          this_INT = tmp;
          TMP = &((*Q).second);
        }
      }
    }
    
    if( ( this_INT > 0) && ( this_INT >= intensity_min_threshold) && (TMP != NULL) ){
      OUT.push_back(*TMP);
    }
    
    // next m/z:
    P++;
  }
  
  return OUT;
}


///////////////////////////////////////////////////////////////////////////////
// get a list of m/z observed in a scan +/- scan_tolerance:
// return the area of the LC elution peaks 
vector<LC_elution_peak> LCMSCData::get_MZ_list(int SCAN, int TOL){
  
  int start_scan = SCAN - TOL;
  int end_scan = SCAN + TOL;
  
  LC_elution_peak* TMP = NULL;
  vector<LC_elution_peak> OUT;  
  
  // go through the structure and find all m/z at this scan:
  MZ_LIST_ITERATOR P = get_DATA_start();
  
  while( P != get_DATA_end()){
    
    double this_INT = 0;
    
    // search around some scan region:
    for( int this_scan = start_scan; this_scan < end_scan; this_scan++){
      
      // find elution peaks by the apex:
      elution_peak_list_ITERATOR Q = (*P).second.find( this_scan );
      
      if( Q != (*P).second.end()){
        double tmp = (*Q).second.get_total_peak_area();
        if(this_INT < tmp){
          this_INT = tmp;
          TMP = &((*Q).second);
        }
      }
    }
    
    if( ( this_INT > 0) && ( this_INT >= intensity_min_threshold) && (TMP != NULL) ){
      OUT.push_back(*TMP);
    }
    
    // next m/z:
    P++;
  }
  
  return OUT;
}

///////////////////////////////////////////////////////////////////////////////
// get all extracted LC peaks:
vector<LC_elution_peak*> LCMSCData::get_ALL_peak(){
  
  LC_elution_peak* TMP = NULL;
  vector<LC_elution_peak*> OUT;  
  
  // go through the structure and find all m/z at this scan:
  MZ_LIST_ITERATOR P = get_DATA_start();
  
  while( P != get_DATA_end()){
        
    // find elution peaks by the apex:
    elution_peak_list_ITERATOR Q = (*P).second.begin();
    
    while( Q != (*P).second.end() ){
      TMP = &((*Q).second);
      // cout<<TMP->get_apex_MZ()<<endl;
      OUT.push_back(TMP);
      Q++;
    }
    
    // next m/z:
    P++;
  }
  
  return OUT;
}

//////////////////////////////////////////////////////
// prints profile of the extracted LC elution peaks:
void LCMSCData::print_profile_of_LC_peaks(ofstream* OUT){

  //////////////////////////////////
  // run through all m/z values:
//  MZ_LIST_ITERATOR P_MZ = get_DATA_start();
//  
//  while( P_MZ != get_DATA_end()){
//        
//    //////////////////////////////////
//    // run through all LC_elution_peaks
//    elution_peak_list_ITERATOR Q_SER = (*P_MZ).second.begin();
//    while( Q_SER != (*P_MZ).second.end()){
//      (*Q_SER).second.print_profile(OUT);
//      Q_SER++;
//    }
//    //////////////////////////////////
//    
//    P_MZ++;
//  }
//  OUT = NULL;
}


/////////////////////////////////////////////////////////////////////////////
// get a vector with all LC peaks ordered by their score:
vector<LC_elution_peak*> LCMSCData::get_ALL_peak_ordered(){
  vector<LC_elution_peak*> DATA = get_ALL_peak();
  return DATA;
}

}
