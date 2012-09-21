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

#include <map>
#include <vector>
#include <string>
#include <list>
#include <iostream>
#include <stdio.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/ms_peak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/CentroidPeak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/CentroidData.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/LC_elution_peak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/BackgroundIntensityBin.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/BackgroundControl.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/Deisotoper.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/LCMSCData.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/simple_math2.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/LC_MS_XML_reader.h>;

// TODO Refactor
// The four includes MS2_feature, featureLCprofile, ms2_info and feature are just 
// there for reading one constant. Get rid of them.
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2_feature.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/featureLCprofile.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/ms2_info.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/feature.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/Process_Data.h>

namespace OpenMS
{

using namespace std;

map<int, float> Process_Data::scan_TR_index;

// to track detected monoistopic mass for debugging:
bool Process_Data::MonoIsoDebugging = false;
double Process_Data::DebugMonoIsoMassMin =318.00; 
double Process_Data::DebugMonoIsoMassMax =  319.00; 

double Process_Data::MS1_intensity_apex_percentil_cutoff;
double Process_Data::MS1_TR_RESOLUTION;
// if data are in centroid form or not:
bool Process_Data::CENTROID_DATA_MODUS = false;



double Process_Data::MZ_TOLERANCE = 10; // overwritten by config
float Process_Data::INTENSITY_THRESHOLD = 0; // overwritten by config
double Process_Data::max_inter_scan_retention_time_distance = 0; // overwritten by config
int Process_Data::min_nb_cluster_members = 0;


////////////////////////////////////////////////
// constructor for the object Process_Data:
Process_Data::Process_Data(){
  
  DATA = new LCMSCData();
  LC_elution_peak_counter = 0;
  
  //////////////////
  // parameter to visualize
  // the data & peak detection:
  // read parameters from "param.def"  
  //read_param* def = new read_param();
  
  //minimal intensity of a ms peak::
  //def->search_tag("FT peak detect MS1 intensity min threshold", &INTENSITY_THRESHOLD);
  
  // m/z tolerance value:
  // def->search_tag("FT peak detect MS1 m/z tolerance", &MZ_TOLERANCE);
 
  // minimal number of cluster members
  //def->search_tag("FT peak detect MS1 min nb peak members", &min_nb_cluster_members);
  //def->search_tag("MS1 max inter scan distance", &max_inter_scan_retention_time_distance );
  TIME_CLUSTERING_BY_RETENTION_TIME = true;
 
  //delete def;
  //def = NULL;
  
  backgroundController = new BackgroundControl();

}

//////////////////////////////////////////////////
// class desctructor of Process_Data
Process_Data::~Process_Data(){
  // empty the raw data:
  pMZ_LIST.clear();
  if(DATA != NULL){
    delete DATA;
    DATA = NULL;
  }
  
  if( backgroundController != NULL ){
    delete backgroundController;
    backgroundController = NULL;
  }
}

//////////////////////////////////////////////////
// class copy constructor of Process_Data
//Process_Data::Process_Data(const Process_Data& tmp){
//}


//////////////////////////////////////////////////
// copy constructor:
//Process_Data& Process_Data::operator=(const Process_Data& tmp){
//  return *this;
//}

///////////////////////////////////////////////////////////////////////////////
// get an observed MZ mass, otherwise end of list iterator
MAIN_ITERATOR Process_Data::get_MZ(double IN_mz){
  return pMZ_LIST.find(IN_mz);
}

///////////////////////////////////////////////////////////////////////////////
// find closest match mz mass in the main structure
MAIN_ITERATOR Process_Data::find_closest_mz_match(double MZ){
  
  MAIN_ITERATOR P = pMZ_LIST.lower_bound(MZ);
  
  if( MZ == (*P).first)
    return P;
  
  double inf = 10000000;
  
  MAIN_ITERATOR P_UP;
  double up = inf;
  MAIN_ITERATOR P_DOWN;
  double down = inf;
  
  if( P != get_MZ_LIST_end()){
    P_UP = P;
    up = fabs((*P_UP).first - MZ);
  }
  
  if( P != get_MZ_LIST_start()){
    P--;
    P_DOWN = P;
    down = fabs(MZ - (*P_DOWN).first);
  }
  
  if(down < up){
    if(down > simple_math2::getMassErrorAtPPMLevel(MZ, MZ_TOLERANCE) )
      printf("\nERROR: no tolerance-match found, eventhough should!!!!\n"); 
    return P_DOWN;
  }
  else{
    if(up > simple_math2::getMassErrorAtPPMLevel(MZ, MZ_TOLERANCE))
      printf("\nERROR: no tolerance-match found, eventhough should!!!!\n"); 
    return P_UP;
  }
  
  
}

///////////////////////////////////////////////////////////////////////////////
// get an observed MZ mass, otherwise end of list iterator
MAIN_ITERATOR Process_Data::get_MZ_lower_bound(double IN_mz){
  return pMZ_LIST.lower_bound(IN_mz);
}

///////////////////////////////////////////////////////////////////////////////
// get end of MZ list:
MAIN_ITERATOR Process_Data::get_MZ_LIST_end(){
  return pMZ_LIST.end();
}
///////////////////////////////////////////////////////////////////////////////
// get start of MZ list:
MAIN_ITERATOR Process_Data::get_MZ_LIST_start(){
  return pMZ_LIST.begin();
}

///////////////////////////////////////////////////////////////////////////////
// erase element in MZ list:
void Process_Data::erase_MZ_LIST_element(MAIN_ITERATOR IN){
  if(IN == pMZ_LIST.end()){
    printf("\nERROR: could not erase end iterator, Process_Data::erase_MZ_LIST_element()!!!!\n");           
  }
  pMZ_LIST.erase(IN);
}

///////////////////////////////////////////////////////////////////////////////
// find element numbers:
map<double, int>::iterator Process_Data::get_nb_MZ_cluster_elements(double IN){

  // map<double, int>::iterator OUT = MZ_CLUSTER.lower_bound(IN);
  // IN = simple_math::ROUND_NUMBER(IN,3);
  map<double, int>::iterator OUT = MZ_CLUSTER.find( IN );
  if(IN == (*OUT).first){
    return OUT;
  }
  
  printf("\nERROR: no match in MZ_CLUSTER found, Process_Data::get_nb_MZ_cluster_elements(double)!!!!\n");       
  return MZ_CLUSTER.end();
  
}


///////////////////////////////////////////////////////////////////////////////
// get the  full summed up intensity
double Process_Data::getPeakIntensitySum(double IN){
  
  double out = 0;
  MAIN_ITERATOR F = pMZ_LIST.find(IN);
  if(F != pMZ_LIST.end()){
    
    MZ_series_ITERATOR p = F->second.begin();
    while( p != F->second.end() ){
      multimap<int, ms_peak >::iterator k = p->begin();
      while( k != p->end() ){
        out += k->second.get_intensity();
        k++;
      }
      p++;
    }
    return out;
  }
  
  printf("\nERROR: no match in MZ_CLUSTER found, Process_Data::getMzAverageAndIntensitySum(double)!!!!\n");       
  return out;
  
}


///////////////////////////////////////////////////////////////////////////////
// erase an element:
void Process_Data::erase_MZ_cluster_element(map<double, int>::iterator IN){
  if(IN == MZ_CLUSTER.end()){
    printf("\nERROR: could not erase end iterator, Process_Data::erase_MZ_cluster_element()!!!!\n");           
  }
  MZ_CLUSTER.erase(IN);
}

///////////////////////////////////////////////////////////////////////////////
// erase an element:
void Process_Data::insert_MZ_cluster_element(double IN, int NB){
  // IN = simple_math::ROUND_NUMBER(IN,3);
  MZ_CLUSTER.insert(make_pair(IN,NB));
}


///////////////////////////////////////////////////////////////////////////////
// find a retention time by the scan number:
double Process_Data::find_retention_time(double IN){

  if(scan_TR_index.size() > 0){
    
    int SCAN = int(ceil(IN));
    map<int, float>::iterator P = scan_TR_index.lower_bound(SCAN);
    
    if( P == scan_TR_index.end()){
      P--;
      return (*P).second;
    }
    
    if( (*P).first == IN){
      return (*P).second;
    }
    else{
      
      double TR_up = (*P).second;
      if(P != scan_TR_index.begin()){
        double SCAN_up = double((*P).first);
        P--;
        double TR_down = (*P).second;
        double SCAN_down = double((*P).first);
        double w_up = (SCAN_up - SCAN_down) / (SCAN_up - IN);
        double w_down = (SCAN_up - SCAN_down) / (IN - SCAN_down);
        return ( w_up * TR_up + w_down*TR_down ) / (w_up + w_down);
      }
      else{
        return TR_up;
      }
    }
  }
  else{
    return 0;
  }
}

///////////////////////////////////////////////////////////////////////////////
// inputs raw /centroided  data into the object:
void Process_Data::add_scan_raw_data(int SCAN, float TR, CentroidData* centroidedData){
  
  Deisotoper dei; 
  
  //////////////////////////////////
  // add the peaks to the background controller:
  list<CentroidPeak> pCentroidPeaks;
  centroidedData->get(pCentroidPeaks);
  backgroundController->addPeakMSScan( TR, &pCentroidPeaks);
  
  dei.go(*centroidedData);
  dei.cleanDeconvPeaks();

  // convert to objects used for mass clustering over retention time
  vector<ms_peak> PEAK_LIST; 
  convert_ms_peaks(SCAN, TR, dei.getDeconvPeaks(), PEAK_LIST);
  
  //  store it:
  this->add_scan_raw_data(PEAK_LIST);
  
  // clear it:
  PEAK_LIST.clear();

  
}
 

///////////////////////////////////////////////////////////////////////////////
// inputs the centroided / deisotoped data into the object:
void Process_Data::add_scan_raw_data(vector<ms_peak> PEAK_LIST){
  
  // add the peaks to the background controller:
  //backgroundController->addPeakMSScan( TR, &PEAK_LIST );

  // iterate through the vector:
  vector<ms_peak>::iterator P = PEAK_LIST.begin();
  while(P != PEAK_LIST.end()){
    
    ms_peak* PEAK = &(*P);
        
    // check if its above the min. intensity:
    if( filterDeisotopicMSPeak( PEAK ) ){
            
      // set the retention time of the ms peak:
      // PEAK->set_retention_time(TR);
    
      // check if this MZ has already been observed:
      MAIN_ITERATOR LCP = check_MZ_occurence( PEAK );
      if(  LCP != get_MZ_LIST_end() ){
        insert_observed_mz(LCP, PEAK);
      }
      else{
        insert_new_observed_mz( PEAK);
      }
    }
    
    PEAK = NULL;    
    P++;
  }
  
}

///////////////////////////////////////////////////////////////////////////////
// check if the ms peak is in the selected mz, z, int range
bool Process_Data::filterDeisotopicMSPeak(ms_peak* PEAK){
  
  //////////////////////////
  // check if its above the min. intensity:
  if( PEAK->get_intensity() < getMinimalIntensityLevel() ){
    return false;
  }
  
  //////////////////////////
  // m/z selection range:
  if( ( PEAK->get_MZ() + simple_math2::getMassErrorAtPPMLevel( PEAK->get_MZ() , feature::PPM_MZ_TOL) < LC_MS_XML_reader::FEATURE_MZ_MIN ) || ( PEAK->get_MZ() - simple_math2::getMassErrorAtPPMLevel( PEAK->get_MZ(), feature::PPM_MZ_TOL ) > LC_MS_XML_reader::FEATURE_MZ_MAX ) ){
    return false;
  }
  
  //////////////////////////
  // charge state selection range:
  if( ( PEAK->get_Chrg() < LC_MS_XML_reader::FEATURE_CHRG_MIN ) || ( PEAK->get_Chrg() > LC_MS_XML_reader::FEATURE_CHRG_MAX ) ){
    return false;
  }

  return true;
}

///////////////////////////////////////////////////////////////////////////////
// insert a newly observed mz into the data structure
void Process_Data::insert_new_observed_mz( ms_peak* PEAK ){
  
  // DEBUGGING
  if( MonoIsoDebugging ){
    if( ( DebugMonoIsoMassMin <= PEAK->get_MZ()) && ( DebugMonoIsoMassMax >= PEAK->get_MZ()) ){
      cout<<endl<<"****"<<endl<<"->Insert new mz cluster from: ";
      PEAK->show_info();
    }
  }
  
  // create first an elution peak:
  elution_peak tmp_TR;
  tmp_TR.insert( make_pair( PEAK->get_Scan(), *PEAK) );
  
  // now make a vector for the mz:
  MZ_series tmp_MZ;
  tmp_MZ.push_back( tmp_TR );

  // into main structure:
  pMZ_LIST.insert( make_pair( PEAK->get_MZ(), tmp_MZ ) );
  
  // insert the mz cluster mean:  
  // insert_MZ_cluster_element( PEAK->get_MZ(), 1 );
  
  // increase the LC_elution_profile counter:
  increase_LC_elution_peak_counter();
  
  PEAK = NULL;
  
}

///////////////////////////////////////////////////////////////////////////////
// insert an already observed mz into the data structure, checks
// if it belongs to an existing LC elution peak or starts a new one:
void Process_Data::insert_observed_mz(MAIN_ITERATOR LCP, ms_peak* PEAK){
  
  ////////////////////////////////////////////////////
  // check if its the same m/z and charge state:
  if( ((*LCP).first == PEAK->get_MZ())  ){
    
    // find the last elution peak cluster: 
    MZ_series_ITERATOR Q = (*LCP).second.end();
    Q--;
      
    // check if this peak should be added to the existing 
    // last elution peak cluster or start a new one:
    if( check_elution_peak_belong(Q, PEAK) ){
      
      // add to this cluster the ms peak:
      (*Q).insert(pair<int, ms_peak>(PEAK->get_Scan(), *PEAK));
      
      // DEBUGGING
      if( MonoIsoDebugging ){
        if( ( DebugMonoIsoMassMin <= PEAK->get_MZ()) && ( DebugMonoIsoMassMax >= PEAK->get_MZ()) ){
          printf("%0.3f: ",(*LCP).first);
          PEAK->show_info();
        }
      }
    
    }
    else{
      
      // add a new one:
      elution_peak tmp_TR;
      tmp_TR.insert( make_pair( PEAK->get_Scan(), *PEAK ) );
      (*LCP).second.push_back(tmp_TR);      
      
      
      // DEBUGGING
      if( MonoIsoDebugging ){
        if( ( DebugMonoIsoMassMin <= PEAK->get_MZ()) && ( DebugMonoIsoMassMax >= PEAK->get_MZ()) ){
          multimap<int, ms_peak>::reverse_iterator q = Q->rbegin();
          int last_scan = (*q).first;
          printf("\n----\n-> Old mz %0.3f to %d, but new TR cluster: ", (*LCP).first, last_scan);
          PEAK->show_info();
          
        }
      }      

      // increase the LC_elution_profile counter:
      increase_LC_elution_peak_counter();
    }    
  }    
  ////////////////////////////////////////////////////
  
  
  else{
    
    // ok, we have the correct m/z:
    double match_mz = (*LCP).first;
    
    // calculate the average mass, get # of observed in the m/z cluster:
    double nb_elements = 1;
    nb_elements = LCP->second.rbegin()->size();    
    
    // calculate the new cluster average mass:
    double peakIntens = getPeakIntensitySum( match_mz );
    double new_mz = peakIntens* match_mz + PEAK->get_MZ()*PEAK->get_intensity();
    new_mz /= ( peakIntens + PEAK->get_intensity() );
    

    //////////////////////////////////
    // erase old and a add new cluster mz:
    // erase_MZ_cluster_element(W);
    // insert_MZ_cluster_element(new_mz, nb_elements);
    
    // now replace the value of the old MZ_SERIES with the new m/z value 
    // and add the input ms peak:
    MZ_series TMP_SER = LCP->second;
    erase_MZ_LIST_element( LCP );
 
    // find the last elution peak cluster: 
    MZ_series_ITERATOR Q = TMP_SER.end();
    Q--;    
    
    // check if this peak should be added to the existing 
    // last elution peak cluster or start a new one:
    if( check_elution_peak_belong(Q, PEAK) ){
      // DEBUGGING
//      if( MonoIsoDebugging ){
//        if( ( DebugMonoIsoMassMin <= PEAK->get_MZ()) && ( DebugMonoIsoMassMax >= PEAK->get_MZ()) ){
//          // get the last MS peak
//          ms_peak* last_peak = &(Q->rbegin()->second);
//          printf("%0.3f (%d):\n", match_mz, (*Q).size());
//          last_peak->show_info();
//          PEAK->show_info();
//          //cout<<endl;
//        }
//      }
      
      // add to this cluster the ms peak:
      (*Q).insert(pair<int, ms_peak>(PEAK->get_Scan(), *PEAK));
      // add it to the mZ cluster:
      pMZ_LIST.insert(pair<double, MZ_series>(new_mz, TMP_SER));
    
    }
    else{
     
      
      // DEBUGGING
      if( MonoIsoDebugging ){
        if( ( DebugMonoIsoMassMin <= PEAK->get_MZ()) && ( DebugMonoIsoMassMax >= PEAK->get_MZ()) ){
          multimap<int, ms_peak>::reverse_iterator q = Q->rbegin();
          int last_scan = (*q).first;
          printf("\n----\n-> Old mz %0.3f to %d, but new TR cluster: ", (*LCP).first, last_scan);
          PEAK->show_info();
        }
      }      
      
      // add a new one:
      elution_peak tmp_TR;
      tmp_TR.insert(pair<int, ms_peak>(PEAK->get_Scan(), *PEAK));
      TMP_SER.push_back(tmp_TR);
      // into main structure:
      pMZ_LIST.insert(make_pair(new_mz, TMP_SER));

      // increase the LC_elution_profile counter:
      increase_LC_elution_peak_counter();
      
    }
  }

  PEAK = NULL;
}


///////////////////////////////////////////////////////////////////////////////
// check if a peak with this scan number belong to this elution cluster:
bool Process_Data::check_elution_peak_belong(MZ_series_ITERATOR P, ms_peak* PEAK){
  
    
  
  // get the last element:
  multimap<int, ms_peak>::reverse_iterator q = P->rbegin();
  int last_scan = (*q).first;
  ms_peak* last_peak = &(q->second);

  // avoid clustering of masses within the same scan:
  if( PEAK->get_Scan() == last_peak->get_Scan()){
    return false;
  }
  
  
  
  if( ! TIME_CLUSTERING_BY_RETENTION_TIME ){

    // compare the scan numbers:
    if( (PEAK->get_Scan() - last_scan) <= getMaxScanDistance() ){
      return true;
    }

  }
  else{
    // compare the scan numbers:
    double deltaTr = PEAK->get_retention_time() - last_peak->get_retention_time();
    if( deltaTr <= max_inter_scan_retention_time_distance ){
      return true;
    }
  }

  return false;
}


///////////////////////////////////////////////////////////////////////////////
// returns the distance to this elution peak:
int Process_Data::getElutionPeakDistance(MZ_series_ITERATOR P, int SCAN){
  
  // get the last element:
  multimap<int, ms_peak>::reverse_iterator q = P->rbegin();
  int last_scan = (*q).first;
  
  return (SCAN - last_scan);
}



///////////////////////////////////////////////////////////////////////////////
// runs through the whole data structure and puts the elution_peaks into
// a proper LC_elution peak object
void Process_Data::extract_elution_peaks(){
  
  backgroundController->processIntensityMaps();
  //backgroundController->plotIntensityMaps();
  // backgroundController->writeIntensityMaps();
  
  double this_MZ = 0;
  // progress_bar bar(get_LC_elution_peak_counter(),"processed");
  
  //////////////////////////////////
  // run through all m/z values:
  MAIN_ITERATOR P_MZ = get_MZ_LIST_start();
  while(P_MZ != get_MZ_LIST_end()){
    
    this_MZ = (*P_MZ).first;
    
    //////////////////////////////////
    // run through all elution peaks
    MZ_series_ITERATOR Q_SER = (*P_MZ).second.begin();
    while( Q_SER != (*P_MZ).second.end()){
      
      // check if this elution peak
      // is accepted as a really LC-elution peak:
      if(check_elution_peak(Q_SER)){
        convert_to_LC_elution_peak(Q_SER, this_MZ);    
      }
      
      Q_SER++;
      
    }
    //////////////////////////////////
    
    P_MZ++;
  }
  
}

///////////////////////////////////////////////////////////////////////////////
// check if this elution peak is accepted as a really LC-elution peak:
bool Process_Data::check_elution_peak(MZ_series_ITERATOR Q_SER){
  
  
  
  if( MonoIsoDebugging ){
    
    ms_peak* PEAK = &((Q_SER->begin())->second);
    
    if( ( DebugMonoIsoMassMin <= PEAK->get_MZ()) && ( DebugMonoIsoMassMax >= PEAK->get_MZ()) ){
      multimap<int, ms_peak>::iterator start = Q_SER->begin();
      multimap<int, ms_peak>::reverse_iterator end = Q_SER->rbegin();
      printf("\n : this peak removed: ");
      start->second.show_info();
      end->second.show_info();
      
    }
  }      
  
  /*
   // define here to how long the LC elution profile should
   // be depending on the ms_peak score:
   double score = 0;
   elution_peak::iterator P = (*Q_SER).begin();
   while( P != (*Q_SER).end() ){
     score += (*P).second.get_score();
     P++;
   }
   score /= (*Q_SER).size();
   int min_cluster_length = int( double(min_nb_cluster_members-1) * (5.0 - score) + 1.0);
   
   */
  
   
  // check if contains more or same than x element: 
  if( int((*Q_SER).size()) >= min_nb_cluster_members){
    return true;
  }
  else{
    
    ////////////////////////////////////////////////
    // check if the peak contains any MS/MS peak which was
    // selected as MS/MS precursor:
    multimap<int, ms_peak >::iterator P = (*Q_SER).begin();
    while( P != (*Q_SER).end() ){
      
      if( P->second.getPrecursorActivation() ){
        return true;
        //P->second.show_info();
      }
      
      P++;
    }
    
    
          
    return false;
  }
}

///////////////////////////////////////////////////////////////////////////////
// convert the MZ_series elution peak element into a LC_elution_peak object
void Process_Data::convert_to_LC_elution_peak(MZ_series_ITERATOR Q_SER, double this_MZ){
    
  // set important processing parameters such as apec cutoff for noise removal
  // and the tr delta steps for peak area integration
  // LC_elution_peak::intensity_apex_percentil_cutoff = Process_Data::MS1_intensity_apex_percentil_cutoff;
  LC_elution_peak::TR_RESOLUTION = Process_Data::MS1_TR_RESOLUTION;
  
  // process MS peaks before addingL:
  processMSPeaks( &(*Q_SER) );
  
  // create the object:
  LC_elution_peak* tmp = new LC_elution_peak(Q_SER, this_MZ); 
  // analyze the peak, i.e. compute parameters
  tmp->analyzeLCElutionPeak();
  
  // save in a new data structure
  DATA->add_LC_elution_peak(this_MZ, tmp);
  
  tmp = NULL;
}


///////////////////////////////////////////////////////////////////////////////
// process a series of MS peaks
// set the signal to noise level:
void Process_Data::processMSPeaks(multimap<int, ms_peak>* in){
  
  multimap<int, ms_peak>::iterator I = in->begin();
  while( I != in->end() ){
    
    ms_peak* peak = &(I->second);
    double bgLevel = backgroundController->getBackgroundLevel( peak->get_MZ(), peak->get_retention_time() );
    double SN = peak->get_intensity() / bgLevel;
    
    /*
    double avSN = 0;
    double totIntens = 0;
    vector<CentroidPeak>::iterator c = peak->get_isotopic_peaks_start();
    while( c != peak->get_isotopic_peaks_end() ){
    
      double SN = c->getOrgIntensity() / bgLevel;
      c->setSignalToNoise(SN);
      avSN += SN *  c->getOrgIntensity();
      totIntens += c->getOrgIntensity();
      c++;
    }
    */
    
    //avSN /= totIntens;
    peak->setSignalToNoise( SN );
    
    I++;
  }

}  

///////////////////////////////////////////////////////////////////////////////
// checks if a mz value has already been seen,
// also look for very close ones and cluster them
MAIN_ITERATOR Process_Data::check_MZ_occurence(ms_peak* PEAK){
  
  if( MonoIsoDebugging ){
    if( ( DebugMonoIsoMassMin <= PEAK->get_MZ()) && ( DebugMonoIsoMassMax >= PEAK->get_MZ()) ){
      if( PEAK->get_Scan() == 6289 ){
        PEAK->show_info();
      }
    }
  }
    
  ////////////////////////////////
  // check first for the possible candidates within the MZ tolerance range:
  double targetMZ = PEAK->get_MZ();
  int targetScan = PEAK->get_Scan();
  
  MAIN_ITERATOR  P = get_MZ_lower_bound(targetMZ); 
  vector<MAIN_ITERATOR> CandidateList;
  
  // go decreasing order
  MAIN_ITERATOR  down = P;
  if( P != get_MZ_LIST_start() ){
    do{ 
      down--;
      int check = compareIteratorToPeak( PEAK, down );
      if( check == 1 ){
        CandidateList.push_back( down );
      }
      else if( check == -1){
        break;
      }
      
      
    }
    while( down != get_MZ_LIST_start() );
  }
  
  
  // go increasing order
  MAIN_ITERATOR up = P;
  while( up != get_MZ_LIST_end() ){
    
    int check = compareIteratorToPeak( PEAK, up );
    if( check == 1 ){
      CandidateList.push_back( up );
    }
    else if( check == -1){
      break;
    }
    
    
    up++;
  }
  
  // here the list of possible candidates for 
  // matching this mass is complete
  // now -> find the best one according o
  // a: closest in m/z
  
  if( CandidateList.empty() ){
    // if its not found at all:
    // return the end:
    P = get_MZ_LIST_end();
  }
  else if( CandidateList.size() == 1 ){
    P = *(CandidateList.begin());
  }
  else{
    
    ///////////
    // check for those with smallest mz:
    // (within the tolerance range:)
    P = get_MZ_LIST_end();
    double smallMZDiff = 1000000;
    int smallScanDiff = 1000000;
    vector<MAIN_ITERATOR>::iterator Z = CandidateList.begin();
    while( Z != CandidateList.end() ){
      
      // get the inter scan distance:
      double MZDiff = fabs(targetMZ - (*Z)->first); 
      // get the difference in the scan numbers:
      MZ_series_ITERATOR x = (*Z)->second.end();
      x--;
      int ScanDiff = getElutionPeakDistance( x, targetScan );

      // store the smallest mz diff:
      if( ( MZDiff < smallMZDiff ) && ( ScanDiff < smallScanDiff )){
        P = *Z;
        smallMZDiff = MZDiff;
      }
      
      // scan difference:
      if( ( ScanDiff < smallScanDiff ) && ( ScanDiff <= getMaxScanDistance()) ) {
        P = *Z;
        smallScanDiff = ScanDiff;
      }
      
      Z++;
      
    }
    
  }
  
  return P;
}




//////////////////////////////////////////////////////////////////
// function which check if a data structure iterator is similar
// to a peak and should be considered 
// returns 1 if ok
// returns 0 if not
// returns -1 if scan range exceeded
int Process_Data::compareIteratorToPeak(ms_peak* peak, MAIN_ITERATOR check){
  
  // check the fragment mass difference:
  double targetMZ = check->first;
  
   
  // compare the precursor Difference:
  // check out of range:
  if( ! simple_math2::compareMassValuesAtPPMLevel( peak->get_MZ(), targetMZ, MZ_TOLERANCE * 4.0 ) ){
    return -1;
  }
  
  // compare the precursor Difference:
  if( ! simple_math2::compareMassValuesAtPPMLevel( peak->get_MZ(), targetMZ, MZ_TOLERANCE) ){
    return 0;
  }
  
  // get the last ms peak in the elution peak:
  
  MZ_series_REV_ITERATOR it = (check->second).rbegin();
  ms_peak* lastPeak = &( it->rbegin()->second);
  // charge state:
  if( peak->get_Chrg() != lastPeak->get_Chrg() ){
    //return 0;
  }
   
  return 1;
}


///////////////////////////////////////////////////////////////////////////////////////
// converts DeconvPeak list to ms_peak vector
void Process_Data::convert_ms_peaks(int SCAN, double TR, list<DeconvPeak>& DECONVPEAK,vector<ms_peak>& MSPEAK){
  
  list<DeconvPeak>::iterator mpi;
  for (mpi=DECONVPEAK.begin();mpi!=DECONVPEAK.end();++mpi) {
    ms_peak peak(SCAN,mpi->getMass(),mpi->getIntensity(),mpi->getCharge(),mpi->getNrIsotopes(),mpi->getScore(),mpi->getIsotopicPeaks());

    if( !mpi->getExtraPeakInfo().empty() ){
      peak.setExtraPeakInfo( mpi->getExtraPeakInfo() );
    }
    
    peak.set_retention_time( TR );
    //peak.show_info();
    MSPEAK.push_back(peak);
  }
}




///////////////////////////////////////////////////////////////////////////////
// go back to the MS1 level and
// find the correct precursor mass by mz and z:
void Process_Data::adjustCorrectToMS1Precursor( double* precursorMZ, int z, int MS1scan, int MS2Scan ){
  
  // if higher isotope picked for MS/MS, then need to start searching 
  // the monoisotopic mass at lower m/z value:
  ms_peak* preCursorPeak = NULL;
  double saveIsotopeDistance = 6;
  double searchMzLowerBound = *precursorMZ - saveIsotopeDistance;
  MAIN_ITERATOR P = pMZ_LIST.lower_bound( searchMzLowerBound );
  
  while( P != pMZ_LIST.end() ){
  
    MZ_series_REV_ITERATOR Pend = (P->second).rbegin();
    ms_peak* myPeak = &(Pend->rbegin()->second);
    // compare the charge states:
    if( myPeak->get_Chrg() == z ){
      // compare the scan numbers:
      int deltaScan = myPeak->get_Scan() - MS1scan;
      if( (int)fabs((double)deltaScan) <= getMaxScanDistance()){      
        if( myPeak->checkIsotopeBelongingAndAdjustMass( *precursorMZ, MZ_TOLERANCE) ){
          // store the precursor peak:
          preCursorPeak = myPeak;
          break;
        }
      }
    }
    
    
    ////////////////////////////
    // break out condition:
    // give tolerance due to wrong isotope picked for MS/MS
    double deltaM = myPeak->get_MZ() - *precursorMZ;
    if( deltaM > simple_math2::getMassErrorAtPPMLevel(myPeak->get_MZ(), 5 * MZ_TOLERANCE) ){
      break;
    }
    
    P++;
  }
  
  if( preCursorPeak != NULL ){
    // store the MS/MS scan number and activate this peak as precursor peak:
    preCursorPeak->activateAsPrecursorPeak( MS2Scan );
    // cout<<*precursorMZ<<"->"<<preCursorPeak->get_MZ()<<" :"<<MS2Scan<<endl;
    *precursorMZ = preCursorPeak->get_MZ();
  }
  else{
    //cout<<*precursorMZ<<" not found"<<endl;
  }
  
}

}
