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

#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <math.h>
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
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/Process_Data.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/ms2_info.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2_feature.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/featureLCprofile.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/feature.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/LC_MS.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS1_feature_merger.h>

namespace OpenMS
{

using namespace std;

double MS1_feature_merger::INTENSITY_APEX_THRESHOLD;
double MS1_feature_merger::MS1_PEAK_AREA_TR_RESOLUTION;

double MS1_feature_merger::INITIAL_TR_TOLERANCE;
double MS1_feature_merger::MS1_FEATURE_MERGING_TR_TOLERANCE;
double MS1_feature_merger::PERCENTAGE_INTENSITY_ELUTION_BORDER_VARIATION;
double MS1_feature_merger::PPM_TOLERANCE_FOR_MZ_CLUSTERING;
bool MS1_feature_merger::MS1_FEATURE_CLUSTERING;
////////////////////////////////////////////////
// constructor for the object ana_summarizer:
MS1_feature_merger::MS1_feature_merger(LC_MS* in){
  lcmsMap = in;
}

//////////////////////////////////////////////////
// class desctructor
MS1_feature_merger::~MS1_feature_merger(){
  lcmsMap = NULL;
}



//////////////////////////////////////////////////
// start the merging process
void MS1_feature_merger::startFeatureMerging(){

  printf( "\t\t -- merging features in LC-MS %s: ", lcmsMap->get_spec_name().c_str());   

  int before = lcmsMap->get_nb_features();
  
  // iterative approach to merge features:
  unsigned int startNbFeatures = -1;
  while( startNbFeatures != lcmsMap->get_nb_features() ){
  
    startNbFeatures = lcmsMap->get_nb_features();
    
    // create mz feature clusters
    createMZFeatureClusters();
    
    // process the mz clusters
    map< double, vector<feature*> >::iterator M = mzClusters.begin( );
    while( M != mzClusters.end() ){
      if( M->second.size() > 1 ){
        processMZFeatureVector( &( M->second ) );
      }
      M++;
    }
    
    // remove the merged features:
    vector<int>::iterator I = idsToRemove.begin();
    while( I != idsToRemove.end() ){
      lcmsMap->remove_feature_by_ID( *I );
      I++;
    }
    
    mzClusters.clear();
    idsToRemove.clear();
  
  }    
    
  int after = lcmsMap->get_nb_features();
  printf( "%d merged\n", before - after);   
  

}


//////////////////////////////////////////////////
// create a distribution of delta Tr for the splited features
void MS1_feature_merger::createMZFeatureClusters(){
  
  vector<LC_MS_FEATURE>::iterator I = lcmsMap->get_feature_list_begin();
  while( I != lcmsMap->get_feature_list_end() ){
    feature* fea = &(*I);
    
    map< double, vector<feature*> >::iterator F = mzClusters.lower_bound( fea->get_MZ() );
    
    if( mzClusters.empty() ){
      vector<feature*> tmp;
      tmp.push_back( fea );
      mzClusters.insert( make_pair( fea->get_MZ(), tmp) );
      
    }
    else if( F == mzClusters.begin() ){
      
      // check below:
      if( compareMZFeatureBeloning( fea, (*F->second.begin()) ) ){
        F->second.push_back( fea );
      }
      else{
        vector<feature*> tmp;
        tmp.push_back( fea );
        mzClusters.insert( make_pair( fea->get_MZ(), tmp) );        
      }
      
    }
    else if( F == mzClusters.end() ){
      F--;      
      // check below:
      if( compareMZFeatureBeloning( fea, (*F->second.begin())) ){
        F->second.push_back( fea );
      }
      else{
        vector<feature*> tmp;
        tmp.push_back( fea );
        mzClusters.insert( make_pair( fea->get_MZ(), tmp) );
      }
    }
    else{
      
      bool found = false;
      if( compareMZFeatureBeloning( fea, (*F->second.begin())) ){
        found = true;
        F->second.push_back( fea );
      }
      F--;
      // check below:
      if( compareMZFeatureBeloning( fea, (*F->second.begin()) )){
        found = true;
        F->second.push_back( fea );
      }      
      
      if( !found ){
        vector<feature*> tmp;
        tmp.push_back( fea );
        mzClusters.insert( make_pair( fea->get_MZ(), tmp) );
      }
      
    }
    
    
    I++; 
}
  
 
}


//////////////////////////////////////////////////
// compare if a feature belongs to another feature
bool MS1_feature_merger::compareMZFeatureBeloning(feature* A, feature* B){
  
  if( ( A->getLCelutionProfile( ) == NULL ) || ( B->getLCelutionProfile( ) == NULL ) ){
    return false;
  }
    
  if( ( A->getLCelutionProfile( )->getNbLCelutionSignals() == 0 ) || ( B->getLCelutionProfile( )->getNbLCelutionSignals() == 0 ) ){
    return false;
  }

  // check mz level:
  if( !simple_math2::compareMassValuesAtPPMLevel( A->get_MZ( ), B->get_MZ( ), MS1_feature_merger::PPM_TOLERANCE_FOR_MZ_CLUSTERING ) ){
    return false;
  }

  // check charge state level:
  if( A->get_charge_state() != B->get_charge_state() ){
    return false;
  }      
  
  return true;
}  

//////////////////////////////////////////////////
// process a vector of m/z features
void MS1_feature_merger::processMZFeatureVector(vector<feature*>* mapF){

  unsigned int length = -1;
  
  // order first the features according to retention time:
  sort(mapF->begin(),mapF->end(), OPERATOR_FEATURE_TR());
  
  // go through the vector to find to be merged features:
  while( length != mapF->size() ){
    
    length = (unsigned int) mapF->size();
    vector<feature*>::iterator SEARCHER = mapF->begin();
    while( SEARCHER != mapF->end() ){
            
      // take the first one as a reference:
      feature* loc = *SEARCHER;
      
      // find features to merge:
      vector<feature*>::iterator F = SEARCHER;
      F++;
      F = findFeaturesToMerge( loc, F, mapF );
      SEARCHER++;
    }
  }
}


//////////////////////////////////////////////////
// find to this feature the features which should be merged
vector<feature*>::iterator MS1_feature_merger::findFeaturesToMerge(feature* search, vector<feature*>::iterator mapCurrent, vector<feature*>* Fmap){
  
  bool log10Intens = true;
  
  // get the elution profile:
  featureLCprofile* searchLC = search->getLCelutionProfile( );
  while( mapCurrent != Fmap->end() ){
  
    // check absolute retention time difference:
    bool toMerge = false;
    feature* mergedTarget = (*mapCurrent);
    double deltaTr = fabs( search->get_retention_time() - mergedTarget->get_retention_time() ); 
    if( deltaTr <= MS1_feature_merger::INITIAL_TR_TOLERANCE){

      // compare the end / start of the elution peak:
      MS1Signal* start;
      MS1Signal* end;
      if( search->get_retention_time() < mergedTarget->get_retention_time() ){
        start = &(searchLC->getLastLCelutionSignal()->second);
        end = &(mergedTarget->getLCelutionProfile( )->getLCelutionSignalsStart()->second);
      }
      else{
        end = &(mergedTarget->getLCelutionProfile( )->getLastLCelutionSignal()->second);
        start = &(searchLC->getLCelutionSignalsStart()->second);          
      }
      double startIntens = start->intensity;
      if(log10Intens){
        startIntens = log10(startIntens);
      }
      double endIntens = end->intensity;
      if (log10Intens){
        endIntens = log10(endIntens);
      }
      
      double deltaIntens = fabs( startIntens - endIntens);
      deltaIntens /= startIntens; 
      deltaTr = fabs( start->TR - end->TR );         
      if( ( deltaTr <= MS1_FEATURE_MERGING_TR_TOLERANCE) && ( deltaIntens <= PERCENTAGE_INTENSITY_ELUTION_BORDER_VARIATION) ){
        toMerge = true;
      }
    
    }
    
    if( toMerge ){
      // search->show_info();(*mapCurrent)->show_info();
      mergeFeatures( search, mergedTarget );
      idsToRemove.push_back( mergedTarget->get_feature_ID() );
      mapCurrent = Fmap->erase( mapCurrent );
      
      if( search->get_peak_area() == 0 ){
        idsToRemove.push_back( search->get_feature_ID() );
        mapCurrent++;
        break;
      }
    }
    else{
      mapCurrent++; 
    }
  }
  
  return mapCurrent;
}
      


/////////////////////////////////////////////////////////////////////////////
// merge the target to the search feature
void MS1_feature_merger::mergeFeatures( feature* target, feature* toMerge ){
  
  
  double TOT_AREA = target->get_peak_area() + toMerge->get_peak_area(); 
  // merge the m/z:
  target->set_MZ( (target->get_peak_area()*target->get_MZ( ) + toMerge->get_peak_area()*toMerge->get_MZ() ) / TOT_AREA );
  // merge S/N:
  target->setSignalToNoise( (target->getSignalToNoise()*target->get_peak_area( ) + toMerge->getSignalToNoise()*toMerge->get_peak_area() ) / TOT_AREA );
  // merge score:
  target->set_peak_score( (target->get_peak_score()*target->get_peak_area( ) + toMerge->get_peak_score()*toMerge->get_peak_area() ) / TOT_AREA );
  
  
  // merge first the elution profiles:  
  featureLCprofile* targetLC = target->getLCelutionProfile( );
  featureLCprofile* mergeLC = toMerge->getLCelutionProfile( );

  // add points of the toMerge to the target:
  map<int, MS1Signal>::iterator LC = mergeLC->getLCelutionSignalsStart();
  while( LC != mergeLC->getLCelutionSignalsEnd() ){
    targetLC->addMS1elutionSignal( &(LC->second));
    LC++; 
  }

  // possible extra info:
  if( target->getFeatureExtraInformation( ).empty() ){
    target->setFeatureExtraInformation( toMerge->getFeatureExtraInformation( ) );
  }
  
  // compute new parameters
  computeNewMS1FeatureParameters( target );
  
  // copy MS/MS information:
  if( toMerge->get_MS2_info( -3.0 ) ){
    target->add_MS2_info( toMerge->get_MS2_SCAN_MAP( ) );
  }
    
}



//////////////////////////////////////////////////////////////////
// Compute a varietiy of parameters for the LC elution peak
void MS1_feature_merger::computeNewMS1FeatureParameters( feature* in ){
  
  
  featureLCprofile* lcProfile = in->getLCelutionProfile( );
  
  // define the apex treshold:
  double maxIntens = -1;
  map<int, MS1Signal>::iterator LC = lcProfile->getLCelutionSignalsStart();
  while( LC != lcProfile->getLCelutionSignalsEnd() ){
    if( maxIntens < (*LC).second.intensity ){
      maxIntens = (*LC).second.intensity;
    }
    LC++; 
  }
  
  // get the MS peak above noise to copmute:
  double THRESHOLD = maxIntens / in->getSignalToNoise();
  vector< MS1Signal* > computeMap;
  LC = lcProfile->getLCelutionSignalsStart();
  
  in->set_scan_start( (*LC).second.scan );
  in->set_retention_time_START( (*LC).second.TR );
  
  while(LC != lcProfile->getLCelutionSignalsEnd()){
    if( (*LC).second.intensity >= THRESHOLD){
      computeMap.push_back( &(LC->second) );
    }
    LC++;
  }
  LC--;
  in->set_scan_end( (*LC).second.scan );
  in->set_retention_time_END( (*LC).second.TR );
  

  
  if( !computeMap.empty() ){
  
    vector< MS1Signal* >::iterator P = computeMap.begin();
    
    double TOT_AREA = 0;
    double start_TR = 0;
    double start_int = 0;
    double apexScan = 0;
    double apexTr = 0;
    double end_TR = 0;  
    double end_int = 0;
    
    start_TR = (*P)->TR;
    start_int = (*P)->intensity;
    P++;
    // go through all peaks in the LC elution profile:
    while(P != computeMap.end()){
      
      if( (*P)->intensity >= THRESHOLD ){
        
        end_TR = (*P)->TR; 
        end_int = (*P)->intensity;
        
        // compute an area between local start / end ms peak:
        double area = computeDeltaArea(start_TR, start_int - THRESHOLD, end_TR, end_int - THRESHOLD);
        TOT_AREA += area;
        apexScan += (double)((*P)->scan) * area;
        apexTr += start_TR * area;
        
        // next scan:
        start_TR = end_TR;
        start_int = end_int;
      }
      
      P++;
    }
    
    // if contained only one peak!
    if(computeMap.size() == 1){
      in->set_peak_area( (float) start_int );
      in->set_retention_time ( in->get_retention_time_START() );
      in->set_scan_number( in->get_scan_start() );
    }
    else{    
      in->set_peak_area( (float) TOT_AREA);
      apexScan /= TOT_AREA;
      in->set_scan_number( (int) apexScan );
      apexTr /= TOT_AREA;
      in->set_retention_time ( apexTr );
    }
    
    // set the apex ms peak:
    LC = lcProfile->getLCelutionSignalMap()->lower_bound( in->get_scan_number() );
    in->set_apex_peak_intensity( (*LC).second.intensity);
  }
  else{
    // no good peak above threshold, so reset all the features parameters to remove the feature
    in->set_peak_area( 0 );
    in->set_scan_number( 0 );
    in->set_retention_time ( 0 );
    
  }
  
}

/////////////////////////////////////////////////////////////////////
// computes the area of between 2 peaks:
double MS1_feature_merger::computeDeltaArea(double START_TR, double START_INT, double END_TR, double END_INT){

  double AREA = 0;
  if( ( START_INT > 0 ) && ( END_INT > 0 ) && (START_TR <= END_TR) ){
    double x = (END_TR - START_TR) / MS1_feature_merger::MS1_PEAK_AREA_TR_RESOLUTION;
    double y = (END_INT - START_INT);
    
    if( ( y != 0 ) && ( x != 0 ) ){
      double m = y/x;
      double INT = START_INT;
      double count = 0;
      while( count <= x){
        AREA += INT;  
        INT += m;
        count++;
      }
      AREA += INT;
    }    
  }
  return AREA;
}

  






}
