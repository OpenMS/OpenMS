///////////////////////////////////////////////////////////////////////////
//
//  PEAK DETECTION OF FOURIER TRANSFORME MS INSTRUMENT DATA
//
//  written by Markus Mueller, markus.mueller@imsb.biol.ethz.ch
//  ( and Lukas Mueller, Lukas.Mueller@imsb.biol.ethz.ch)
//  October 2005
//  
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
//  
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
// 
//

#include <list>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/CentroidPeak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/ms_peak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/BackgroundIntensityBin.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/LC_MS_XML_reader.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/BackgroundControl.h>

using namespace std;

BackgroundControl::BackgroundControl(){
  init();
}

BackgroundControl::~BackgroundControl(){
  intensityBinMap.clear();
}

void BackgroundControl::init(){
  
  // create a vector of intensity bin objects
  
  // first in the tr dimension:
  double trStart = LC_MS_XML_reader::TR_MIN;
  while( trStart <= LC_MS_XML_reader::TR_MAX ){
  
    // inner loop is the mzBins:
    map<double, BackgroundIntensityBin> mzArray;
    double mzStart = LC_MS_XML_reader::FEATURE_MZ_MIN;
    while( mzStart <= LC_MS_XML_reader::FEATURE_MZ_MAX ){
      
      BackgroundIntensityBin* bin = new BackgroundIntensityBin(mzStart,trStart);
      mzArray.insert( make_pair( mzStart, *bin ) );
      delete bin;
      bin = NULL;
    
      mzStart += BackgroundIntensityBin::MZ_BINS;
    }
    
    intensityBinMap.insert( make_pair (trStart, mzArray) );
    trStart += BackgroundIntensityBin::TR_BINS;
  }
  
}

void BackgroundControl::addPeakMSScan( double TR, list<CentroidPeak>* peakList ){
  
  map<double, map< double, BackgroundIntensityBin>  >::iterator F = findTrKey( TR );
  if( F != intensityBinMap.end() ){
    
    // find the mz bins:
    map< double, BackgroundIntensityBin>* mzMap = &(F->second);
    
    list<CentroidPeak>::iterator mpi;
    for (mpi=peakList->begin();mpi!=peakList->end();++mpi) {
      
      map< double, BackgroundIntensityBin>::iterator F_mz = findMzKey( mpi->getMass(), mzMap );
      if(  F_mz != mzMap->end() ){
        F_mz->second.addIntensity( mpi->getIntensity() );
      }
    }
  }
}

map< double, BackgroundIntensityBin>::iterator BackgroundControl::findMzKey( double mz, map< double, BackgroundIntensityBin>* mzMap ){
  
  double constraint = BackgroundIntensityBin::MZ_BINS / 2.0;
  map<double, map< double, BackgroundIntensityBin>::iterator > outMap;
  
  map< double, BackgroundIntensityBin>::iterator F = mzMap->lower_bound( mz );
  if( F != mzMap->end() ){
    double delta = fabs( F->first - mz );
    if( delta <= constraint ){ 
      outMap.insert( make_pair( delta, F) );
    }
  }
  
  
  if( F != mzMap->begin() ){
    F--;
    double delta = fabs( mz - F->first );    
    if( delta <= constraint ){
      outMap.insert( make_pair( delta, F) );
    }
  }
  
  if( !outMap.empty() ){
    return outMap.begin()->second;
  }
  
  return mzMap->end();
}

map<double, map< double, BackgroundIntensityBin>  >::iterator BackgroundControl::findTrKey( double Tr ){
  
  double constraint = BackgroundIntensityBin::TR_BINS * 2;
  
  map<double, map<double, map< double, BackgroundIntensityBin>  >::iterator > outMap;
  map<double, map< double, BackgroundIntensityBin>  >::iterator F = intensityBinMap.lower_bound( Tr );
  if( F != intensityBinMap.end() ){
    double delta = fabs( Tr - F->first );    
    if( delta <= constraint ){
      outMap.insert( make_pair( delta, F) );
    }
    
  }

  if( F != intensityBinMap.begin() ){
    F--;

    double delta = fabs( Tr - F->first );    
    if( delta <= constraint ){
      outMap.insert( make_pair( delta, F) );
    }
  }
  
  if( !outMap.empty() ){
    return outMap.begin()->second;
  }

  return intensityBinMap.end();
}

double BackgroundControl::getBackgroundLevel( ms_peak* in){
  return getBackgroundLevel( in->get_MZ(), in->get_retention_time()  );
}

double BackgroundControl::getBackgroundLevel( double mz, double tr){
  // find the corresponding retention time bin:
  map<double, map< double, BackgroundIntensityBin> >::iterator F = findTrKey( tr  );
  if( F != intensityBinMap.end() ){
    map< double, BackgroundIntensityBin>::iterator F2 = findMzKey( mz, &(F->second) );
    if( F2 != F->second.end() ){
      return F2->second.getMean();
    }
  }
  return -1.0;
}

void BackgroundControl::processIntensityMaps(  ){
  
  map<double, map< double, BackgroundIntensityBin> >::iterator P1 = intensityBinMap.begin();
  while( P1 != intensityBinMap.end() ){
    
    map< double, BackgroundIntensityBin> ::iterator P2 = P1->second.begin();
    while( P2 != P1->second.end() ){
      
      P2->second.processIntensities( );
      P2++;
    }
    
    P1++;
  }
}

