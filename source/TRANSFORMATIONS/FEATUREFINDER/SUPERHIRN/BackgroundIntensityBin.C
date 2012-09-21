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

#include <string>
#include <vector>
#include <map>
#include <math.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/ms_peak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/BackgroundIntensityBin.h>

using namespace std;

double BackgroundIntensityBin::TR_BINS = 2.0;
double BackgroundIntensityBin::MZ_BINS = 50.0;
double BackgroundIntensityBin::INTENS_BINS = 50;
int BackgroundIntensityBin::MIN_BIN_COUNT = 1;

BackgroundIntensityBin::BackgroundIntensityBin(double mz, double tr){
  mzCoord = mz;
  trCoord = tr;
  zCoord = -1;
}

// check if a peak belongs to this intenity bin
bool BackgroundIntensityBin::checkBelonging(ms_peak* peak){
  
  // check charge state:
  if( zCoord != -1 ){
    if( peak->get_charge_state() != zCoord ){
      return false;
    }
  }
  
  // check tr:
  double tr = peak->get_retention_time();
  if( ( tr < (trCoord - TR_BINS /2.0) ) || ( tr > (trCoord + TR_BINS /2.0) ) ){
    return false;
  }
  
  double mz = peak->get_MZ();

  if( ( mz < (mzCoord - MZ_BINS /2.0) ) || ( mz > (mzCoord + MZ_BINS /2.0) ) ){
    return false;
  }
  
  addIntensity( peak->get_intensity() );
  peak = NULL;
  return true;
}

void BackgroundIntensityBin::addMSPeak( ms_peak* peak ){
  addIntensity( peak->get_intensity() );
  peak = NULL;
}

void BackgroundIntensityBin::addIntensity( double intens ){
  IntensityMap.push_back( intens );
}

// copied from simple_math
double simple_math_WEIGHTED_AVERAGE(map< double, double >* IN){
  
  double AVERAGE = 0;
  double TOT_WEIGHT = 0;  
  
  if( IN->size() > 1 ){
    
    map< double, double >::iterator START = IN->begin();  
    while( START != IN->end() ){
      TOT_WEIGHT += (*START).second;
      AVERAGE += ( (*START).first * (*START).second );
      START++;
    }
    
    return AVERAGE / TOT_WEIGHT;  
  }
  else{
    return (*IN->begin()).first;
    IN = NULL;
  }
}

// process collected intensities in the map
void BackgroundIntensityBin::processIntensities( ){
  
  computeIntensityHist( );
  
  if(!IntensityHist.empty()){
    mean = simple_math_WEIGHTED_AVERAGE( &IntensityHist );
  }
  else{
    mean = 0;
  }
}

// copmute an intensity histogram
void BackgroundIntensityBin::computeIntensityHist( ){
  
  double constraint = BackgroundIntensityBin::INTENS_BINS;
  
  // insert into the histogram map
  vector<double>::iterator P = IntensityMap.begin();
  while( P != IntensityMap.end() ){
    
    // intensity to bin:
    double intens = (*P);
    
    // find a key:
    map<double, double>::iterator F = IntensityHist.lower_bound( intens ) ;
    if( F != IntensityHist.end() ){
      
      // check this one:
      map<double, double>::iterator check = F;      
      double mainLow = fabs( check->first - intens );
      double deltaHigh = 1000000;
      if( check != IntensityHist.begin() ){
        check--;
        deltaHigh = fabs( check->first - intens );
        if( mainLow > deltaHigh ){ 
          mainLow = deltaHigh;
          F = check;
        }
      }
      if( mainLow > constraint){
        F = IntensityHist.end();
      }
      else{
        F->second += 1.0;
      }
    }  
    
    if( F == IntensityHist.end() ){
      IntensityHist.insert( make_pair( intens, 1.0 ) );
    }
          
    P++;
  }
 
  // filter out bins of only 1 counts:
  map<double, double>::iterator F = IntensityHist.begin( ) ;
  while( F != IntensityHist.end() ){
    
    if( F->second == MIN_BIN_COUNT ){
      IntensityHist.erase( F++ );
    }
    else{
      F++;
    }
  }
}

BackgroundIntensityBin::~BackgroundIntensityBin(){
  IntensityMap.clear();
}

