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


#ifndef USE_BACKGROUND_INTENSITY_BIN_H
#define USE_BACKGROUND_INTENSITY_BIN_H

#include <string>
#include <vector>
#include <map>

using namespace std;


class BackgroundIntensityBin {

private:

  // mz and tr coordinates of the bin:
  double mzCoord;
  double trCoord;
  double zCoord;  
  
  vector<double> IntensityMap;  
  map<double, double> IntensityHist;

  double mean;
  double median;
  void computeIntensityHist();

public:
  
  static double TR_BINS;
  static double MZ_BINS;
  static double INTENS_BINS;
  static int MIN_BIN_COUNT;
  
  ~BackgroundIntensityBin();

  BackgroundIntensityBin(double, double);
  
  // check if a peak belongs to this intenisty bin
  bool checkBelonging(ms_peak*);
  // add intensity to BackgroundIntensityBin
  void addIntensity( double );
  // add peak to BackgroundIntensityBin
  void addMSPeak( ms_peak* );
  // process collected intensities in the map
  void processIntensities( );

  vector<double>* getIntensityMap(){ return &IntensityMap;};
  map<double, double>* getIntensityHist(){ return &IntensityHist;};
  double getMean(){return mean;};
};

#endif

    
