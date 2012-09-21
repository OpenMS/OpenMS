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
//
// **********************************************************************//
// CLASS ms_peak:
//
// variable description:
//
//
// function description:
//
//
// **********************************************************************//


#ifndef MS_PEAK_H
#define MS_PEAK_H

#include <string>
#include <vector>
#include <map>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/ms_peak.h>

class CentroidPeak;

class ms_peak{

    
  ////////////////////////////////////////////////
  // declaration of the private members:
  
  double precursorMZ;
  double MZ;
  float INTENSITY;
  int SCAN;
  double TR;
  unsigned int CHRG;
  unsigned int NRISOTOPES;
  float SCORE;
  
  std::string extraMSPeakInfo;
  
  // child scan options:
  bool precursorMass;
  int childScan;
  
  double SignalToNoise;
  std::vector<CentroidPeak> ISOPEAKS;
  
private:

    ////////////////////////////////////////////////
    // declaration of the public members:
    
public:
    
    // class destructor
    ~ms_peak();

    // class constructor
    ms_peak(int,double,float);
    ms_peak();
    ms_peak(int,double,float, unsigned int,unsigned int, float,std::vector<CentroidPeak>);
    ms_peak(const ms_peak&);
    ms_peak(const ms_peak*);

    //////////////////////////////////////////////////
    // overload operators:
    ms_peak& operator=(const ms_peak&);
    ms_peak& operator<=(const ms_peak&);
    ms_peak& operator>=(const ms_peak&);
    ms_peak& operator<(const ms_peak&);
    ms_peak& operator>(const ms_peak&);
    
    // show content of peak:
  void show_info();
  
  // store the MS/MS scan number and activate this peak as precursor peak:
  void activateAsPrecursorPeak( int );
    
    
    
  
  //////////////////////////////////////////////////
  // check if the input mass matches one of the isotopic masses
  bool checkIsotopeBelongingAndAdjustMass(double, double );
    
  
  ///////////////////////////////
  // start here all the get / set
  // function to access the
  // variables of the class
  
  std::vector<CentroidPeak>& get_isotopic_peaks() {return ISOPEAKS;}
  std::vector<CentroidPeak>::iterator get_isotopic_peaks_start() {return ISOPEAKS.begin();}
  std::vector<CentroidPeak>::iterator get_isotopic_peaks_end() {return ISOPEAKS.end();}
  
  void setExtraPeakInfo(std::string in){extraMSPeakInfo = in;};
  std::string getExtraPeakInfo(){return extraMSPeakInfo;};
  
  // precursor mass of the ms2 scane:
  void setPrecursorMZ(double in){ precursorMZ = in;};
  double getPrecursorMZ(){ return precursorMZ;};
  
  // precursor mass charge state:
  void setPrecursorCHRG(int in){ CHRG = in;};
  int getPrecursorCHRG(){ return CHRG;};
  
  // check if this peak has been determined as precursor:
  bool getPrecursorActivation( ){return precursorMass;};
                               
  int get_Chrg(){return CHRG;};
  void set_Chrg(int z){CHRG = z;};
  int get_Scan(){return SCAN;};
  float get_intensity(){return INTENSITY;};
  double get_MZ(){return MZ;};
  int get_scan_number(){return SCAN;};
  void set_retention_time(double IN){TR = IN;};
  double get_retention_time(){return TR;};
  unsigned int get_charge_state(){return CHRG;};
  unsigned int get_nr_isotopes(){return NRISOTOPES;};
  float get_score(){return SCORE;};
  
  
  double getSignalToNoise(){return SignalToNoise;};
  void setSignalToNoise(double in){SignalToNoise = in;};
};
#endif
