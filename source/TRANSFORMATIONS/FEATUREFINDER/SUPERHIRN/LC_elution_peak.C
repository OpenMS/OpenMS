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

#include <math.h>
#include <stdio.h>
#include <iostream>
#include <map>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/consensIsotopePattern.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/ms_peak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/CentroidPeak.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/LC_elution_peak.h>

namespace OpenMS
{

using namespace std;


// resolution of the retention time, for peak area copmuting:
float LC_elution_peak::TR_RESOLUTION;

// cut off, where everything small than this precentile of the
// apex is discarded
// float LC_elution_peak::intensity_apex_percentil_cutoff;

// parameters to debug a ceratain mass range
double LC_elution_peak::DEBUG_MASS_START = -1;
double LC_elution_peak::DEBUG_MASS_END = -1;


////////////////////////////////////////////////
// constructor for the object LC_elution_peak:
LC_elution_peak::LC_elution_peak(){

  f_observed_Mass = 0;
  fIsotopMass = 0;
  fMonoMass = 0;
  fVolume = 0;
  fCharge = 0;
  fNrIsotopes = 0;
  fScanNumberStart = 0;
  fScanNumberApex = 0;
  fScanNumberEnd = 0;
  fpeak_area = 0;
  fapex_intensity = 0;
  fRT = 0;
  fStartTR = 0;
  fEndTR = 0;
  isotopePattern = NULL;
  
}

////////////////////////////////////////////////
// constructor for the object LC_elution_peak:
LC_elution_peak::LC_elution_peak(MZ_series_ITERATOR data, double MZ){
  
  f_observed_Mass = MZ;
  fIsotopMass = 0;
  intens_signals = (*data);
  fMonoMass = 0;
  fVolume = 0;
  fCharge = 0;
  fNrIsotopes = 0;
  fScanNumberStart = 0;
  fScanNumberApex = 0;
  fScanNumberEnd = 0;
  fapex_intensity = 0;
  fpeak_area = 0;
  fRT = 0;
  fStartTR = 0;
  fEndTR = 0;
  isotopePattern = NULL;
}


//////////////////////////////////////////////////
// class desctructor of LC_elution_peak
LC_elution_peak::~LC_elution_peak(){
  
  intens_signals.clear();
  CHRG_MAP.clear();  
  if( isotopePattern != NULL ){
    delete isotopePattern;
    isotopePattern = NULL;
  }
}

////////////////////////////////////////////////
// constructor for the object feature:
LC_elution_peak::LC_elution_peak(const LC_elution_peak* tmp){
  CHRG_MAP = tmp->CHRG_MAP;
  f_observed_Mass = tmp->f_observed_Mass;
  fpeak_area = tmp->fpeak_area;
  fapex_intensity = tmp->fapex_intensity;
  fIsotopMass = tmp->fIsotopMass;
  fMonoMass = tmp->fMonoMass;
  fVolume = tmp->fVolume;
  fCharge = tmp->fCharge;
  fNrIsotopes = tmp->fNrIsotopes;
  fScanNumberStart = tmp->fScanNumberStart;
  fScanNumberApex = tmp->fScanNumberApex;
  fScanNumberEnd = tmp->fScanNumberEnd;
  fRT = tmp->fRT;
  fStartTR = tmp->fStartTR;
  fEndTR = tmp->fEndTR;
  intens_signals = tmp->intens_signals;
  fSignalToNoise = tmp->fSignalToNoise;
  fSNIntensityThreshold = tmp->fSNIntensityThreshold;
  isotopePattern = new consensIsotopePattern( tmp->isotopePattern );
  elutionPeakExtraInfo = tmp->elutionPeakExtraInfo;
}


//////////////////////////////////////////////////
// class copy constructor of LC_elution_peak
LC_elution_peak::LC_elution_peak(const LC_elution_peak& tmp){
  
  CHRG_MAP = tmp.CHRG_MAP;
  f_observed_Mass = tmp.f_observed_Mass;
  fpeak_area = tmp.fpeak_area;
  fapex_intensity = tmp.fapex_intensity;
  fIsotopMass = tmp.fIsotopMass;
  fMonoMass = tmp.fMonoMass;
  fVolume = tmp.fVolume;
  fCharge = tmp.fCharge;
  fNrIsotopes = tmp.fNrIsotopes;
  fScanNumberStart = tmp.fScanNumberStart;
  fScanNumberApex = tmp.fScanNumberApex;
  fScanNumberEnd = tmp.fScanNumberEnd;
  fRT = tmp.fRT;
  fStartTR = tmp.fStartTR;
  fEndTR = tmp.fEndTR;
  fSNIntensityThreshold = tmp.fSNIntensityThreshold;
  intens_signals = tmp.intens_signals;
  fSignalToNoise = tmp.fSignalToNoise;
  isotopePattern = new consensIsotopePattern( tmp.isotopePattern );
  elutionPeakExtraInfo = tmp.elutionPeakExtraInfo;

}


//////////////////////////////////////////////////
// copy constructor:
LC_elution_peak& LC_elution_peak::operator=(const LC_elution_peak& tmp){
  
  CHRG_MAP = tmp.CHRG_MAP;
  f_observed_Mass = tmp.f_observed_Mass;
  fpeak_area = tmp.fpeak_area;
  fapex_intensity = tmp.fapex_intensity;
  fIsotopMass = tmp.fIsotopMass;
  fMonoMass = tmp.fMonoMass;
  fVolume = tmp.fVolume;
  fCharge = tmp.fCharge;
  fNrIsotopes = tmp.fNrIsotopes;
  fScanNumberStart = tmp.fScanNumberStart;
  fScanNumberApex = tmp.fScanNumberApex;
  fScanNumberEnd = tmp.fScanNumberEnd;
  fRT = tmp.fRT;
  fStartTR = tmp.fStartTR;
  fEndTR = tmp.fEndTR;
  intens_signals = tmp.intens_signals;
  fSNIntensityThreshold = tmp.fSNIntensityThreshold;
  fSignalToNoise = tmp.fSignalToNoise;
  isotopePattern = new consensIsotopePattern( tmp.isotopePattern );
  elutionPeakExtraInfo = tmp.elutionPeakExtraInfo;

  return *this;
}


//////////////////////////////////////////////////////////////////
// get the original M/Z of a ms_peak
double LC_elution_peak::get_MZ(int IN){
  
  SIGNAL_iterator P = intens_signals.lower_bound(IN);
  if((*P).first == IN){
    return (*P).second.get_MZ();
  }
  
  if(P == get_signal_list_end()){
    P--;
    return (*P).second.get_MZ();
  }

  if(P == get_signal_list_start()){
    return (*P).second.get_MZ();
  }
  
  double SCAN_UP = (*P).first;
  P--;
  double SCAN_DOWN = (*P).first;
  
  if( (SCAN_UP - IN) <= (IN - SCAN_DOWN) )
    P++;
  
  return (*P).second.get_MZ();
  
}


////////////////////////////////////////////////////////////////////
// find the closest existing mz peak in the elution profile:
ms_peak* LC_elution_peak::find_true_peak(float SCAN){
  
  int int_SCAN = int(floor(SCAN));
  
  SIGNAL_iterator P = intens_signals.upper_bound(int_SCAN);
  if( P == intens_signals.end() ){
    P--;
    return &((*P).second);
  }
  else if( P == intens_signals.begin() ){
    return &((*P).second);
  }
  else{
    float dis_UP = float((*P).first) - SCAN;
    P--;
    float dis_down = SCAN - float((*P).first);
    if( dis_UP > dis_down){
      return &((*P).second);
    }
    else{
      P++;
      return &((*P).second);      
    }
  } 
}

/////////////////////////////////////////////
// print the elution profile from a peak:
//void LC_elution_peak::print_profile(ofstream* OUT){
//
//  char buffer[256];
//  
//  //  sprintf(buffer,"\n#peak[%f]\n", get_apex_MZ());  
//  // OUT->write(buffer, strlen(buffer));
//  SIGNAL_iterator P = get_signal_list_start();
//  while(P != get_signal_list_end()){
//    sprintf(buffer, "%f %d %d %f\n",get_apex_MZ(), (*P).first, (*P).second.get_charge_state(),(*P).second.get_intensity());
//    OUT->write(buffer, strlen(buffer));
//    P++;
//  }
//}



/////////////////////////////////////////////
// print the elution profile from a peak:
void LC_elution_peak::show_info(){

  /*
  cout<<endl;
  
  SIGNAL_iterator P = get_signal_list_start();
  while(P != get_signal_list_end()){
  cout<<get_apex_MZ()<<" "<<(*P).first<<" "<<(*P).second.get_intensity()<<endl;;
  P++;
  }
  */
  printf("scan:[%d,%d,%d], TR:[%0.2f,%0.2f,%0.2f],m/z=%0.4f(+%d),area=%0.2e(%0.2f),S/N=%0.2f\n",fScanNumberStart,fScanNumberApex,fScanNumberEnd,fStartTR,fRT,fEndTR,get_apex_MZ(), get_charge_state(),get_total_peak_area(),get_apex_intensity(), getSignalToNoise());
}


//////////////////////////////////////////////////////////////////
// determine the intensity background baseline based on S/N 
// value:
void LC_elution_peak::setSNIntensityThreshold(){
  
  fSignalToNoise = 0;
  fSNIntensityThreshold = 0;
  
  double TotArea = 0;
  SIGNAL_iterator P = get_signal_list_start();
  while(P != get_signal_list_end()){
    ms_peak* peak = &(P->second);
    fSignalToNoise += peak->getSignalToNoise() *  peak->get_intensity(); 
    fSNIntensityThreshold +=  peak->get_intensity() * ( peak->get_intensity()/ peak->getSignalToNoise()) ;
    TotArea += peak->get_intensity(); 
    P++;
  }
    
  // set the signal to noise:
  fSignalToNoise /= TotArea;
  // set the noise threshold:
  fSNIntensityThreshold /= TotArea;

}


/*
//////////////////////////////////////////////////////////////////
// define the apex into the elution profile::
void LC_elution_peak::define_apex(){
  
  // for apex calculation:
  double TOT_TR = 0;
  double TOT_INT = 0;
  double TOT_SCAN = 0;
  
  double start_TR = 0;
  double start_int = 0;
  double mid_TR = 0;  
  double mid_int = 0;	
  double mid_scan = 0;
  double end_TR = 0;  
  double end_int = 0;
      
  SIGNAL_iterator P = get_signal_list_start();
  
  // set start et first scan in the LC_peak:
  start_TR = (*P).second.get_retention_time();
  start_int = (*P).second.get_intensity();
  mid_TR = (*P).second.get_retention_time();
  mid_int = (*P).second.get_intensity();
  mid_scan = (*P).second.get_scan_number();


  P++;
  
  /////////////////////////////////////////////////////
  // go through all peaks in the LC elution profile:
  while(P != get_signal_list_end()){
    
    end_TR = (*P).second.get_retention_time(); 
    end_int = (*P).second.get_intensity();
      
    ///////////////////////////////////////////////////
    // copmute the area around the mid ms peak:
    // compute area between end/mid:
    double area_A = compute_delta_area( mid_TR, mid_int, end_TR, end_int);
    
    // compute area between mid/start:
    double area_B = compute_delta_area(start_TR, start_int, mid_TR, mid_int);
    
    // add up for apex calculation:
    TOT_TR += mid_TR * ( area_A + area_B );
    TOT_INT += ( area_A + area_B );  
    TOT_SCAN += mid_scan * ( area_A + area_B );
    
    // next scan:
    start_TR = mid_TR;
    start_int = mid_int;
    mid_TR = end_TR;
    mid_int = end_int;
    mid_scan = (*P).second.get_scan_number();
    P++;
  }
  
  // if contained only one peak!
  if(get_nb_ms_peaks() == 1){
    TOT_INT = start_int;
    TOT_TR = start_int * start_TR;
    TOT_SCAN = start_int * mid_scan;
  }
  // do for the last peak as well
  else{
    // compute area between mid/start:
    double area = compute_delta_area(start_TR, start_int, mid_TR, mid_int);
    // add up for apex calculation:
    TOT_TR += mid_TR * area;
    TOT_INT += area;  
    TOT_SCAN += mid_scan * area;
  }
  
  //////////////////////////////////////////////
  // compute the apex:
  // get averaged scan number:
  TOT_TR /= TOT_INT;
  TOT_SCAN /= TOT_INT;
  
  fRT = TOT_TR;   
  // set the apex ms peak:
  ms_peak* APEX = find_true_peak(TOT_SCAN);
  if ( ! APEX->getExtraPeakInfo().empty()) {
    setElutionPeakExtraInfo( APEX->getExtraPeakInfo() );
  }
  
  
  // find retention time and intensity of apex:
  fScanNumberApex = APEX->get_scan_number();
  fapex_intensity = APEX->get_intensity();

}
*/


//////////////////////////////////////////////////////////////////
// Compute a varietiy of parameters for the LC elution peak
void LC_elution_peak::computeLCElutionPeakParameters( ){
  
  
  double TOT_AREA = 0;
  double apexScan = 0;
  double apexTr = 0;
  ms_peak* endPeak = NULL;
  ms_peak* startPeak = NULL;

  // find the first peaks above the background intensity:
  SIGNAL_iterator P = get_signal_list_start();
  fScanNumberStart = (*P).second.get_scan_number();
  fStartTR = (*P).second.get_retention_time();
  
  // set start et first scan in the LC_peak:
  while(P != get_signal_list_end()){
    if( (*P).second.get_intensity() >= fSNIntensityThreshold){
      break;
    }
    P++;
  }
  
  // FLO: On windows, there is an error when we try to de-reference P when  
  // P refers to get_signal_list_end() - the case is quite obvious when we 
  // see that P is incremented afterwards. Without knowing I just introduce
  // a check whether P is equal to get_signal_list_end().
  if (P != get_signal_list_end()) {
	  startPeak = &((*P).second);
	  
	  // to compute some other parameters at the same time:
	  update_CHRGMAP( &(*P).second ); 
	  
	  P++;
  }
  
  // go through all peaks in the LC elution profile:
  while(P != get_signal_list_end()){
    
    if( (*P).second.get_intensity() >= fSNIntensityThreshold ){
      if( startPeak != NULL ){
        endPeak = &((*P).second);
      }
      else{
        startPeak = &((*P).second);
      }
    }
    else{
      endPeak = NULL;
      startPeak = NULL;
    }
    
    if( ( endPeak != NULL ) && ( startPeak != NULL ) ){
      
      // to compute some other parameters at the same time:
      update_CHRGMAP( endPeak ); 
      
      ///////////////////////////////////////////////////
      // compute an area between local start / end ms peak:
      double area = compute_delta_area( 
                                       startPeak->get_retention_time(), startPeak->get_intensity() - fSNIntensityThreshold, 
                                       endPeak->get_retention_time(), endPeak->get_intensity() - fSNIntensityThreshold
      );
      
      TOT_AREA += area;
      apexScan += (double)(P->first) * area;
      apexTr += startPeak->get_retention_time() * area;
      
      // next scan:
      startPeak = endPeak;
    }
    
    P++;
  }
  
  // if contained only one peak!
  if(get_nb_ms_peaks() == 1){
    TOT_AREA = startPeak->get_intensity();
    fScanNumberEnd = fScanNumberStart;
    fEndTR = startPeak->get_retention_time();
  }
  else{
    
    P--;
    fScanNumberEnd = (*P).second.get_scan_number();
    fEndTR = (*P).second.get_retention_time();
    fpeak_area = TOT_AREA;
    apexScan /= TOT_AREA;
    apexTr /= TOT_AREA;
    fRT = apexTr;
  }
  
  // set the apex ms peak:
  ms_peak* APEX = find_true_peak( apexScan );
  if ( ! APEX->getExtraPeakInfo().empty()) {
    setElutionPeakExtraInfo( APEX->getExtraPeakInfo() );
  }
    
  // find retention time and intensity of apex:
  fScanNumberApex = APEX->get_scan_number();
  fapex_intensity = APEX->get_intensity();

}

/*
//////////////////////////////////////////////////////////////////
// compute the total peak area:
void LC_elution_peak::compute_LC_peak_area(){
  
  //////////////////////
  double TOT_AREA = 0;
  double THRESHOLD = get_apex_intensity() * intensity_apex_percentil_cutoff;
  
  double start_TR = 0;
  double start_int = 0;
  double mid_TR = 0;  
  double mid_int = 0;
  double end_TR = 0;  
  double end_int = 0;
  fScore = 0;
  
  SIGNAL_iterator P = get_signal_list_start();
  
  // set start et first scan in the LC_peak:
  start_TR = (*P).second.get_retention_time();
  start_int = (*P).second.get_intensity();
  mid_TR = (*P).second.get_retention_time();
  mid_int = (*P).second.get_intensity();
  
  // to compute some other parameters at the same time:
  update_CHRGMAP( &(*P).second ); 
  update_SCOREMAP( &(*P).second ); 
  
  P++;
  
  /////////////////////////////////////////////////////
  // go through all peaks in the LC elution profile:
  while(P != get_signal_list_end()){
    
    if( P == get_signal_list_end())
      break;
    
    end_TR = (*P).second.get_retention_time(); 
    end_int = (*P).second.get_intensity();
    
    // to compute some other parameters at the same time:
    update_CHRGMAP( &(*P).second ); 
    update_SCOREMAP( &(*P).second ); 
    
    ///////////////////////////////////////////////////
    // compute an area around the mid ms peak:
    
    // compute area between end/mid:
    TOT_AREA += compute_delta_area( mid_TR, mid_int - THRESHOLD, end_TR, end_int - THRESHOLD);
    
    // compute area between mid/start:
    TOT_AREA += compute_delta_area(start_TR, start_int - THRESHOLD, mid_TR, mid_int - THRESHOLD);
    
    // next scan:
    start_TR = mid_TR;
    start_int = mid_int;
    mid_TR = end_TR;
    mid_int = end_int;
    
    P++;
  }
  
  // if contained only one peak!
  if(get_nb_ms_peaks() == 1){
    TOT_AREA = start_int;
    // to compute some other parameters at the same time:
    P--;
  }
  // do for the last peak as well
  else{
    // compute area between mid/start:
    TOT_AREA += compute_delta_area(start_TR, start_int, mid_TR, mid_int);
  }
  
  // set the area:
  fpeak_area = TOT_AREA;
  
  // define start end scan:
  P = get_signal_list_start();
  fScanNumberStart = (*P).second.get_scan_number();
  fStartTR = (*P).second.get_retention_time();
  
  P = get_signal_list_end();
  P--;
  fScanNumberEnd = (*P).second.get_scan_number();
  fEndTR = (*P).second.get_retention_time();
  
}
*/


/////////////////////////////////////////////////////////////////////
// compute the charge state of the LC peak
void LC_elution_peak::compute_CHRG(){
  
  bool view = false;
  double mass = get_apex_MZ();
  
  if( ( DEBUG_MASS_START <= mass) && ( mass <= DEBUG_MASS_END) ){
    view = true;
    show_info();
  }
  
  int maxCount = -1; 
  multimap< int, int>::iterator C = CHRG_MAP.begin();
  while( C != CHRG_MAP.end() ){
  
    if( view ){
      cout<<(*C).first<<":"<<(*C).second<<endl;
    }
    
    if( maxCount < (*C).second ){
      fCharge = (*C).first; 
      maxCount = (*C).second;
    }
    C++;
  }
  
  if( view ){
    cout<<fCharge<<endl;
  }
  
  CHRG_MAP.clear();
}




/////////////////////////////////////////////////////////////////////
// computes the area of between 2 peaks:
double LC_elution_peak::compute_delta_area(double START_TR, double START_INT, double END_TR, double END_INT){
  
  double AREA = 0;
  
  if( ( START_INT > 0 ) && ( END_INT > 0 ) && (START_TR <= END_TR) ){
    
    double x = (END_TR - START_TR) / TR_RESOLUTION;
    double y = fabs(END_INT - START_INT);
    
    if( ( x != 0 ) && ( y != 0 ) ){
      
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



////////////////////////////////////////////////////////////////////////////////
// print all monositopic peak cluster along the LC profile:
//void LC_elution_peak::printIsotopClusters(){
//  
//  
//  char outfile[255];
//  sprintf( outfile, "isotopCluster_%0.3f_%d", get_apex_MZ(), get_charge_state());
//  string tmp = outfile;
//  data_plotter* PLOT = new data_plotter(tmp);
//  
//  double delta = 0.002;
//    
//  // now the cluster profiles:
//  map<int, ms_peak>::iterator R = intens_signals.begin();
//  while( R != intens_signals.end() ){
//    
//    ms_peak* peak = &(*R).second;
//    map<double, double> cluster;
//    vector<CentroidPeak>::iterator p = peak->get_isotopic_peaks_start();
//    while( p != peak->get_isotopic_peaks_end() ){
//      cluster.insert( make_pair( (*p).getMass() + delta, (*p).getFittedIntensity()  ) );
//      // cluster.insert( make_pair( (*p).getMass() + delta, (*p).getFittedIntensity()  ) );
//      p++; 
//    }
//    
//    delta += delta;
//    
//    char buffer[255];
//    sprintf( buffer, "scan=%d(%0.4f)", (*R).first, peak->get_MZ());
//    PLOT->add_plot_data_impulses( &cluster , buffer);
//    
//    R++;
//  }
//  
//  PLOT->set_POINT_SIZE( 2 );
//  PLOT->plot_TWOD_data();
//  //PLOT->print_MATRIX_data_to_TXT();
//  delete PLOT;
//  PLOT = NULL;
//  
//}
//
//////////////////////////////////////////////////////////////////////////////////
//// print the consensus isotope pattern:
//void LC_elution_peak::printConsensIsotopPattern(){
//  
//  char outfile[255];
//  sprintf( outfile, "consensIsotope_%0.3f_%d", get_apex_MZ(), get_charge_state());
//  string tmp = outfile;
//  data_plotter* PLOT = new data_plotter(tmp);
//  
//  // now the cluster profiles:
//  map<double, double>::iterator I = isotopePattern->getConsensIsotopeIteratorStart();
//  while( I != isotopePattern->getConsensIsotopeIteratorEnd() ){
//    
//    char buffer[255];
//    sprintf( buffer, "%0.4f", (*I).first);
//    map<double, double> data;
//    data.insert(make_pair( (*I).first, (*I).second) );
//    PLOT->add_plot_data_impulses( &data , buffer);
//    
//    I++;
//  }
//  
//  PLOT->set_POINT_SIZE( 2 );
//  PLOT->plot_TWOD_data();
//  //PLOT->print_MATRIX_data_to_TXT();
//  delete PLOT;
//  PLOT = NULL;
//  
//}


////////////////////////////////////////////////////////////////////////////////
// print all monositopic peak cluster along the LC profile:
void LC_elution_peak::createConsensIsotopPattern(){
  
  //////////////
  // go through the different isotopes and
  // constructe a consensus patterns:
  isotopePattern = new consensIsotopePattern();
  
  multimap<int, ms_peak>::iterator R = intens_signals.begin();
  while( R != intens_signals.end() ){
    
    ms_peak* peak = &(*R).second;
    map<double, double> isotopeCluster;
    
    vector<CentroidPeak>::iterator p = peak->get_isotopic_peaks_start();
    while( p != peak->get_isotopic_peaks_end() ){
      isotopePattern->addIsotopeTrace( (*p).getMass(), (*p).getFittedIntensity()  );
      p++; 
    }
    
    R++;
  }
  
  // create the pattern:
  isotopePattern->constructConsusPattern();

}


  



////////////////////////////////////////////////////////////////////////////////
// define all required peak parameters from a single MS peak:
void LC_elution_peak::defineLCElutionPeakParametersFromMSPeak(){
  
  
  APEX = &(get_signal_list_start()->second);
  
  fMonoMass = APEX->get_MZ();
  fVolume = APEX->get_intensity();
  fCharge = APEX->get_Chrg();
  fScanNumberStart = APEX->get_Scan();
  fScanNumberApex = fScanNumberStart;
  fScanNumberEnd = fScanNumberStart;
  fapex_intensity = APEX->get_intensity();
  fRT = APEX->get_retention_time();
  fStartTR = fRT;
  fEndTR = fRT;
  fpeak_area = APEX->get_intensity();
  fSignalToNoise = APEX->getSignalToNoise();
  createConsensIsotopPattern();

  
}
  
}
