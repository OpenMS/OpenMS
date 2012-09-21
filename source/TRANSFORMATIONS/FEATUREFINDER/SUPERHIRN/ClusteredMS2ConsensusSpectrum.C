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

#include <vector>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2Fragment.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2ConsensusSpectrum.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/ClusteredMS2ConsensusSpectrum.h>

using namespace std;


////////////////////////////////////////////////
// constructor for the object ClusteredMS2ConsensusSpectrum:
ClusteredMS2ConsensusSpectrum::ClusteredMS2ConsensusSpectrum( MS2ConsensusSpectrum* in):MS2ConsensusSpectrum(in){
  
  precursorMZ = in->getPrecursorMZ();
  TR = in->getTR();
  z = in->getPrecursorChrg();
  apexScan = in->getApexScan();
  
  this->addMS2ConsensusSpectrum( in );
}

////////////////////////////////////////////////
// constructor for the object MS2ConsensusSpectrum:
ClusteredMS2ConsensusSpectrum::ClusteredMS2ConsensusSpectrum(double iPrecursorMZ, double iTR, int iChrg, int iApexScan) : MS2ConsensusSpectrum(iPrecursorMZ, iTR, iChrg, iApexScan ){}


////////////////////////////////////////////////
// constructor for the object ClusteredMS2ConsensusSpectrum:
ClusteredMS2ConsensusSpectrum::ClusteredMS2ConsensusSpectrum( MS2Fragment* in):MS2ConsensusSpectrum(in){
  // store the scan numbers:
  MS2Scans.push_back( in->getApexScan() );
}


//////////////////////////////////////////////////
// class desctructor of ClusteredMS2ConsensusSpectrum
ClusteredMS2ConsensusSpectrum::~ClusteredMS2ConsensusSpectrum(){
  // ClusteredSpectra.clear();
}

//////////////////////////////////////////////////
// class copy constructor of ClusteredMS2ConsensusSpectrum
ClusteredMS2ConsensusSpectrum::ClusteredMS2ConsensusSpectrum(const ClusteredMS2ConsensusSpectrum& tmp):MS2ConsensusSpectrum(tmp) {
  MS2Scans = tmp.MS2Scans;
}

//////////////////////////////////////////////////
// class copy constructor of ClusteredMS2ConsensusSpectrum
ClusteredMS2ConsensusSpectrum::ClusteredMS2ConsensusSpectrum(const ClusteredMS2ConsensusSpectrum* tmp):MS2ConsensusSpectrum(tmp){
  MS2Scans = tmp->MS2Scans;
}


//////////////////////////////////////////////////
// copy constructor:
ClusteredMS2ConsensusSpectrum& ClusteredMS2ConsensusSpectrum::operator=(const ClusteredMS2ConsensusSpectrum& tmp){
  MS2ConsensusSpectrum::operator=(tmp);
  MS2Scans.clear();
  MS2Scans = tmp.MS2Scans;
  return *this;
}

//////////////////////////////////////////////////
// remove outlier fragments based on their:
// MS2Fragment::OutlierAttribute = ...
// 1: retention time
// 2: precursor mass
// etc.
void ClusteredMS2ConsensusSpectrum::removeOutlierFragments(){

  // store in this vector the fragments:
  vector< pair<double, void*> > ValueVector;
  
  /*
  // store fragments by the desired attribute as outlier detection value:
  map<double, MS2Fragment>::iterator P = MS2FragmentPeaks.begin();
  while( P != MS2FragmentPeaks.end() ){
    
    // get the attribute:
    double value = (*P).second.getOutlierDetectionAttribute();
    
    // store in the vector
    ValueVector.push_back( pair<double, void*>( value, &(*P).second ) );     
   
    P++;
  }
  
  // start the iterative outlier removal by the set attribut of a ms2 fragment:
  simple_math myMath;
  myMath.ITERATIVE_OUTLIER_DETECTION_BY_DIXON( &ValueVector );

  // convert back after the oulier detected vector:
  map<double, MS2Fragment> newFragments;
  vector< pair<double, void*> >::iterator I = ValueVector.begin();
  while( I != ValueVector.end() ){
    pair<double, void*> p = (*I);
    MS2Fragment* frag = (MS2Fragment*) p.second;
    newFragments.insert( make_pair( frag->getFragmentPeakArea(), *frag ) );
    I++;
  }
  
  // copy back:
  MS2FragmentPeaks.clear();
  MS2FragmentPeaks = newFragments;
  newFragments.clear();
   */
  
}


//////////////////////////////////////////////////
// trace the fragments across MS/MS spectra:
void ClusteredMS2ConsensusSpectrum::constructClusteredConsenusSpectraFragments(MS2ConsensusSpectrum* newMSMS){
  
  /*
  // stores the individual consens spectra:
  map<int, MS2ConsensusSpectrum* >::iterator C = ConsensSpectra.begin();
  while( C != ConsensSpectra.end() ){
    extractFragmentsFromConsensSpectra( (*C).first, (*C).second );
    C++; 
  }
  
  // process extracted aligned fragment traces:
  processAlignedFragmentTraces();
  */
}



//////////////////////////////////////////////////
// extracts fragments from a MS/MS spectra and inserts
// them into the Clustered MS/MS spectrum:
void ClusteredMS2ConsensusSpectrum::extractFragmentsFromSpectra( MS2ConsensusSpectrum* in  ){
  
  
  // go through the MS2 fragments and find the common one:
  multimap<double, MS2Fragment>::iterator P = in->getMS2FragmentPeakStart();
  while( P != in->getMS2FragmentPeakEnd() ){
    
    // fragment search mass:
    MS2Fragment* frag = &(P->second);
    double searchMZ = frag->getFragmentMz();
    
    MS2Fragment* matchedFrag = this->findMS2Fragment( searchMZ );    
    
    // in case where the fragment has been detected in a previous MS/MS
    // merge thes in one MS2 fragment
    if( matchedFrag != NULL ){
      this->mergeMS2Fragments( matchedFrag, frag );
    }
    else{
      this->addMS2Fragment( frag );
    }
    
    P++;
  }
}






//////////////////////////////////////////////////
// add a MS2 fragment:
void ClusteredMS2ConsensusSpectrum::addMS2ConsensusSpectrum( MS2ConsensusSpectrum* in ){
  
  // extract directly the MS/MS clustered spectrum
  extractFragmentsFromSpectra( in );
  
  // store the scan numbers:
  MS2Scans.push_back( in->getApexScan() );
  
  
}


//////////////////////////////////////////////////////
// merge a MS2 fragment into the target MS2 fragment:
void ClusteredMS2ConsensusSpectrum::mergeMS2Fragments(MS2Fragment* target, MS2Fragment* toMerge ){
  
  // sum up intensities:
  target->setFragmentPeakArea( target->getFragmentPeakArea() + toMerge->getFragmentPeakArea());
  
  // average m/z:
  target->setFragmentMz( ( target->getFragmentMz() + toMerge->getFragmentMz() ) / 2.0);
  
  // average m/z:
  target->setPrecursorMZ( ( target->getPrecursorMZ() + toMerge->getPrecursorMZ() ) / 2.0);
  
  
}

//////////////////////////////////////////////////////
// plot all the consensus MS2 spectrum in one plot:
void ClusteredMS2ConsensusSpectrum::plotCombinedSpectra( ){
  /*
  
  char buffer[255];
  sprintf( buffer, "CombinedMS2ConsSpec%0.2f", precursorMZ);
  string tmp = buffer;
  data_plotter* PLOT = new data_plotter( tmp ); 
  
  map<int, MS2ConsensusSpectrum* >::iterator C = ConsensSpectra.begin();
  while( C != ConsensSpectra.end() ){

    sprintf( buffer, "LC-MS %d", C->first);
    
    map<double, double> data;
    map<double, MS2Fragment>::iterator I = C->second->getMS2FragmentPeakStart();
    while( I != C->second->getMS2FragmentPeakEnd() ){
      data.insert( make_pair( I->second.getFragmentMz(), I->second.getFragmentPeakArea()) );
      I++;
    }
    
    PLOT->add_plot_data_impulses( &data, buffer);
  
    data.clear();
    C++;
  }
  
  
  PLOT->plot_TWOD_data();
  
  delete PLOT;
  PLOT = NULL;
  */
  
}

