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

#include <map>
#include <vector>
#include <string>
#include <math.h>
#include <stdio.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/simple_math2.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2Fragment.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2ConsensusSpectrum.h>

namespace OpenMS
{

using namespace std;

// mass to charge tolerance for MS2 trace level:
double MS2ConsensusSpectrum::MS2_MZ_TOLERANCE;



////////////////////////////////////////////////
// constructor for the object MS2ConsensusSpectrum:
MS2ConsensusSpectrum::MS2ConsensusSpectrum(){
}

////////////////////////////////////////////////
// constructor for the object MS2ConsensusSpectrum:
MS2ConsensusSpectrum::MS2ConsensusSpectrum(MS2Fragment* in){  
  addMS2Fragment( in );
}

////////////////////////////////////////////////
// constructor for the object MS2ConsensusSpectrum:
MS2ConsensusSpectrum::MS2ConsensusSpectrum(double iPrecursorMZ, double iTR, int iChrg, int iApexScan){
  
  precursorMZ = iPrecursorMZ;
  TR = iTR;
  startTR = TR;
  endTR = TR;
  z = iChrg;
  apexScan = iApexScan;
  
}



//////////////////////////////////////////////////
// class desctructor of MS2ConsensusSpectrum
MS2ConsensusSpectrum::~MS2ConsensusSpectrum(){
  MS2FragmentPeaks.clear();
}

//////////////////////////////////////////////////
// class copy constructor of MS2ConsensusSpectrum
MS2ConsensusSpectrum::MS2ConsensusSpectrum(const MS2ConsensusSpectrum& tmp){
  TR = tmp.TR;
  startTR = tmp.startTR;
  endTR = tmp.endTR;
  z = tmp.z;
  apexScan = tmp.apexScan;
  startScan = tmp.startScan;
  endScan = tmp.endScan;
  precursorMZ = tmp.precursorMZ;
  MS2FragmentPeaks.clear();
  MS2FragmentPeaks = tmp.MS2FragmentPeaks;
}

//////////////////////////////////////////////////
// class copy constructor of MS2ConsensusSpectrum
MS2ConsensusSpectrum::MS2ConsensusSpectrum(const MS2ConsensusSpectrum* tmp){
  TR = tmp->TR;
  startTR = tmp->startTR;
  endTR = tmp->endTR;
  z = tmp->z;
  apexScan = tmp->apexScan;
  startScan = tmp->startScan;
  endScan = tmp->endScan;
  precursorMZ = tmp->precursorMZ;
  MS2FragmentPeaks.clear();
  MS2FragmentPeaks = tmp->MS2FragmentPeaks;
}


//////////////////////////////////////////////////
// copy constructor:
MS2ConsensusSpectrum& MS2ConsensusSpectrum::operator=(const MS2ConsensusSpectrum& tmp){
  TR = tmp.TR;
  startTR = tmp.startTR;
  endTR = tmp.endTR;
  z = tmp.z;
  apexScan = tmp.apexScan;
  startScan = tmp.startScan;
  endScan = tmp.endScan;
  precursorMZ = tmp.precursorMZ;
  MS2FragmentPeaks.clear();
  MS2FragmentPeaks = tmp.MS2FragmentPeaks;
  return *this;
}

//////////////////////////////////////////////////
// remove outlier fragments based on their:
// MS2Fragment::OutlierAttribute = ...
// 1: retention time
// 2: precursor mass
// etc.
void MS2ConsensusSpectrum::removeOutlierFragments(){

  // store in this vector the fragments:
  vector< pair<double, void*> > ValueVector;
  
  // store fragments by the desired attribute as outlier detection value:
  multimap<double, MS2Fragment>::iterator P = MS2FragmentPeaks.begin();
  while( P != MS2FragmentPeaks.end() ){
    
    // get the attribute:
    double value = (*P).second.getOutlierDetectionAttribute();
    // store in the vector
    ValueVector.push_back( pair<double, void*>( value, &(*P).second ) );     
   
    P++;
  }
  
  // start the iterative outlier removal by the set attribut of a ms2 fragment:
  simple_math2 myMath;
  myMath.ITERATIVE_OUTLIER_DETECTION_BY_DIXON( &ValueVector );

  // convert back after the oulier detected vector:
  multimap<double, MS2Fragment> newFragments;
  vector< pair<double, void*> >::iterator I = ValueVector.begin();
  while( I != ValueVector.end() ){
    pair<double, void*> p = (*I);
    MS2Fragment* frag = (MS2Fragment*) p.second;
    newFragments.insert( make_pair( frag->getFragmentMz(), *frag ) );
    I++;
  }
  
  // copy back:
  MS2FragmentPeaks.clear();
  MS2FragmentPeaks = newFragments;
  newFragments.clear();
  
}


//////////////////////////////////////////////////
// process the stored fragments:
void MS2ConsensusSpectrum::processConsenusSpectraFragments(){
  
  
  if( MS2FragmentPeaks.size() > 1 ){
    //////////////////////////
    // remove outlier fragments based on their:
    // MS2Fragment::OutlierAttribute = ...
    // 1: retention time
    // 2: precursor mass
    // etc.
    MS2Fragment::OutlierAttribute = 1;
    removeOutlierFragments();
    
    
    /*
    //////////////////////////////////////////
    // compute now average values of the MS2 consensus spectrum:
    // - precursor mz
    // - apex TR
    // - apex scan
    // - charge state:
    
    // convert back after the oulier detected vector:
    vector<double> TRs;
    vector<double> MZs;
    vector<double> Scans;
    vector<double> Zs;
    map< double, MS2Fragment>::iterator I = MS2FragmentPeaks.begin();
    while( I != MS2FragmentPeaks.end() ){
      TRs.push_back( (*I).second.getTR() );
      MZs.push_back( (*I).second.getPrecursorMZ() );
      Scans.push_back( (double)(*I).second.getApexScan() );
      Zs.push_back( (double)(*I).second.getCHRG() );
      I++;
    }
    
    // recompute the retention time:
    simple_math myMath;
    TR = myMath.MEDIAN( &TRs );
    precursorMZ = myMath.MEDIAN( &MZs );
    apexScan = myMath.MEDIAN( &Scans );
    z = myMath.MEDIAN( &Zs );
    */
    
    computeMS2SpectrumParameters();
  }
 }



//////////////////////////////////////////////////
// compute MS2 parameters
void MS2ConsensusSpectrum::computeMS2SpectrumParameters(){
  
  if( MS2FragmentPeaks.size() > 1 ){
    
    
    
    double totArea = 0;
    TR = 0;
    startTR = 0;
    endTR = 0;
    precursorMZ = 0;
    
    double iz = 0;
    double iapexScan = 0;
    double istartScan = 0;
    double iendScan = 0;
    
    multimap< double, MS2Fragment>::iterator I = MS2FragmentPeaks.begin();
    while( I != MS2FragmentPeaks.end() ){
      
      double thisArea = (*I).second.getFragmentPeakArea();
      totArea += thisArea;
      TR += thisArea * (*I).second.getTR();
      startTR += thisArea * (*I).second.getStartTR();      
      endTR += thisArea * (*I).second.getEndTR();
      precursorMZ += thisArea * (*I).second.getPrecursorMZ();
      istartScan += thisArea * (*I).second.getStartScan();
      iendScan += thisArea * (*I).second.getEndScan();
      iapexScan += thisArea * (*I).second.getApexScan();
      iz += thisArea * (*I).second.getCHRG();
            
      I++;
    }

    TR /= totArea;
    startTR /= totArea;
    endTR /= totArea;
    precursorMZ /= totArea;
    
    istartScan /= totArea;
    startScan = (int)istartScan;
    
    iendScan /= totArea;
    endScan = (int)iendScan;
    
    iz /= totArea;
    z = (int)iz;
    
    iapexScan /= totArea;
    apexScan = (int)iapexScan;
    
    
    // rearange 

  }
  else{
    
  
    MS2Fragment* in = &(MS2FragmentPeaks.begin()->second);
    
    // start / end scane:
    startScan = in->getStartScan();
    endScan = in->getEndScan();
    
    // start /end TR borders:
    startTR = in->getStartTR();
    endTR = in->getEndTR();
    
    // precursor and apex:
    precursorMZ = in->getPrecursorMZ();
    TR = in->getTR();
    z = in->getCHRG();
    apexScan = in->getApexScan();

  }
}


//////////////////////////////////////////////////
// add a MS2 fragment:
void MS2ConsensusSpectrum::addMS2Fragment( MS2Fragment* in ){
  
  // make a mz map:
  MS2FragmentPeaks.insert( make_pair( in->getFragmentMz(), *in ) );
  
  // compute online the average retention time
  // and the precursor mass:
  computeMS2SpectrumParameters();
  
}

//////////////////////////////////////////////////////
// plot the consensus MS2 spectrum:
void MS2ConsensusSpectrum::plotSpectrum( ){

  char buffer[255];
  sprintf( buffer, "MS2ConsSpec%0.2f", precursorMZ);
  plotSpectrum( buffer );
  
}

//////////////////////////////////////////////////////
// plot the consensus MS2 spectrum:
void MS2ConsensusSpectrum::plotSpectrum( string fileName ){
  
 // deleted
  
}

//////////////////////////////////////////////////////
// show MS2 spectrum info:
void MS2ConsensusSpectrum::show_info( ){
  
  printf( "\tMS2 consenus spectrum: m/z=%0.3f,Tr=%0.2f,scan=%d,#fragments=%d\n",
          precursorMZ, TR, apexScan, MS2FragmentPeaks.size());

  /*
  map<double, MS2Fragment>::iterator P = MS2FragmentPeaks.begin();
  while( P != MS2FragmentPeaks.end() ){
    P->second.show_info();
    P++;
  }
   */
  
  
  
}




////////////////////////////////////////////////////////
// find a corresponding MS2 fragment
MS2Fragment* MS2ConsensusSpectrum::findMS2Fragment( double mass ){
  
  ///////////////////////
  // collect a list of iterators with potential candidates:
  map< double, multimap< double, MS2Fragment >::iterator > candidates;
  
  // scan lower mass tolerance:
  multimap<double, MS2Fragment>::iterator F = MS2FragmentPeaks.lower_bound( mass );
  multimap< double, MS2Fragment >::iterator I = F;
  if( I != MS2FragmentPeaks.begin() ){
    I--;
  }
  
  while( simple_math2::compareMassValuesAtPPMLevel(I->second.getFragmentMz(), mass, MS2ConsensusSpectrum::MS2_MZ_TOLERANCE ) ){
    
    candidates.insert( make_pair( fabs( I->second.getFragmentMz() - mass ), I) );
    if( I == MS2FragmentPeaks.begin() ){
      break;
    }
    
    // next:
    I--;
  }
  
  
  // scan upper mass tolerance:
  I = F;
  if( ( I != MS2FragmentPeaks.end() ) && ( I != MS2FragmentPeaks.begin() ) ){

    while( simple_math2::compareMassValuesAtPPMLevel(I->second.getFragmentMz(), mass, MS2ConsensusSpectrum::MS2_MZ_TOLERANCE ) ){
      
      candidates.insert( make_pair( fabs( I->second.getFragmentMz() - mass ), I) );
      I++;
      if( I == MS2FragmentPeaks.end() ){
        break;
      }
    }
  }
  
  
  /////////////////////////////////////////////////////
  // find now the one with the best match:
  // i.e. take the one with the smallest Mz difference:
  if( !candidates.empty() ){
    return &((candidates.begin())->second->second);
  }

  return NULL;
  
}


//////////////////////////////////////////////////
// remove H2O loss region of the MS2 spectra
void MS2ConsensusSpectrum::removeWaterLossRegion( ){
  
  // define water loss region:
  // 3 times H2O / by charge state:
  // max = precursor mass:
  double minLossMZRegion = precursorMZ - 30;
  double maxLossMZRegion = precursorMZ;
    
  multimap<double, MS2Fragment>::iterator I = getMS2FragmentPeakStart();  
  while( I != getMS2FragmentPeakEnd() ) {
    
    
    if( ( I->second.getFragmentMz() >= minLossMZRegion ) && ( I->second.getFragmentMz() < maxLossMZRegion ) ){
      getMS2FragmentMap()->erase( I++ );
    }
    else{
      I++;
    }
  
  }  
  
  
}



//////////////////////////////////////////////////////
// copmute the similarity of the elution shape of the
// MS2 fragment to this MS2 consensus spectrum
double MS2ConsensusSpectrum::getLCElutionPeakSimilarity( MS2Fragment* frag ){
  
  double startTR = frag->getStartTR();
  if( startTR > getStartTR() ){
    startTR = getStartTR();
  }
  
  double totLCSpec = getEndTR() - startTR;
  double startLCSpec = getTR() - startTR;
  double corSpec = startLCSpec / totLCSpec;

  double totLCMS2 = frag->getEndTR() - startTR;
  double startLCMS2 = frag->getTR() - startTR;
  double corMS2 = startLCMS2 / totLCMS2;

  ///////////
  double av = fabs(getEndTR() - frag->getEndTR());
  av += fabs(getTR() - frag->getTR() );
  av += fabs(getStartTR() - frag->getStartTR() );
  return av;
  
  return corMS2/corSpec;
}

}
