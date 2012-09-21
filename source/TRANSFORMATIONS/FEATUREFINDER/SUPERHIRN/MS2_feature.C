///////////////////////////////////////////////////////////////////////////
//
//  written by Lukas N Mueller, 30.3.05
//  Lukas.Mueller@imsb.biol.ethz.ch
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
//
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
//
// **********************************************************************//
// CLASS MS2_feature:
//
//
// function description:
//
//
// **********************************************************************//

#include <stdio.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2Fragment.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/ClusteredMS2ConsensusSpectrum.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2_feature.h>


////////////////////////////////////////////////
// constructor for the object MS2_feature:
MS2_feature::MS2_feature(MS2Fragment* in):ClusteredMS2ConsensusSpectrum(in){
  ID = -1;
}

////////////////////////////////////////////////
// constructor for the object MS2ConsensusSpectrum:
MS2_feature::MS2_feature(double iPrecursorMZ, double iTR, int iChrg, int iApexScan) : ClusteredMS2ConsensusSpectrum(iPrecursorMZ, iTR, iChrg, iApexScan ){
  ID = -1;
}

//////////////////////////////////////////////////
// class desctructor of MS2_feature
MS2_feature::~MS2_feature(){
}

//////////////////////////////////////////////////
// class copy constructor of MS2_feature
MS2_feature::MS2_feature(const MS2_feature& tmp):ClusteredMS2ConsensusSpectrum(tmp){
  ID = tmp.ID;
}

//////////////////////////////////////////////////
// class copy constructor of MS2_feature
MS2_feature::MS2_feature(const MS2_feature* tmp):ClusteredMS2ConsensusSpectrum(tmp){
  ID = tmp->ID;
}

//////////////////////////////////////////////////
// copy constructor:
MS2_feature& MS2_feature::operator=(const MS2_feature& tmp){
  ClusteredMS2ConsensusSpectrum::operator=(tmp);
  ID = tmp.ID;
  return *this;
}


/////////////////////////////////////////////
// show info 
void MS2_feature::show_info(){
  
  printf("DELETED");
  
  //((ClusteredMS2ConsensusSpectrum)this).show_info();
}

    