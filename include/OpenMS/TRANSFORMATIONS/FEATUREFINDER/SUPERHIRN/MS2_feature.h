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
// variable description:
//
//
// function description:
//
//
// **********************************************************************//


#ifndef MS2_FEATURE_H
#define MS2_FEATURE_H

#include "MS2Fragment.h"
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/ClusteredMS2ConsensusSpectrum.h>


class MS2_feature : public ClusteredMS2ConsensusSpectrum {

  
  using ClusteredMS2ConsensusSpectrum::operator=;
    
  ////////////////////////////////////////////////
  // declaration of the private members:
  
private:
  
  int ID;
  
  ////////////////////////////////////////////////
  // declaration of the public members:
  
public:
  
  
  // class destructor
  ~MS2_feature();
  // class constructor
  MS2_feature();
  MS2_feature(MS2Fragment*);
  MS2_feature(double iPrecursorMZ, double iTR, int iChrg, int iApexScan);
  
  // class copy constructor
  MS2_feature(const MS2_feature&);
  // class copy constructor
  MS2_feature(const MS2_feature*);
  
  
  //////////////////////////////////////////////////
  // overload operators:
  MS2_feature& operator=(const MS2_feature&);
  bool operator==(const MS2_feature&);
  MS2_feature& operator<=(const MS2_feature&);
  MS2_feature& operator>=(const MS2_feature&);
  MS2_feature& operator<(const MS2_feature&);
  MS2_feature& operator>(const MS2_feature&);
  
  
  // show info 
  void show_info();

  
  ///////////////////////////////
  // start here all the get / set
  // function to access the
  // variables of the class

  void setID( int in ){ID = in;};
  int getID( ){return ID;};

  
};

#endif

    
