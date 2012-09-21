///////////////////////////////////////////////////////////////////////////
//
//  written by Lukas N Mueller, 30.3.05
//  Lukas.Mueller@imsb.biol.ethz.ch
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
//  
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
//
//
// **********************************************************************//
// CLASS featureLCprofile:
//
// variable description:
//
//
// function description:
//
//
// **********************************************************************//


#ifndef _FEATURE_LC_PROFILE
#define _FEATURE_LC_PROFILE_H

#include <OpenMS/CONCEPT/Types.h>

#include <map>

namespace OpenMS
{

// this structure stores the monoisotopic
// signals of LC elution peak:
struct MS1Signal{
  double mass;
  double TR;
  double intensity;
  int scan;
  int charge;
};



class OPENMS_DLLAPI featureLCprofile{

    
  ////////////////////////////////////////////////
  // declaration of the private members:
private:
  
  std::map<int, MS1Signal> LCelutionSignals;
  std::map<int, MS1Signal> outsideLCelutionSignals;

  
  // elution area
  double LCelutionArea;
  
  // apex signale:
  MS1Signal apexMS1Signal;
  
  ////////////////////////////////////////////////
  // declaration of the public members:
  
public:
  
  // class destructor
  ~featureLCprofile();
  
  // class constructor
  featureLCprofile();
  featureLCprofile(double , double, int, double );
  featureLCprofile(double , double , double , int , int, double );
  
  // class copy constructor
  featureLCprofile(const featureLCprofile&);
  // class copy constructor
  featureLCprofile(const featureLCprofile*);
  
  
  //////////////////////////////////////////////////
  // overload operators:
  featureLCprofile& operator=(const featureLCprofile&);
  bool operator==(const featureLCprofile&);
  featureLCprofile& operator<=(const featureLCprofile&);
  featureLCprofile& operator>=(const featureLCprofile&);
  featureLCprofile& operator<(const featureLCprofile&);
  featureLCprofile& operator>(const featureLCprofile&);
  
  // change all elution times by a factor:
  void changeElutionTimesByFactor(double);
  
  
  ///////////////////////////////
  // start here all the get / set
  // function to access the
  // variables of the class

  // add / get signals:
  void addMS1elutionSignal(  double , double , int , int, double );
  void addOutsideMS1elutionSignal(  double , double , int , int, double );
  void addMS1elutionSignal( MS1Signal* in);

  std::map<int, MS1Signal>* getLCelutionSignalMap(){ return &LCelutionSignals;};
  std::map<int, MS1Signal>::iterator getLCelutionSignalsStart(){ return LCelutionSignals.begin();};
  std::map<int, MS1Signal>::reverse_iterator getLastLCelutionSignal(){ return LCelutionSignals.rbegin();};
  std::map<int, MS1Signal>::iterator getLCelutionSignalsEnd(){ return LCelutionSignals.end();};
  int getNbLCelutionSignals(){ return LCelutionSignals.size();};

  
};

} // ns

#endif

    
