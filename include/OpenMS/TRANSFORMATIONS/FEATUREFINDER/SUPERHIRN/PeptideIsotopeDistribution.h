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


#ifndef PEPTIDE_ISOTOPE_DISTRIBUTION_H
#define PEPTIDE_ISOTOPE_DISTRIBUTION_H

#include <OpenMS/CONCEPT/Types.h>

#include <string>
#include <vector>

namespace OpenMS
{

class OPENMS_DLLAPI PeptideIsotopeDisribution{

    
  ////////////////////////////////////////////////
  // declaration of the private members:
  
private:
  
  std::vector<double> mass;
  std::vector<double> intens;
  double* intensArray;
  
  std::string sq;
  std::string name;
  std::string summary;
  int chargeState;
  int id;
  double RtSegment;
  double RtEnd;
  double RtStart;
  
  ////////////////////////////////////////////////
  // declaration of the public members:

public:
  
  // class destructor
  ~PeptideIsotopeDisribution();
  
  // class constructor
  PeptideIsotopeDisribution();
  PeptideIsotopeDisribution(std::vector<double>, std::vector<double>);
  PeptideIsotopeDisribution(std::vector<double>, std::vector<double>, int, std::string, std::string, int);
  PeptideIsotopeDisribution(std::vector<double>, std::vector<double>, int, std::string, std::string, int, double);
  // class copy constructor
  PeptideIsotopeDisribution(const PeptideIsotopeDisribution&);
  // class copy constructor
  PeptideIsotopeDisribution(const PeptideIsotopeDisribution*);
  
  // show info:
  void show_info();

  // construct info in for or a string about the PeptideIsotopeDisribution:
  void constructSummaryString();
  
  
  //////////////////////////////////////////////////
  // overload operators:
  PeptideIsotopeDisribution& operator=(const PeptideIsotopeDisribution&);
  bool operator==(const PeptideIsotopeDisribution&);
  PeptideIsotopeDisribution& operator<=(const PeptideIsotopeDisribution&);
  PeptideIsotopeDisribution& operator>=(const PeptideIsotopeDisribution&);
  PeptideIsotopeDisribution& operator<(const PeptideIsotopeDisribution&);
  PeptideIsotopeDisribution& operator>(const PeptideIsotopeDisribution&);

  // return info in form or a string about PeptideIsotopeDisribution:
  std::string getIsotopeDistInfo();


  ///////////////////////////////
  // start here all the get / set
  // function to access the
  // variables of the class

  // get the intensity values in form of an array
  double* getIntensityArray( );
  std::string getName(){return name;};
  std::string getSequence(){return sq;};
  int getChargeState(){return chargeState;};
  double getMonoMass(){ return mass[0];};
  
  int getID(){return id;};
  std::string getSummary(){return summary;};
  double getRTSegment(){return RtSegment;};
  
  
  void setRTStart(double in){RtStart=in;};
  double getRTStart(){return RtStart;};
  void setRTEnd(double in){RtEnd=in;};
  double getRTEnd(){return RtEnd;};

};

} // ns

#endif

    
