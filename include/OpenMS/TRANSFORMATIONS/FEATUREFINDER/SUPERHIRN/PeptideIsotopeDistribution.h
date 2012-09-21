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

#include <string>
#include <vector>

using namespace std;

class PeptideIsotopeDisribution{

    
  ////////////////////////////////////////////////
  // declaration of the private members:
  
private:
  
  vector<double> mass;
  vector<double> intens;
  double* intensArray;
  
  string sq;
  string name;
  string summary;
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
  PeptideIsotopeDisribution(vector<double>, vector<double>);
  PeptideIsotopeDisribution(vector<double>, vector<double>, int, string, string, int);
  PeptideIsotopeDisribution(vector<double>, vector<double>, int, string, string, int, double);
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
  string getIsotopeDistInfo();


  ///////////////////////////////
  // start here all the get / set
  // function to access the
  // variables of the class

  // get the intensity values in form of an array
  double* getIntensityArray( );
  string getName(){return name;};
  string getSequence(){return sq;};
  int getChargeState(){return chargeState;};
  double getMonoMass(){ return mass[0];};
  
  int getID(){return id;};
  string getSummary(){return summary;};
  double getRTSegment(){return RtSegment;};
  
  
  void setRTStart(double in){RtStart=in;};
  double getRTStart(){return RtStart;};
  void setRTEnd(double in){RtEnd=in;};
  double getRTEnd(){return RtEnd;};

};

#endif

    
