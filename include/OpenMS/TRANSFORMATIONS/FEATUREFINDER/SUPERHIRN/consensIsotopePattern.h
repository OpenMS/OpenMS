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
// CLASS consensIsotopePattern:
//
// variable description:
//
//
// function description:
//
//
// **********************************************************************//


#ifndef CONSENS_ISOTOPE_PATTERN_H
#define CONSENS_ISOTOPE_PATTERN_H

#include <map>
#include <vector>

class consensIsotopePattern{

  ////////////////////////////////////////////////
  // declaration of the private members:
  
private:
  
  // stores teh consenus pattern:
  std::map<  double , double > isotopesTrace; 
  std::vector<  double > mzIsotopesStDev; 
  std::vector<  double > intensIsotopesStDev; 
  
  // stores the detected patterns by retention time
  std::map<double, std::pair< std::vector<double>, std::vector<double> > > rawIsotopes;
  
  ////////////////////////////////////////////////
  // declaration of the public members:
  
  
public:
    
    static double FT_MZ_TOLERANCE;
  
  // class destructor
  ~consensIsotopePattern();
  
  // class constructor
  consensIsotopePattern();
  // class copy constructor
  consensIsotopePattern(const consensIsotopePattern&);
  // class copy constructor
  consensIsotopePattern(const consensIsotopePattern*);
  
  
  //////////////////////////////////////////////////
  // overload operators:
  consensIsotopePattern& operator=(const consensIsotopePattern&);
  bool operator==(const consensIsotopePattern&);
  consensIsotopePattern& operator<=(const consensIsotopePattern&);
  consensIsotopePattern& operator>=(const consensIsotopePattern&);
  consensIsotopePattern& operator<(const consensIsotopePattern&);
  consensIsotopePattern& operator>(const consensIsotopePattern&);
  
  
  // construc the consenus pattern:
  void constructConsusPattern( );
  // order an isotope trace in the correct cluster:
  void addIsotopeTrace( double, double );
  // condens the pattern, make averge peaks from the traces:
  void condensIsotopePattern( std::pair< std::vector<double>, std::vector<double> >*);

  
  ///////////////////////////////
  // start here all the get / set
  // function to access the
  // variables of the class

  std::map<double, double>::iterator getConsensIsotopeIteratorStart(){return isotopesTrace.begin();}; 
  std::map<double, double>::iterator getConsensIsotopeIteratorEnd(){return isotopesTrace.end();}; 


};

#endif

    
