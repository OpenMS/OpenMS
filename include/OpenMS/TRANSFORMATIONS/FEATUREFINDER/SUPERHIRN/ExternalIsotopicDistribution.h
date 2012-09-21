///////////////////////////////////////////////////////////////////////////
//
//
//  by Lukas Mueller, Lukas.Mueller@imsb.biol.ethz.ch
//  October 2005
//  
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
//
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
// 
//

#ifndef _USE_EXTERNAL_ISOTOPIC_DISTRIBUTION_H
#define _USE_EXTERNAL_ISOTOPIC_DISTRIBUTION_H

#include <OpenMS/CONCEPT/Types.h>

namespace OpenMS
{

class OPENMS_DLLAPI ExternalIsotopicDistribution{

  
  ////////////////////////////////////////////////
  // declaration of the private members:
  
private:

  static std::multimap< double, PeptideIsotopeDisribution> allExternalPepdistributions;
  
  ////////////////////////////////////////////////
  // declaration of the public members:

public:
  
  static bool EXTERNAL_ISOTOPIC_PROFILES;  
  static std::string XMLInputFile;
  static double EXTERNAL_DISTRIBUTION_MONO_ISOTOPE_PPM_TOLERANCE;

  //static multimap< double, double> allExternalPepdistributions;

  // class destructor
  ~ExternalIsotopicDistribution();
  
  // class constructor
  ExternalIsotopicDistribution();
  // class copy constructor
  ExternalIsotopicDistribution(const ExternalIsotopicDistribution&);
  // class copy constructor
  ExternalIsotopicDistribution(const ExternalIsotopicDistribution*);
  
  
  //////////////////////////////////////////////////
  // overload operators:
  ExternalIsotopicDistribution& operator=(const ExternalIsotopicDistribution&);
  bool operator==(const ExternalIsotopicDistribution&);
  ExternalIsotopicDistribution& operator<=(const ExternalIsotopicDistribution&);
  ExternalIsotopicDistribution& operator>=(const ExternalIsotopicDistribution&);
  ExternalIsotopicDistribution& operator<(const ExternalIsotopicDistribution&);
  ExternalIsotopicDistribution& operator>(const ExternalIsotopicDistribution&);
  
  
  ///////////////////////////////
  // start here all the get / set
  // function to access the
  // variables of the class
  
  
  
  
  //////////////////////////////////////////////////
  // function to extract external isotopic profiles
  // by an input monoisotopic mass
  static PeptideIsotopeDisribution* extractExternalIsotopicProfile( double, int , double);
  //  check if two masses to be the same
  static bool checkMonoIsotopOverlap( double , int, double, PeptideIsotopeDisribution* );
  
  //////////////////////////////////////////////////
  // parse the external isotopic profiles
  // from xml file
  //  static void parseExternalIsotopicProfiles( );
  // extract isotopic profiles
  //  static void extractIsotopicProfiles( TiXmlDocument* );

  // init the retention time segments
  static void initRetentionTimeSegments( double, double);

  
};

} // ns

#endif
