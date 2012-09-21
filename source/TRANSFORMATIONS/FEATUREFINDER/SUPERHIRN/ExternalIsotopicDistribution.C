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

#include <string>
#include <vector>
#include <math.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/PeptideIsotopeDistribution.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/ExternalIsotopicDistribution.h>

using namespace std;

bool ExternalIsotopicDistribution::EXTERNAL_ISOTOPIC_PROFILES = false;
string ExternalIsotopicDistribution::XMLInputFile;
double ExternalIsotopicDistribution::EXTERNAL_DISTRIBUTION_MONO_ISOTOPE_PPM_TOLERANCE;
multimap< double, PeptideIsotopeDisribution> ExternalIsotopicDistribution::allExternalPepdistributions;



////////////////////////////////////////////////
// constructor for the object ExternalIsotopicDistribution:
ExternalIsotopicDistribution::ExternalIsotopicDistribution(){
}

//////////////////////////////////////////////////
// class desctructor of ExternalIsotopicDistribution
ExternalIsotopicDistribution::~ExternalIsotopicDistribution(){
}

//////////////////////////////////////////////////
// class copy constructor of ExternalIsotopicDistribution
ExternalIsotopicDistribution::ExternalIsotopicDistribution(const ExternalIsotopicDistribution& tmp){
}

//////////////////////////////////////////////////
// class copy constructor of ExternalIsotopicDistribution
ExternalIsotopicDistribution::ExternalIsotopicDistribution(const ExternalIsotopicDistribution* tmp){
}


//////////////////////////////////////////////////
// copy constructor:
ExternalIsotopicDistribution& ExternalIsotopicDistribution::operator=(const ExternalIsotopicDistribution& tmp){
  return *this;
}




//////////////////////////////////////////////////
// function to extract external isotopic profiles
// by an input monoisotopic mass
PeptideIsotopeDisribution* ExternalIsotopicDistribution::extractExternalIsotopicProfile( double monoMass, int charge, double RT){

  multimap<double, PeptideIsotopeDisribution>::iterator F = allExternalPepdistributions.lower_bound( monoMass );
  multimap<double, PeptideIsotopeDisribution>::iterator P = F;
  
  // search down:
  do{
  
    if( checkMonoIsotopOverlap( monoMass, charge, RT, &(P->second) ) ){
      return &( P->second );
    }

    P--;    
  }
  while( P != allExternalPepdistributions.begin() );
  
  
  // search up:
  while( F != allExternalPepdistributions.end() ){
    if( checkMonoIsotopOverlap( monoMass, charge, RT, &(P->second) ) ){
      return &( P->second );
    }
    F++;
  }

  return NULL;

}

// COPIED FROM SIMPLE_MATH
//////////////////////////////////////////////////////
// compare to masses at the PPM value and decided
// if they fall into the m/z tolerance window
bool simple_math_compareMassValuesAtPPMLevel( double mzA, double mzB, double PPM_TOLERANCE ){
  
  // take the average mass:
  double avMass = (mzA + mzB) / 2.0;
  
  // define the parts per million:
  double ppmValue = avMass / 1000000.00;
  double ppmDeltaTol = ppmValue * PPM_TOLERANCE;
  
  double deltaMass = fabs( mzA - mzB);
  if( deltaMass > ppmDeltaTol ){
    return false;
  }
  
  return true;
}


//////////////////////////////////////////////////
//  check if two masses to be the same
bool ExternalIsotopicDistribution::checkMonoIsotopOverlap( double MeasMonoMass, int TestCharge, double TR, PeptideIsotopeDisribution* dist ){
  
  // check the charge state:
  if( TestCharge != dist->getChargeState()){
    return false;
  }
  
  // check the mass tolerance:
  if( !simple_math_compareMassValuesAtPPMLevel( dist->getMonoMass(), MeasMonoMass, ExternalIsotopicDistribution::EXTERNAL_DISTRIBUTION_MONO_ISOTOPE_PPM_TOLERANCE  ) ){
    return false;
  }
  
  if( ( TR < dist->getRTStart() ) || ( TR > dist->getRTEnd() ) ){
    return false;
  }

  return true;  
}



//////////////////////////////////////////////////
// parse the external isotopic profiles
// from xml file
//void ExternalIsotopicDistribution::parseExternalIsotopicProfiles( ){
//  
//  TiXmlDocument* PARSER = new TiXmlDocument( XMLInputFile.c_str() );
//  if( PARSER->LoadFile() ){
//    printf("\t\t\t parsing External Isotopic Distribution XML file '%s' ...\n", XMLInputFile.c_str() );    
//  }
//  else{    
//    
//    if( PARSER != NULL ){
//      delete PARSER;
//      PARSER = NULL;
//    }
//    
//    if( !XMLInputFile.empty() ){
//      printf("\n\t** Could not open External Isotopic Distribution XML file %s ** \n", XMLInputFile.c_str() );
//    }
//  }
//  
//  // parsing was successful, then go an extract the patterns
//  if( PARSER != NULL ){
//    extractIsotopicProfiles( PARSER );
//    
//    if( ! allExternalPepdistributions.empty() ){
//      ExternalIsotopicDistribution::EXTERNAL_ISOTOPIC_PROFILES = true;
//    }
//    
//  }
//  
//  delete PARSER;
//  PARSER = NULL;
//
//  
//}


//////////////////////////////////////////////////
// extract isotopic profiles
//void ExternalIsotopicDistribution::extractIsotopicProfiles( TiXmlDocument* PARSER ){
//  
//  
//  // get the element of the ms search summary and go through each of them:
//  TiXmlNode* IsotopProfiles = PARSER->FirstChild( "PeptidePins" ); 
//  
//  TiXmlNode* profile;
//  for( profile = IsotopProfiles->FirstChild("PeptidePin"); profile; profile = profile->NextSibling() ){
//    
//    TiXmlElement* profileE = profile->ToElement();
//    int id = atoi( profileE->Attribute("id") );
//    string pinName = profileE->Attribute("name");        
//    string pinSQ = profileE->Attribute("peptideSequence");        
//    double rtSegment = atof( profileE->Attribute("rt_segment") );        
//        
//    if( profile != NULL ){
//      
//      TiXmlNode* myIsotopes = profile->FirstChild( "Isotopes" ); 
//      TiXmlElement* isotopes_E = myIsotopes->ToElement();
//      double chargeState = atof( isotopes_E->Attribute("charge") );        
//      vector<double> masses;
//      vector<double> intens;
//      
//      double totIntens = 0;
//      
//      TiXmlNode* isotope;
//      for( isotope = myIsotopes->FirstChild("Isotope"); isotope; isotope = isotope->NextSibling() ){
//        TiXmlElement* isotope_E = isotope->ToElement();
//        double mass = atof( isotope_E->Attribute("mz") );
//        double intensity = atof( isotope_E->Attribute("Intensity") );        
//        
//        // store:
//        masses.push_back( mass );
//        intens.push_back(intensity );
//        totIntens += intensity;
//      }
//      
//      // normalize the intensities:
//      vector<double>::iterator I = intens.begin();
//      while( I != intens.end()){
//        *I /= totIntens;
//        I++;
//      }
//      
//      
//      PeptideIsotopeDisribution* newIsoPep = new PeptideIsotopeDisribution( masses, intens, chargeState, pinName, pinSQ, id , rtSegment);
//      newIsoPep->show_info();
//      
//      // store the profile:
//      allExternalPepdistributions.insert( make_pair(  masses[0], *newIsoPep ) );
//      delete newIsoPep;
//      newIsoPep = NULL;
//    }
//    
//  }
//  
//}




//////////////////////////////////////////////////
// init the retention time segments
void ExternalIsotopicDistribution::initRetentionTimeSegments( double start, double end){
  
  if( !allExternalPepdistributions.empty() ){
    
    double maxS = 0;
    multimap<double, PeptideIsotopeDisribution>::iterator P = allExternalPepdistributions.begin(  );
    while( P != allExternalPepdistributions.end() ){
      
      if( P->second.getRTSegment() > maxS ){
        maxS = P->second.getRTSegment(); 
      }
      
      P++;
    }
    
    // divide into segments:
    double Ssize = (end - start) / maxS;
    
    P = allExternalPepdistributions.begin(  );
    while( P != allExternalPepdistributions.end() ){
      
      double s = P->second.getRTSegment();
      
      double lStart = (s - 1.0) * Ssize;
      P->second.setRTStart( lStart );
      double lEnd = (s) * Ssize;
      P->second.setRTEnd( lEnd );
      
      P++;
    }
  }
  
  
}
