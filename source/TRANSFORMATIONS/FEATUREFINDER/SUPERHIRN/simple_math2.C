///////////////////////////////////////////////////////////////////////////
//
//  written by Lukas N Mueller, 30.3.05
//  Lukas.Mueller@imsb.biol.ethz.ch
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
//
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
// 

#include <math.h>
#include <vector>
#include <string>
#include <algorithm>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/simple_math2.h>

using namespace std;

// significance values:
const double simple_math2::T_TEST_001[12]={.999,.964,.895,.822,.763,.716,.675,.647,.544,.491,.455,.430};
const double simple_math2::T_TEST_002[12]={.998,.949,869,.792,.731,.682,.644,.614,.515,.464,.430,.407};
const double simple_math2::T_TEST_01[12]={.886,.679,.559,.484,.433,.398,.370,.349,.284,.251,.230,.216};
const double simple_math2::T_TEST_02[12]={.782,.561,.452,.387,.344,.314,.291,.274,.220,.193,.176,.165};
const double simple_math2::T_TEST_05[12]={.782,.561,.452,.387,.344,.314,.291,.274,.220,.193,.176,.165};
string simple_math2::ALPHA_VALUE="0.05";

simple_math2::simple_math2(){
  LOW_CHECK = true;
  HIGH_CHECK = true;
}

//////////////////////////////////////////////////////
// compare to masses at the PPM value and decided
// if they fall into the m/z tolerance window
bool simple_math2::compareMassValuesAtPPMLevel( double mzA, double mzB, double PPM_TOLERANCE ){
  
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

//////////////////////////////////////////////////////
// get the masse error at the PPM value 
double simple_math2::getMassErrorAtPPMLevel( double mz, double PPM_TOLERANCE ){
  double ppmValue = mz / 1000000.00;
  return ppmValue * PPM_TOLERANCE;
}
//////////////////////////////////////////////////
// iterative detection of outliers by the DIXON's test (Gibbons, 1994):
void simple_math2::ITERATIVE_OUTLIER_DETECTION_BY_DIXON(vector<double>* IN){
  
  vector < pair<double, double> > TMP;
  vector<double>::iterator P = IN->begin();
  while( P != IN->end() ){
    TMP.push_back( make_pair( *P, 0) );
    P++;
  }
  IN->clear();
  
  
  ////////////////
  // start here teh outlier detection iteratively:
  unsigned int size;
  do{
    size = TMP.size();
    // make a outlier detection:
    OUTLIER_DETECTION_BY_DIXON( &TMP );
  }
  while( size != TMP.size() );
  
  
  vector< pair<double, double> >::iterator T = TMP.begin();
  while( T != TMP.end() ){
    IN->push_back( (*T).first );
    T++;
  }   
  
}

//////////////////////////////////////////////////
// iterative detection of outliers by the DIXON's test (Gibbons, 1994):
void simple_math2::ITERATIVE_OUTLIER_DETECTION_BY_DIXON(vector< pair<double, void*> >* IN){
  
  ////////////////
  // start here teh outlier detection iteratively:
  unsigned int size;
  do{
    size = IN->size();
    // make a outlier detection:
    OUTLIER_DETECTION_BY_DIXON( IN );
  }
  while( size != IN->size() );
  
}

//////////////////////////////////////////////////
// detection of outliers by the DIXON's test (Gibbons, 1994):
void simple_math2::OUTLIER_DETECTION_BY_DIXON(vector<double>* IN){
  
  vector < pair<double, double> > TMP;
  vector<double>::iterator P = IN->begin();
  while( P != IN->end() ){
    TMP.push_back( make_pair( *P, 0) );
    P++;
  }
  IN->clear();
  
  // to the actual outlier detection:
  OUTLIER_DETECTION_BY_DIXON( &TMP );
  
  vector< pair<double, double> >::iterator T = TMP.begin();
  while( T != TMP.end() ){
    IN->push_back( (*T).first );
    T++;
  }  
}

//////////////////////////////////////////////////
// detection of outliers by the DIXON's test (Gibbons, 1994):
void simple_math2::OUTLIER_DETECTION_BY_DIXON(vector< pair<double, double> >* IN){
  
  int SAMPLE_SIZE = IN->size();
  
  if( (SAMPLE_SIZE > 2) && (SAMPLE_SIZE < 2000) ){
    
    vector< pair<double, double> >::iterator START = IN->begin();
    vector< pair<double, double> >::iterator END = IN->end();
    END--;
    
    double TAU_HIGH;
    double TAU_LOW;
    
    if( SAMPLE_SIZE < 8 ){
      
      double XN = (*END).first;
      END--;
      double XN_1 = (*END).first;
      double X1 = (*START).first;
      START++;
      double X2 = (*START).first;
      
      TAU_HIGH = ( XN - XN_1) / ( XN - X1 );
      TAU_LOW = ( X2 - X1) / ( XN - X1 );      
    }
    else if( SAMPLE_SIZE < 11 ){
      
      double XN = (*END).first;
      END--;
      double XN_1 = (*END).first;
      double X1 = (*START).first;
      START++;
      double X2 = (*START).first;
      
      TAU_HIGH = ( XN - XN_1) / ( XN - X2 ); 
      TAU_LOW = ( X2 - X1) / ( XN_1 - X1 ); 
      
    }
    else if( SAMPLE_SIZE < 14 ){
      
      double XN = (*END).first;
      END--;
      double XN_1 = (*END).first;
      END--;
      double XN_2 = (*END).first;
      double X1 = (*START).first;
      START++;
      double X2 = (*START).first;
      START++;
      double X3 = (*START).first;
      
      TAU_HIGH = ( XN - XN_2) / ( XN - X2 ); 
      TAU_LOW = ( X3 - X1) / ( XN_1 - X1 ); 
      
    }
    else{
      
      double XN = (*END).first;
      END--;
      END--;
      double XN_2 = (*END).first;
      double X1 = (*START).first;
      START++;
      START++;
      double X3 = (*START).first;
      
      TAU_HIGH = ( XN - XN_2) / ( XN - X3 ); 
      TAU_LOW = ( X3 - X1) / ( XN_2 - X1 ); 
      
    }
    
    
    // delete here outliers if detected:
    if( check_T_TEST( TAU_HIGH, SAMPLE_SIZE ) && HIGH_CHECK ){
      END = IN->end();
      END--;
      IN->erase( END );
    }
    
    if( check_T_TEST( TAU_LOW, SAMPLE_SIZE )  && LOW_CHECK ){
      START = IN->begin();
      IN->erase( START );
    }
    
  }
}

///////////////////////////////////////////////
// check T-TEST:
bool simple_math2::check_T_TEST( double IN , int SAMPLE_NB){
  
  // do only if more than 3 values;
  if( SAMPLE_NB < 3 ){
    return false;
  }
  
  const double* SIGN = NULL;
  if( ALPHA_VALUE == "0.001"){ 
    SIGN = T_TEST_001;
  }
  else if( ALPHA_VALUE == "0.002"){
    SIGN = T_TEST_002;
  }
  else if( ALPHA_VALUE == "0.01"){
    SIGN = T_TEST_01;
  }
  else if( ALPHA_VALUE == "0.02"){
    SIGN = T_TEST_02;
  }
  else if( ALPHA_VALUE == "0.05"){
    SIGN = T_TEST_05;
  }
  else{
    SIGN = T_TEST_01;
  }
  
  if( SAMPLE_NB < 11 ){
    if( SIGN[SAMPLE_NB - 3] <= IN ){
      return true;
    }
  }
  else if( SAMPLE_NB < 16 ){
    if( SIGN[8] <= IN ){
      return true;
    }   
  }
  else if( SAMPLE_NB < 21 ){
    if( SIGN[9] <= IN ){
      return true;
    }    
  }      
  else if( SAMPLE_NB < 25 ){
    if( SIGN[10] <= IN ){
      return true;
    }    
  }
  else if( SAMPLE_NB < 31 ){
    if( SIGN[10] <= IN ){
      return true;
    }    
  }
  else{
    if( SIGN[10] <= IN ){
      return true;
    }    
  }
  
  return false;
}


//////////////////////////////////////////////////
// detection of outliers by the DIXON's test (Gibbons, 1994):
void simple_math2::OUTLIER_DETECTION_BY_DIXON(vector< pair<double, void*> >* IN){
  
  // sort the vector:
  sort( IN->begin(), IN->end(), VECTOR_OPERATOR() );
  
  int SAMPLE_SIZE = IN->size();
  
  if( (SAMPLE_SIZE > 2) && (SAMPLE_SIZE < 2000) ){
    
    vector< pair<double, void*> >::iterator START = IN->begin();
    vector< pair<double, void*> >::iterator END = IN->end();
    END--;
    
    double TAU_HIGH;
    double TAU_LOW;
    
    if( SAMPLE_SIZE < 8 ){
      
      double XN = (*END).first;
      END--;
      double XN_1 = (*END).first;
      double X1 = (*START).first;
      START++;
      double X2 = (*START).first;
      
      TAU_HIGH = ( XN - XN_1) / ( XN - X1 );
      TAU_LOW = ( X2 - X1) / ( XN - X1 );      
    }
    else if( SAMPLE_SIZE < 11 ){
      
      double XN = (*END).first;
      END--;
      double XN_1 = (*END).first;
      double X1 = (*START).first;
      START++;
      double X2 = (*START).first;
      
      TAU_HIGH = ( XN - XN_1) / ( XN - X2 ); 
      TAU_LOW = ( X2 - X1) / ( XN_1 - X1 ); 
      
    }
    else if( SAMPLE_SIZE < 14 ){
      
      double XN = (*END).first;
      END--;
      double XN_1 = (*END).first;
      END--;
      double XN_2 = (*END).first;
      double X1 = (*START).first;
      START++;
      double X2 = (*START).first;
      START++;
      double X3 = (*START).first;
      
      TAU_HIGH = ( XN - XN_2) / ( XN - X2 ); 
      TAU_LOW = ( X3 - X1) / ( XN_1 - X1 ); 
      
    }
    else{
      
      double XN = (*END).first;
      END--;
      END--;
      double XN_2 = (*END).first;
      double X1 = (*START).first;
      START++;
      START++;
      double X3 = (*START).first;
      
      TAU_HIGH = ( XN - XN_2) / ( XN - X3 ); 
      TAU_LOW = ( X3 - X1) / ( XN_2 - X1 ); 
      
    }
    
    
    // delete here outliers if detected:
    if( check_T_TEST( TAU_HIGH, SAMPLE_SIZE ) ){
      END = IN->end();
      END--;
      IN->erase( END );
    }
    if( check_T_TEST( TAU_LOW, SAMPLE_SIZE ) ){
      START = IN->begin();
      IN->erase( START );
    }
    
  }
}