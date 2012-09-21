///////////////////////////////////////////////////////////////////////////
//
//  written by Lukas N Mueller, 30.3.05
//  Lukas.Mueller@imsb.biol.ethz.ch
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
//
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
// 


#ifndef _USE_simple_math2
#define _USE_simple_math2

class simple_math2 {
  
private:
  
  bool HIGH_CHECK;
  bool LOW_CHECK;

public:
  
  simple_math2();
  
  static const double T_TEST_001[12];
  static const double T_TEST_002[12];
  static const double T_TEST_01[12];
  static const double T_TEST_02[12];
  static const double T_TEST_05[12];
  static std::string ALPHA_VALUE;
  
  // this structure provides the function to compare
  // in the sorting algorithm:
  struct VECTOR_OPERATOR{
    // provide the compare function for sort:
    bool operator()(const std::pair<double,  void*> A,const std::pair<double,  void*> B) const{
      // check if they have same mass
      if(A.first == B.first){
        return false;
      }
      else{
        return A.first > B.first;
      }
    }
  };
  
  // if they fall into the m/z tolerance window
  static bool compareMassValuesAtPPMLevel( double, double , double );
  // get the masse error at the PPM value     
  static double getMassErrorAtPPMLevel( double , double );
  
  void ITERATIVE_OUTLIER_DETECTION_BY_DIXON(std::vector<double>* IN);
  
  void ITERATIVE_OUTLIER_DETECTION_BY_DIXON(std::vector< std::pair<double, void*> >* IN);
  
  void OUTLIER_DETECTION_BY_DIXON(std::vector<double>* IN);
  
  void OUTLIER_DETECTION_BY_DIXON(std::vector< std::pair<double, double> >* IN);
  
  void OUTLIER_DETECTION_BY_DIXON(std::vector< std::pair<double, void*> >* IN);
  
  bool check_T_TEST( double IN , int SAMPLE_NB);
    
};


#endif
