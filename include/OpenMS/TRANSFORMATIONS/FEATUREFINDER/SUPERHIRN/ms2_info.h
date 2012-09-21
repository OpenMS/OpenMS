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
// CLASS ms2_info:
//
// variable description:
//
//
// function description:
//
//
// **********************************************************************//


#ifndef MS2_INFO_H
#define MS2_INFO_H

#include <string>
#include <vector>
#include <map>

using namespace std;

class ms2_info{

    
  ////////////////////////////////////////////////
  // declaration of the private members:
  
private:
  
  int ID;
  
  string SQ;
  string FULL_SQ;
  string PREV_AA;
  vector<string> AC;
  string ORIGINAL_INTERACT_FILE;
  string MS2_TYPE_TAG;
  
  // peptide prophet analysis:
  float PEP_PROB;
  
  // sorcerer search results:
  double XCORR;
  double DELTA_CN;
  
  double MONO_MZ;
  double THEO_MZ;
  double NEUTRAL_MR;
  
  int CHRG;
  int SCAN_START;
  int SCAN_END;
  
  double TR;
    
  map< int, double> MOD_LIST;  
  
  // static values:
  static const double _MONO_H;
  static const double _MONO_O;

  ////////////////////////////////////////////////
  // declaration of the public members:
  
  
public:
    static const double mono_mass[26];
  static const char AA[20];   
  static double MS2_TR_TOL;
  static bool THEO_MATCH_MODUS;
  static double MS2_MZ_PPM_TOLERANCE;
  
  // class destructor
  ~ms2_info();
  
  // class constructor
  ms2_info();
  ms2_info(int);
  ms2_info(string, string, float);
  ms2_info(string, string, int, float);
  ms2_info(string, string, float, int, int);
  // class copy constructor
  ms2_info(const ms2_info&);
  // class copy constructor
  ms2_info(const ms2_info*);
  
  
  //////////////////////////////////////////////////
  // overload operators:
  ms2_info& operator=(const ms2_info&);
  bool operator==(const ms2_info&);
  ms2_info& operator<=(const ms2_info&);
  ms2_info& operator>=(const ms2_info&);
  ms2_info& operator<(const ms2_info&);
  ms2_info& operator>(const ms2_info&);
  
  // add modification
  void add_modification(int, double);
  map< int, double>::iterator get_Modification_list_start(){return MOD_LIST.begin();};
  map< int, double>::iterator get_Modification_list_end(){return MOD_LIST.end();};
  map< int, double>::iterator find_Modification(int pos){return MOD_LIST.find(pos);};
  map< int, double>* get_Modification_list(){return &(MOD_LIST);};
  bool check_MODIFICATION(){ return !MOD_LIST.empty();};

  // calculates the theoretical mass from a sequence:
  void set_THEO_MASS_from_SQ();
  double get_THEO_MZ(){return THEO_MZ;};
  // sets modificatied SQ:
  void set_FULL_SQ();
  void set_SQ(string IN){SQ = IN; set_THEO_MASS_from_SQ();set_FULL_SQ();};
  
  // show info:
  void show_info();
  
  // check whethere proteotypic peptide:
  bool get_PROTEO_TYPE();
  // check the tryptic state:
  // 2: full tryptic
  // 1: semi tryptic
  // 0: non tryptic
  int get_TRYPTIC_STATE();
    
  
  // AC functions:
  // check if this AC or not:
  bool compare_AC( string );
  // search a pattern in the  AC list:
  bool search_AC_pattern( string );


  ///////////////////////////////
  // start here all the get / set
  // function to access the
  // variables of the class
  string get_SQ(){return SQ;};
  string get_MOD_SQ(){return FULL_SQ;};
  string get_TOTAL_SQ(){return (get_PREV_AA( )+ "." + get_MOD_SQ());};
  string get_AC(){return *(AC.begin());};
  vector<string> get_ALL_AC(){return AC;};
  vector<string>::iterator get_ALL_AC_START(){return AC.begin();};
  vector<string>::iterator get_ALL_AC_END(){return AC.end();};
  bool find_AC( string );
  void set_AC(string);
  float get_PEP_PROB(){return PEP_PROB;};
  void set_PEP_PROB(float IN){PEP_PROB = IN;};
  
  double get_MONO_MZ(){return MONO_MZ;};
  void set_MONO_MZ(double);

  double get_NEUTRAL_MR(){return NEUTRAL_MR;};
  void set_NEUTRAL_MR(double);

  int get_CHRG(){return CHRG;};
  void set_CHRG(int IN){CHRG = IN;};

  int get_SCAN(){return SCAN_START;};
  int get_SCAN_START(){return SCAN_START;};
  void set_SCAN_START(int IN){SCAN_START = IN;};

  int get_SCAN_END(){return SCAN_END;};
  void set_SCAN_END(int IN){SCAN_END = IN;};

  int get_ID(){return ID;};
  
  double get_DELTA_CN(){ return DELTA_CN;};
  void set_DELTA_CN( double IN){ DELTA_CN = IN;};

  double get_XCORR(){ return XCORR;};
  void set_XCORR( double IN){ XCORR = IN;};
  
  void set_MS2_TYPE_TAG(string IN){MS2_TYPE_TAG = IN;};
  string get_MS2_TYPE_TAG(){return MS2_TYPE_TAG;};

  string get_ORIGINAL_INTERACT_FILE(){ return ORIGINAL_INTERACT_FILE;};
  void set_ORIGINAL_INTERACT_FILE( string IN){ ORIGINAL_INTERACT_FILE = IN;};
  
  // access the retentino parameter:
  double getRetentionTime(){return TR;};
  void setRetentionTime(double IN){TR = IN;};

  double get_MONO_AA_MASS(int);

  string get_PREV_AA( ){return PREV_AA;};
  void set_PREV_AA( string IN ){PREV_AA = IN;};
};

#endif

    
