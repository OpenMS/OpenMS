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
// CLASS FEATURE:
// FEATURE of spectra containing: AC,
//
// variable description:
//
//
// function description:
//
//
// **********************************************************************//



#ifndef _FEATURE_H
#define _FEATURE_H

class feature{
  
  
  ////////////////////////////////////////////////
  // declaration of the private members:
  
private:
  
  /////////////////////////////////////////////
  // IDENTIFICATION PARAMETERS:
  // name of the spectra:
  map<double, vector< ms2_info> > MS2_SCANS;
  /////////////////////////////////////////////
  
  /////////////////////////////////////////////
  // raw MS peak paramaters:
  int scan_apex;
  int scan_start;
  int scan_end;
  double total_peak_area;
  double apex_peak_intensity;
  double PEAK_SCORE;
  double SignalToNoise;
  double BackgroundNoise;

  ////////////////////////////////////////////
  // Analysis parameters: 
  double alignment_error_up;
  double alignment_error_down;
  double SCORE_HOLDER;
  bool feature_match_status;
  double PI;
  
  ////////////////////////////////////////////
  // LC/MS run ID parameters:
  int spectrum_ID;
  int MASTER_ID;
  
  ///////////////////////////////////////////
  // string to store ms1 feature extra information:
  string featureExtraInformation;
  
  ///////////////////////////////////////////
  // LC elution profile:
  featureLCprofile* LCprofile;
  
  // static values:
  static double _MONO_H;
  static double _MONO_O;
  
  //////////////////////////////////////////
  // LC/MS matching things:
  map<int, feature> matched_feature_list;
  
  
  // ranges of m/z and tr:
  double TR_APEX;
  double MONO_MZ_START;
  double MONO_MZ_END;
  double MONO_MZ_ORIGINAL;
  
  //////////////////////////////////////////
  // associated MS2 feature:
  MS2_feature* MS2TraceFeature;
  
  ////////////////////////////////////////////////
  // declaration of the public members:
  
public:
    
    // tolerance in m/z and TR:
  static double PPM_MZ_TOL;
  //static double MZ_TOL;
  static double TR_TOL;
  static double PEPTIDE_PROBABILITY_THRESHOLD;
  static bool PRINT_ALL_ACs;
  static bool STORE_ALL_LOW_PROBABILITY_MS2_SCANS;
  
  
  
  double TR;   
  double MONO_MZ;
  double TR_START;   
  double TR_END;
  int charge_state;
  int feature_ID;

  // class destructor
  ~feature();
  
  // copy constructor:
  feature(const feature&);
  
  // class constructor
  feature(const feature*);
  // constructor for the object feature:
  feature(double, double, int, int, int,int, float, float, float);
  // constructor for the object feature:
  feature(float, int, int);
  feature( MS2_feature* );
  feature();
  // copy constructor:
  feature& operator=(const feature&);
  
  // show the content of the spectra
  void show_info();
  // show MS/MS spectra info:
  void showMS2consensSpectraInfo( );

  
  // writes out teh feature to a file:
  void print_2_file(ofstream*);
  
  //////////////////////////////////
  // comparision operators:
  bool operator==(const feature&);
  /*
  bool operator<(const feature&);
  bool operator>(const feature&);
  bool operator<=(const feature&);
  bool operator>=(const feature&);
  */
  
  // writes out the important information:
  void print_content(ofstream*,bool);
  // add MS/MS info to the feature:
  void add_MS2_info(ms2_info*);
  void add_MS2_info( map<double, vector<ms2_info> >*);
  bool get_MS2_info();
  bool get_MS2_info( double );
  bool check_MS2_empty(){return MS2_SCANS.empty();};
  void removeAllMS2Information(){return MS2_SCANS.clear();};
  int get_MS2_SCANS_SIZE(){return MS2_SCANS.size();};
  map<double, vector<ms2_info> >* get_MS2_SCAN_MAP(){return &MS2_SCANS;};
  map<double, vector<ms2_info> >::iterator get_MS2_SCANS_START(){return MS2_SCANS.begin();};
  map<double, vector<ms2_info> >::iterator get_MS2_SCANS_END(){return MS2_SCANS.end();};
  // get the best ms2 scan == closest to the apex:
  ms2_info* get_best_MS2_SCAN();
  ms2_info* get_best_MS2_SCAN( double );
  
  void setFeatureExtraInformation( string in){featureExtraInformation=in;};
  string getFeatureExtraInformation(){return featureExtraInformation;};
  
  
  // functions to set/access machted features:
  void add_matched_feature(feature*);
  map<int,feature>* get_match_list_REFERENCE(){return &matched_feature_list;};
  map<int,feature> get_match_list(){return matched_feature_list;};
  map<int,feature>::iterator get_match_list_start(){return matched_feature_list.begin();};
  map<int,feature>::iterator get_match_list_end(){return matched_feature_list.end();};
  map<int,feature>::iterator find_match_by_id(int ID){return matched_feature_list.find(ID);};

  // get feature at a certain LC-MS by LC_MS id
  feature* get_feature( int );
  
  // get the total peak are over all matched features:
  double get_MATCHED_peak_area();
  bool check_match_by_id(int);
  void erase_match_list(){matched_feature_list.clear();};
  // get the profile over all matched features:
  map<int, double> get_feature_profile();

  // return number of times this feature has been seen = nb_replicates in list plus 1!
  int get_replicate_match_nb(){return (matched_feature_list.size() + 1);};
  int get_matching_nb(){return get_replicate_match_nb();};
  // return the sum of all intensities over replicates:
  double get_replicate_intensity_sum();
  
  ///////////////////////////////
  // start here all the get / set
  // function to access the
  // variables of the class
  
  // access the parent mass of feature, calculated from the SQ
  double get_MZ(){return MONO_MZ;};
  void set_MZ(double in){MONO_MZ = in;};
  double get_MZ_START(){return MONO_MZ_START;};
  void set_MZ_START(double IN){MONO_MZ_START=IN;};
  double get_MZ_END(){return MONO_MZ_END;};
  void set_MZ_END(double IN){MONO_MZ_END=IN;};
  
  double get_THEO_MZ(){return get_best_MS2_SCAN()->get_MONO_MZ();};
  double get_THEO_MZ( double T){return get_best_MS2_SCAN( T )->get_MONO_MZ();};
  string get_AC(){return get_best_MS2_SCAN()->get_AC();};
  string get_AC(double T){return get_best_MS2_SCAN(T)->get_AC();};
  bool check_AC( string IN){return get_best_MS2_SCAN()->compare_AC(IN);};
  bool check_AC( string IN, double T){return get_best_MS2_SCAN( T )->compare_AC(IN);};
  string get_SQ(){return get_best_MS2_SCAN()->get_SQ();};
  string get_SQ(double T){return get_best_MS2_SCAN(T)->get_SQ();};
  string get_TOTAL_SQ(){return get_best_MS2_SCAN()->get_TOTAL_SQ();};
  string get_TOTAL_SQ(double T){return get_best_MS2_SCAN(T)->get_TOTAL_SQ();};
  string get_MOD_SQ(){return get_best_MS2_SCAN()->get_MOD_SQ();};
  string get_MOD_SQ(double T){return get_best_MS2_SCAN(T)->get_MOD_SQ();};
  double get_pep_prob(){return get_best_MS2_SCAN()->get_PEP_PROB();};
  double get_pep_prob(double T){return get_best_MS2_SCAN(T)->get_PEP_PROB();};
  string get_MS2_TYPE_TAG(){return get_best_MS2_SCAN()->get_MS2_TYPE_TAG();};
  string get_MS2_TYPE_TAG(double T){return get_best_MS2_SCAN(T)->get_MS2_TYPE_TAG();};
  int get_MS2_scan(){return get_best_MS2_SCAN()->get_SCAN_START();};
  int get_MS2_scan(double T){return get_best_MS2_SCAN(T)->get_SCAN_START();};
  map<double, vector<ms2_info> >* get_MS2_SCAN_LIST(){return &(MS2_SCANS);};
  map<double, vector<ms2_info> >::iterator get_MS2_SCAN_LIST_START(){return MS2_SCANS.begin();};
  map<double, vector<ms2_info> >::iterator get_MS2_SCAN_LIST_END(){return MS2_SCANS.end();};
  
  int get_scan_number(){return scan_apex;};
  void set_scan_number(int IN){scan_apex = IN;};
  int get_scan_start(){return scan_start;};
  void set_scan_start(int IN){scan_start = IN;};
  int get_scan_end(){return scan_end;};
  void set_scan_end(int IN){scan_end = IN;};
  int get_charge_state(){return charge_state;};
  void set_charge_state(int IN){charge_state = IN;};
  void set_peak_area(float IN){total_peak_area = IN;};
  float get_peak_area(){return total_peak_area;};
  // get peak area at a certain LC/MS:
  double get_peak_area( int );
  float get_apex_peak_intensity(){return apex_peak_intensity;};
  void set_apex_peak_intensity(double in){apex_peak_intensity=in;};
  void normalize_peak_area_by_factor(double factor){ total_peak_area *= factor;};
  
  double get_alignment_error_up(){return alignment_error_up;};
  void set_alignment_error_up(double IN){alignment_error_up = IN;};
  double get_alignment_error_down(){return alignment_error_down;};
  void set_alignment_error_down(double IN){alignment_error_down = IN;};
  
  void set_SCORE_HOLDER(double IN){SCORE_HOLDER = IN;};
  double get_SCORE_HOLDER(){ return SCORE_HOLDER;};
  
  double get_retention_time(){return TR;};
  void set_retention_time(double IN){TR = IN;};
  double get_retention_time_START(){return TR_START;};
  void set_retention_time_START(double IN){TR_START=IN;};
  double get_retention_time_END(){return TR_END;};
  void set_retention_time_END(double IN){TR_END=IN;};

  // original mz and Tr coordinates
  double get_raw_retention_time_apex(){return TR_APEX;};
  void set_raw_retention_time_apex(double IN){TR_APEX=IN;};
  double get_raw_MZ(){return MONO_MZ_ORIGINAL;};
  void set_raw_MZ(double IN){MONO_MZ_ORIGINAL = IN;};
  
  
  // feature ID:
  void set_feature_ID(int IN){feature_ID = IN;};
  int get_feature_ID(){return feature_ID;};
  
  void set_spectrum_ID(int IN){
    spectrum_ID = IN;
    /*
    if( MS2TraceFeature != NULL){
      MS2TraceFeature.set
    }
     */
  };
  int get_spectrum_ID(){return spectrum_ID;};

  void set_MASTER_ID(int IN){MASTER_ID = IN;};
  int get_MASTER_ID(){return MASTER_ID;};
  
  // check how many matches
  int get_nb_common_match();
  
  // get/set the peak score
  double get_peak_score(){return PEAK_SCORE;};
  void set_peak_score(double in){PEAK_SCORE = in;};
  
  // get the molecular mass of the corrsponding peptide!
  double get_Molecular_Mass();
  
  // fetaure PI:
  double get_FEATURE_PI(){return PI;};
  void set_FEATURE_PI(double IN){PI = IN;};
  
  // check charge states, in cases where a feature was
  // created based on a MS2 trace, charge state is unknown ( = -1 )
  // -> derivce the charge state from the matched feature (if this is 
  // also not -1 
  void deriveChargeStates( feature*  );
    
  
  // LC elution profile
  void setLCelutionProfile( featureLCprofile* IN ){LCprofile = IN;};
  featureLCprofile* getLCelutionProfile( ){return LCprofile;};
  	
  //////////////////////////////////////////////
  // parameters computed over matched features:
  double get_profile_retention_time();
  double get_profile_Molecular_Mass();

  /////////////////////////////////////////////
  // status if feature has been matched:
  bool get_feature_match_status(){ return feature_match_status; };
  void set_feature_match_status( bool IN){ feature_match_status = IN; };
  
  
  ///////////////////////////////////////////
  // access the MS2 feature
  void addMS2Feature( MS2_feature* in){ MS2TraceFeature = new MS2_feature( in );};
  void removeMS2Feature( ){ delete MS2TraceFeature; MS2TraceFeature=NULL;};
  MS2_feature* getMS2Feature(){ return MS2TraceFeature;};  
  
  double getSignalToNoise(){return SignalToNoise;};
  void setSignalToNoise(double in){SignalToNoise = in;};

  double getBackgroundNoiseLevel(){return BackgroundNoise;};
  void setBackgroundNoiseLevel(double in){BackgroundNoise = in;};

  
  //////////////////////////////////////////////
  // get static members:
  //static double get_MZ_TOL(){return MZ_TOL;};
  static double get_TR_TOL(){return TR_TOL;};
  static double get_MONO_H(){return _MONO_H;};
  
  // compare to masses at the PPM value and decided
  // if they fall into the m/z tolerance window
  static bool compareFeatureMassValuesAtPPMLevel( double , double );
    
  // get the masse error at the PPM value 
  static double getFeatureMassErrorAtPPMLevel( double);
  
};

#endif

