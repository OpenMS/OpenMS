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
// CLASS FT_PeakDetectController:
//
// variable description:
//
//
// function description:
//
//
// **********************************************************************//


#ifndef USE_FT_PEAK_DETECT_CONTROLLER_H
#define USE_FT_PEAK_DETECT_CONTROLLER_H

namespace OpenMS
{

class OPENMS_DLLAPI FT_PeakDetectController{

    
    ////////////////////////////////////////////////
    // declaration of the private members:

private:

  ////////////////////////////////////////////////
  // declaration of the public members:
  
  // LCMS runs
  //LC_MS* THIS_LC_MS;  
  std::vector<feature> fakeFeatureList;
  std::vector<LC_MS> LC_MS_RUNS;  
  
  // paths:
  std::string targetMzXML;
  std::string SOURCE_DIR;
  std::string OUTPUT_DIR;
  
  
public:
  
  LC_MS* THIS_LC_MS;
  
  typedef std::map<double, RawData*> Map;
  typedef std::vector<Map> Vec;

    static bool CREATE_FEATURE_ELUTION_PROFILES;
  static bool LCelutionPeakDebugging;
  static double LCelutionPeakMassMin;
  static double LCelutionPeakMassMax;
  
  static MS2_feature* SearchedM2Feature;

  static bool FEATURE_FAKE_INSERTION_BASED_ON_MS2_FEATURE;


  // class destructor
  ~FT_PeakDetectController();
  
  // class constructor
  FT_PeakDetectController();
  // class copy constructor
  FT_PeakDetectController(const FT_PeakDetectController&);
  
  
  /////////////////////////////////////////////////////////
  // function for batch processing of mzXML data
  // parses LC-MS from runs from a directory or file of raw mzXML data:
  void parseMzXMLData();
  
  
  //////////////////////////////////////////////////
  // mzXML parsing functions for a single MzXML file:
  // start the scan parsing of a mzXML file:
  void start_scan_parsing_of_mzXML_file(Vec datavec);
  
  
  // **** for the MS1 level post processing:
  // process MS1 level data
  void process_MS1_level_data_structure( FT_PEAK_DETEC_mzXML_reader*);
  // adds an elution peak to the LC/MS run:
  void add_raw_peak_to_LC_MS_run( LC_elution_peak* );
  // function to add the elution profile to the feature:
  void addLCelutionProfile( feature* , LC_elution_peak* );
  
  /////////////////////////////////////////////////////////////
  // reads already paths of existing LC-MS runs in xml format into the
  // memory
  // for now, open file system for every check, but otherwise could eb done
  // in the constructor
  bool checkIfFeatureExtractionExists( std::string );
    

  // **** for the MS2 level post processing:
  // process MS2 level data
  void process_MS2_level_data_structure( FT_PEAK_DETEC_mzXML_reader*);
  // processes the tracted signals on teh MS2 level
  void extract_MS2_elution_features();
  // combine the MS2 feature trace data to the MS1 features:
  void associateMS2FeatureToMS1Feature( MS2_feature*);
  // add an observed MS2 feature to the MS1 feature
  // if an observation is already there, then
  // construct a merged MS2 feature
  void addMS2FeatureToMS1Feature( MS2_feature*, feature* );
    
  // construct here fake ms1 features based on a observed MS2 feature
  // which however could not be matched to a exiting ms1 feature
  void constructMS1FeatureFromMS2Feature( MS2_feature* );
    

  
  /////////////////////////////////////////////////////////
  // write a parsed LC/MS into directory:
  void write_out_parsed_LC_MS(LC_MS*);
  // add fake MS/MS information for the MS1 feature:
  void addFakeMSMSToFeature( feature* );

  
  
  
  // this structure provides the function to find
  // the correct MS1 fetaure for a target MS2 feature:
  class MS2ToMS1Comparer{
    
public:
    
    
    // with unary predicator to provide the compare function for sort:
    bool operator()( const feature MS1 ) const{

      
      /*
      // use here the m/z tolerance of the MS1 since precursors are compared
      double deltaMZ = fabs( MS1.MONO_MZ - FT_PeakDetectController::SearchedM2Feature->getPrecursorMZ() );
      if( deltaMZ > feature::MZ_TOL ){
        return false;
      }
       */

      // use here the m/z tolerance of the MS1 since precursors are compared
      if( !feature::compareFeatureMassValuesAtPPMLevel( MS1.MONO_MZ, FT_PeakDetectController::SearchedM2Feature->getPrecursorMZ()) ){
        return false;
      }
      
      // charge state:
      if( MS1.charge_state != FT_PeakDetectController::SearchedM2Feature->getPrecursorChrg() ){
        return false;
      }
      
      // compare teh retention time borders:
      double deltaTR = MS1.TR_START - FT_PeakDetectController::SearchedM2Feature->getStartTR();
      if( deltaTR > feature::TR_TOL ){
        return false;
      }
      
      deltaTR = FT_PeakDetectController::SearchedM2Feature->getEndTR() - MS1.TR_END;
      if( deltaTR > feature::TR_TOL ){
        return false;
      }
      
      
      return true;
    }
    
  };
  
  
  
  
  
  //////////////////////////////////////////////////
  // overload operators:
  FT_PeakDetectController& operator=(const FT_PeakDetectController&);
  FT_PeakDetectController& operator<=(const FT_PeakDetectController&);
  FT_PeakDetectController& operator>=(const FT_PeakDetectController&);
  FT_PeakDetectController& operator<(const FT_PeakDetectController&);
  FT_PeakDetectController& operator>(const FT_PeakDetectController&);
  
  
  ///////////////////////////////
  // start here all the get / set
  // function to access the
  // variables of the class
  
  // target file:
  void set_target_file(std::string IN){targetMzXML = IN;};
  std::string get_target_file(){return targetMzXML;};
    
  // get the vector of LC/MS runs:
  std::vector<LC_MS> get_parsed_DATA(){return LC_MS_RUNS;};
  bool getParsedDataEmpty(){return LC_MS_RUNS.empty();};
  std::vector<LC_MS>::iterator get_parsed_DATA_START(){return LC_MS_RUNS.begin();};
  std::vector<LC_MS>::iterator get_parsed_DATA_END(){return LC_MS_RUNS.end();};
};

} // ns

#endif

    
