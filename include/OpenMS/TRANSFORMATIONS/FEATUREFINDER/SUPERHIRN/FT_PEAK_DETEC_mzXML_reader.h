///////////////////////////////////////////////////////////////////////////
//
//  PEAK DETECTION OF FOURIER TRANSFORME MS INSTRUMENT DATA
//
//  written by Markus Mueller, markus.mueller@imsb.biol.ethz.ch
//  and Lukas Mueller, Lukas.Mueller@imsb.biol.ethz.ch
//  October 2005
//
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
//
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
// 
//


#ifndef FT_PEAK_DETEC_MZ_XML_READER_H
#define FT_PEAK_DETEC_MZ_XML_READER_H

// **********************************************************************//
// CLASS file reader:
// provides function to open /reade / modify and close text files
//
// variable description:
//
//
// function description:
//
//
// **********************************************************************//


// file structure remapping
// for the mzXML parser ramp:
// - for old ramp use FILE*
// - for new cramp use RAMPFILE*
typedef FILE* MZXML_FILE;
//typedef RAMPFILE RAMP_FILE;



class FT_PEAK_DETEC_mzXML_reader{

    
    ////////////////////////////////////////////////
    // declaration of the private members:
    
private:
  
  Process_Data* MS1_LC_MS_DATA_PROCESSOR;
  // MS2_Process_Data* MS2_LC_MS_DATA_PROCESSOR;

  string current_file;
  //MZXML_FILE file_handler;
  //RAMP_FILE* ramp_file_Struct;
  
  off_t index_offset;
  off_t* scan_index;

  // mzXML parameters
  int total_scan;  
  map<int,float> scan_TR_index;
  double minRT;
  double maxRT;
  
  
  
  int SCAN_MIN;
  int SCAN_MAX;
    
  int nbMS2Scans;
  
  ////////////////////////////////////////////////
  // declaration of the public members:
  
public:
  
  
  typedef map<double, RawData*> Map;
  typedef vector<Map> Vec;
    
    static int	sfReportMonoPeaks; // 1 if info about monoisotopic peaks should be written to mono_peaks.txt
  static string sfDebugDirectory; // Directory where peak detection debug files are written
  static int	sfReportScanNumber; // if sfReportMonoPeaks is set to 1, details about this spectrum will be written to debug files 

  static int MS1_base_inter_scan_distance;
  static int MS2_base_inter_scan_distance;
  static bool MS2_PEAK_PROCESSING;
  static double TR_MIN;
  static double TR_MAX;

  
  static vector<double> PEAK_EXTRACTION_SCAN_LEVELS;
  static vector<double> FRAGMENT_MASS_SCAN_LEVELS;
    
  // class destructor
  ~FT_PEAK_DETEC_mzXML_reader();
  // class constructor
  FT_PEAK_DETEC_mzXML_reader();
  
  //////////////////////////////////////////////////
  // file open / read / close function
  // set the current mzXML file;
  //void open_mzxml_file();
  // close the current mzXML file;
  //void close_mzxml_file();
  // set indexes of the current mzXML file;
  void set_current_indexes(double pminrt, double pmaxrt);
  
  // reads the ms data from a mzXML file opened by teh handler
  void read_mzXML_DATA(Vec datavec);

  // get a MS scan at a given scan number within
  // a mass range and stores it into the ms region structure
  void get_MS_scan(off_t IN, double TR, RawData* data);

  
  ///////////////////////////////////////////////////////
  // MS1 level processing routine functions:
  // process the MS1 level input data:
  // - construct a RawData object with input peaks
  // - centoid them
  // - add them to Process_Data Structure:
  void processMS1InputData(int , float, RawData* data);
 
  
  ///////////////////////////////////////////////////////
  // MS2 level processing routine functions:
  // process the MS2 level input data:
  // - construct a RawData object with input peaks
  // - centoid them
  // - add them to Process_Data Structure:
  //void processMS2InputData(int , double, float, int,  MY_RAMP_PEAKS* );
  
  

  // set the maximal inter-monoisotopic distance
  // for the same LC-elution peak
  int setInterMonoIsotopicLCDistance( int , int, int);
  //int setInterMonoIsotopicLCDistance( int , int, int, double);
  
  // get a previous scan at the MS level from the current scan
  //int getPreviousMSXScan( int, int );

  
  ///////////////////////////////////////////////////////////////////////////////
  // go back to the MS1 level and
  // find the correct precursor mass by mz and z:
  //void findMS1PrecursorData( double* , int , int, int);
    

  ///////////////////////////////////////////////////////////////////////////////
  // check if the scan number should be processed by MS Precursor Mass extraction
  bool checkMSPrecursorMassScan( int );
  // check if the scan number should be processed as MSn FragmentMass spectrum
  bool checkMSFragmentMassScan( int );

  
  // print the scan header:
  //void print_scan_header( ScanHeaderStruct* );  
  
  ///////////////////////////////
  // start here all the get / set
  // function to access the
  // variables of the class
  
  // set the current mzXML file handler;
  //void set_file_handler(MZXML_FILE IN){file_handler = IN;};
  //MZXML_FILE get_file_handler(){return file_handler;};
  //RAMP_FILE* get_Ramp_file_handler(){return ramp_file_Struct;};
  //MZXML_FILE get_Ramp_file_handler(){return file_handler;};

  // total scan access:
  int get_total_scan(){return total_scan;};
  off_t get_index_off_set(){return index_offset;};
  
  // get a scan with the index:
  off_t get_scan(int i){return scan_index[i];};
  
  
  //string get_current_file(){return current_file;};
  //void set_current_file(string IN){current_file = IN;};

  // get the MS1 processed data:
  Process_Data* get_processed_MS1_data_structure(){return MS1_LC_MS_DATA_PROCESSOR;};
  // get the MS2 processed data:
  //MS2_Process_Data* get_processed_MS2_data_structure(){return MS2_LC_MS_DATA_PROCESSOR;};

  // build up an index scan vs retention time:
  void insert_into_scan_TR_index(int IN, float TR){Process_Data::insert_into_scan_TR_index(IN, TR);};
  
  // converts DeconvPeak list to ms_peak vector
  //void convert_ms_peaks(int, double, list<DeconvPeak>&,vector<ms_peak>&);

  // delete existing debug files before writing in append mode
  //void delete_existing_debug_file();

};

#endif

