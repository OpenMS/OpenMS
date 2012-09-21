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


#ifndef PROCESS_DATA_H
#define PROCESS_DATA_H

#include <map>
#include <vector>
#include <string>
#include <list>

namespace OpenMS
{

typedef std::multimap<int, ms_peak > elution_peak;
typedef std::vector< elution_peak > MZ_series;
typedef std::vector< elution_peak >::iterator MZ_series_ITERATOR;
typedef std::vector< elution_peak >::reverse_iterator MZ_series_REV_ITERATOR;
typedef std::multimap< double, MZ_series> MAIN_DATA_STRUCTURE;
typedef MAIN_DATA_STRUCTURE::iterator MAIN_ITERATOR;


class OPENMS_DLLAPI Process_Data{

    
    ////////////////////////////////////////////////
    // declaration of the private members:

private:
  
  
protected:
  
  ////////////////////////////////////////////////
  // declaration of the public members:
  // map that tracts all observed masses in keys
 
  
  // clustering of ms peak along time axis:
  bool TIME_CLUSTERING_BY_RETENTION_TIME;
  
  // max_distance from next elution peak member in scan numbers:
  int max_inter_scan_distance;
  
  // cluster data structure:
  LCMSCData* DATA;
  
  // and stores as contents vectors of type 
  MAIN_DATA_STRUCTURE pMZ_LIST;
  
  // tracks the number of observed mz cluster
  // elements:
  std::map<double, int> MZ_CLUSTER;
  unsigned int LC_elution_peak_counter;
  
  
  ///////////////////////////////////////////////
  // data processing classes:
  BackgroundControl* backgroundController;
  

  
public:
  
  // minimal intensity level:
  static float INTENSITY_THRESHOLD;
  
  // m/z tolerance value:
  static double MZ_TOLERANCE;
  
  // max_distance from next elution peak member in min.:
  static double max_inter_scan_retention_time_distance;
  
  // define minimal number of members in LC elution peaks cluster
  static int min_nb_cluster_members;
    
    static std::map<int, float> scan_TR_index;
  
  // to track detected monoistopic mass for debugging:
  static bool MonoIsoDebugging;
  static double DebugMonoIsoMassMin; 
  static double DebugMonoIsoMassMax; 
  static double MS1_intensity_apex_percentil_cutoff;
  static double MS1_TR_RESOLUTION;
  // if data are in centroid form or not:
  static bool CENTROID_DATA_MODUS;

  // class destructor
  ~Process_Data();
  
  // class constructor
  Process_Data( );
  // class copy constructor
  Process_Data(const Process_Data&);
  
  
  ///////////////////////////////////////////////////////////////////////////////
  // inputs raw /centroided  data into the object:
  void add_scan_raw_data(int , float , CentroidData* );    
  // inputs raw data into the object:
  void add_scan_raw_data(std::vector<ms_peak>);
  
  //////////////////////////////////////////////////
  // overload operators:
  Process_Data& operator=(const Process_Data&);
  Process_Data& operator<=(const Process_Data&);
  Process_Data& operator>=(const Process_Data&);
  Process_Data& operator<(const Process_Data&);
  Process_Data& operator>(const Process_Data&);

  // insert an already observed mz into the data structure, checks
  // if it belongs to an existing LC elution peak or starts a new one:
  void insert_observed_mz(MAIN_ITERATOR, ms_peak*);
  // insert a newly observed mz into the data structure
  void insert_new_observed_mz( ms_peak* );
  
  
  // converts DeconvPeak list to ms_peak vector
  void convert_ms_peaks(int ,double , std::list<DeconvPeak>& ,std::vector<ms_peak>& );

  // check if the ms peak is in the selected mz, z, int range
  bool filterDeisotopicMSPeak(ms_peak* );
  
  //////////////////////////////////////////////////////////////////
  // function which check if a data structure iterator is similar
  // to a peak and should be considered 
  // returns 1 if ok
  // returns 0 if not
  // returns -1 if scan range exceeded
  int compareIteratorToPeak(ms_peak* , MAIN_ITERATOR );
  // checks if a mz value has already been seen,
  // also look for very close ones and cluster them
  MAIN_ITERATOR check_MZ_occurence(ms_peak*);
    
  
  ///////////////////////////////////////////////////////////////////////////////
  // process a series of MS peaks
  // set the signal to noise level:
  void processMSPeaks(std::multimap<int, ms_peak>* );
    
  
  ///////////////////////////////////////////////////////////////////////////////
  // get the  full summed up intensity
  double getPeakIntensitySum(double );
  
  
  // check if a peak with this scan number belong to this elution cluster:
  bool check_elution_peak_belong(MZ_series_ITERATOR, ms_peak*);
  // returns the distance to this elution peak:
  int getElutionPeakDistance(MZ_series_ITERATOR, int);

  // runs through the whole data structure and puts the elution_peaks into
  // a proper LC_elution peak object
  void extract_elution_peaks();
  
  // check if this elution peak is accepted as a really LC-elution peak:
  bool check_elution_peak(MZ_series_ITERATOR);
  
  // convert the MZ_series elution peak element into a LC_elution_peak object
  void convert_to_LC_elution_peak(MZ_series_ITERATOR, double);
  
  // find a retention time by the scan number:
  double find_retention_time(double);
  
  // find closest match mz mass in the main structure
  MAIN_ITERATOR find_closest_mz_match(double);
  
  // go back to the MS1 level and
  // find the correct precursor mass by mz and z:
  void adjustCorrectToMS1Precursor( double* , int, int, int);
    
  
  ///////////////////////////////////////////////////////////
  // access methods to the object variables:
  
  // get an observed MZ mass, otherwise end of list iterator
  MAIN_ITERATOR get_MZ(double);
  // get an observed MZ mass, otherwise end of list iterator
  MAIN_ITERATOR get_MZ_lower_bound(double);
  // get end of MZ list:
  MAIN_ITERATOR get_MZ_LIST_end();
  // get start of MZ list:
  MAIN_ITERATOR get_MZ_LIST_start();
  // erase element in MZ list:
  void erase_MZ_LIST_element(MAIN_ITERATOR); 
  int getNbMSTraces(){ return pMZ_LIST.size();};
  
  
  double getMinimalIntensityLevel(){return INTENSITY_THRESHOLD;};
  


  // access the MZ_CLUSTER:
  // find element numbers:
  std::map<double, int>::iterator get_nb_MZ_cluster_elements(double);
  // erase an element:
  void erase_MZ_cluster_element(std::map<double, int>::iterator);
  // insert an element:
  void insert_MZ_cluster_element(double, int);
  
  // add the scan vs TR index to the data structure:
  void add_scan_TR_index(std::map<int, float> IN){scan_TR_index = IN;};

  // get the processed data:
  LCMSCData* get_processed_data(){return DATA;};
    
  // increase the LC_elution_profile counter:
  void increase_LC_elution_peak_counter(){LC_elution_peak_counter++;};
  unsigned int get_LC_elution_peak_counter(){return LC_elution_peak_counter;};
  
  // get the maximal scan distance between two same monoisotopic masses
  int getMaxScanDistance(){ return max_inter_scan_distance;};
  void setMaxScanDistance(int in){ max_inter_scan_distance = in;};
  
  
  // build up an index scan vs retention time:
  static void insert_into_scan_TR_index(int IN, float TR){scan_TR_index.insert(std::pair<int, float>(IN,TR));};

};

} // ns

#endif

    
