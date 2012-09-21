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


#ifndef _LCMSCData_h_
#define _LCMSCData_h_

#include <map>

using namespace std;

typedef map<int, LC_elution_peak> elution_peak_list;
typedef map<int, LC_elution_peak>::iterator elution_peak_list_ITERATOR;
typedef map< double, elution_peak_list> MZ_LIST;
typedef MZ_LIST::iterator MZ_LIST_ITERATOR;


class LCMSCData
{
public:
  
		LCMSCData();
		LCMSCData(const LCMSCData&);
		LCMSCData(const LCMSCData*);  
		virtual ~LCMSCData();

  //////////////////////////////////////////////////
  // overload operators:
  LCMSCData& operator=(const LCMSCData&);
  
  // prints profile of the extracted LC elution peaks:
  void print_profile_of_LC_peaks(ofstream*);
 
  
  ///////////////////////////////
  // start here all the get / set
  // function to access the
  // variables of the class
  
  // start / end of the mz list:
  MZ_LIST_ITERATOR get_DATA_start(){return DATA.begin();};
  MZ_LIST_ITERATOR get_DATA_end(){return DATA.end();};
  
  // find data of a specific m/z:
  MZ_LIST_ITERATOR get_MZ_by_iterator(double);
  // add data into the structure:
  void add_LC_elution_peak(double , LC_elution_peak*);
  
  // get a list of m/z observed in a scan +/- scan_tolerance:
  // return the area of the LC elution peaks 
  vector<LC_elution_peak> get_MZ_list(int);
  vector<LC_elution_peak> get_MZ_list(int, int);
  
  // get all extracted LC peaks:
  vector<LC_elution_peak*> get_ALL_peak();
  // get a vector with all LC peaks ordered by their score:
  vector<LC_elution_peak*> get_ALL_peak_ordered();
  
  
  /*
  //////////////////////////////////////////////
  // operator which compares LC elution peaks 
  // by their score:
  struct OPERATOR_LC_PEAK{
    bool operator()(LC_elution_peak* A,LC_elution_peak* B) const{
      return A->get_LC_score() < B->get_LC_score();
    }
  };
*/
  
  // minimal MS peak intensity allowed:
  static float intensity_min_threshold;
  

 
protected:
    
    // main data struture:
    MZ_LIST DATA;
  
  
};
#endif
