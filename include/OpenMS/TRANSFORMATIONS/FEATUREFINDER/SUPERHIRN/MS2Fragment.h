///////////////////////////////////////////////////////////////////////////
//
//  PEAK DETECTION OF FOURIER TRANSFORME MS INSTRUMENT DATA
//
//  written by Markus Mueller, markus.mueller@imsb.biol.ethz.ch
//  ( and Lukas Mueller, Lukas.Mueller@imsb.biol.ethz.ch)
//  October 2005
//
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
//
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
// 
//


#ifndef USE_MS2_FRAGMENT_H
#define USE_MS2_FRAGMENT_H



class MS2Fragment{

    
    ////////////////////////////////////////////////
    // declaration of the private members:

private:
  
  ////////////////////////////////////////////////
  // declaration of the public members:
  
  // AMRT tag
  double precursorMZ;
  int precursorCHRG;
  double TR;
  int scan;
  int z;
  
  double fragmentMZ;
  double intensityArea;
  
  // scan and TR ranges:
  int scanStart;
  int scanEnd;
  double trStart;
  double trEnd;
  
  
public:
    
    static int OutlierAttribute;
  
  // class destructor
  ~MS2Fragment();
  
  // constructor for the object MS2Fragment:
  MS2Fragment(double iPrecursorMZ, int iPrecursorCHRG, double iTR, int iScan, int iZ, double iFragmentMZ, double iIntensityArea,
                           int iScanStart, int iScanEnd, double iTrStart, double iTrEnd);    
  MS2Fragment(double iPrecursorMZ, int iPrecursorCHRG, double iTR, int iScan, int iZ, double iFragmentMZ, double iIntensityArea );

  
    // class constructor
  MS2Fragment();
    // class copy constructor
  MS2Fragment(const MS2Fragment&);
  // class copy constructor
  MS2Fragment(const MS2Fragment*);
  
  // show info of the MS2 fragment
  void show_info();

  // get the attribute of the fragment
  // according to which outliers are removed
  double getOutlierDetectionAttribute();
    
  
  //////////////////////////////////////////////////
  // overload operators:
  MS2Fragment& operator=(const MS2Fragment&);
  bool operator==(const MS2Fragment&);
  MS2Fragment& operator<=(const MS2Fragment&);
  MS2Fragment& operator>=(const MS2Fragment&);
  MS2Fragment& operator<(const MS2Fragment&);
  MS2Fragment& operator>(const MS2Fragment&);
  
  
  ///////////////////////////////
  // start here all the get / set
  // function to access the
  // variables of the class
  
  
  // get hte averaged precurso mass:
  double getPrecursorMZ(){ return precursorMZ;};
  void setPrecursorMZ(double iMZ){ precursorMZ=iMZ;};
  // get hte averaged precurso chrg:
  int getPrecursorCHRG(){ return precursorCHRG;};
  // retention time:
  double getTR(){return TR;};
  // start TR:
  double getStartTR(){return trStart;};
  // end TR:
  double getEndTR(){return trEnd;};
  // get the Fragment MZ:
  double getFragmentMz(){return fragmentMZ;};
  void setFragmentMz(double iMz){fragmentMZ=iMz;};
  // get teh charge state:
  int getCHRG(){return z;};
  // get the apex scan:
  int getApexScan(){return scan;};
  // get the apex scan:
  int getStartScan(){return scanStart;};
  // get the apex scan:
  int getEndScan(){return scanEnd;};
  
  // get the integrated peak area:
  double getFragmentPeakArea(){return intensityArea;};
  void setFragmentPeakArea(double iIntens){intensityArea=iIntens;};


};

#endif

    
