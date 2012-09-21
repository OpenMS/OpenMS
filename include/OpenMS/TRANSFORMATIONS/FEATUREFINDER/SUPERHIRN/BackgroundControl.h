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


#ifndef USE_BACKGROUND_CONTROL_H
#define USE_BACKGROUND_CONTROL_H



class BackgroundControl{

    
    ////////////////////////////////////////////////
    // declaration of the private members:

private:

  
  map<double, map< double, BackgroundIntensityBin> > intensityBinMap;
  
  

  ////////////////////////////////////////////////
  // declaration of the public members:
  
public:
  
  // class destructor
  ~BackgroundControl();

  // class constructor
  BackgroundControl();
  // class copy constructor
  BackgroundControl(const BackgroundControl&);
  // class copy constructor
  BackgroundControl(const BackgroundControl*);
  
  
  ////////////////////////////////////////////////
  // initialization of the BackgroundControl classe:
  void init();
    
  //////////////////////////////////////////////////
  // add a peak to the BackgroundControl
  void addPeak( ms_peak* );
  // add peaks of a ms scan:
  void addPeakMSScan( double , vector<ms_peak>* );
  // add peaks of a ms scan:
  void addPeakMSScan( double , list<CentroidPeak>* );

  
  //////////////////////////////////////////////////
  // get the background intensity level for a peak:
  double getBackgroundLevel( ms_peak* );
  double getBackgroundLevel( double mz, double tr);
  
  // find a key in the intensity map:
  map<double, map< double, BackgroundIntensityBin>  >::iterator findTrKey( double );
  // find a key in the m/z map:
  map< double, BackgroundIntensityBin>::iterator findMzKey( double mz, map< double, BackgroundIntensityBin>* );

    
  //////////////////////////////////////////////////
  // process the intensity maps:
  void processIntensityMaps(  );    
  // write out intensity maps:
  //void writeIntensityMaps( );
  // gnuplot intensity maps:
  //void plotIntensityMaps(  );

    
  //////////////////////////////////////////////////
  // overload operators:
  BackgroundControl& operator=(const BackgroundControl&);
  bool operator==(const BackgroundControl&);
  BackgroundControl& operator<=(const BackgroundControl&);
  BackgroundControl& operator>=(const BackgroundControl&);
  BackgroundControl& operator<(const BackgroundControl&);
  BackgroundControl& operator>(const BackgroundControl&);
  
  
  ///////////////////////////////
  // start here all the get / set
  // function to access the
  // variables of the class
  
};

#endif

    
