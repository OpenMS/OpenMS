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
// CLASS Spec MERGER:
// merges 2 spectra which have been preprocessed
// by spec_merge-> common lc_peaks are marked!!!!
//
// variable description:
//
//
// function description:
//
//
// **********************************************************************//


#ifndef MS1_FEATURE_MERGER_H
#define MS1_FEATURE_MERGER_H

namespace OpenMS
{

class OPENMS_DLLAPI MS1_feature_merger{

    
  ////////////////////////////////////////////////
  // declaration of the private members:
  
private:
  
  //////
  // the common lc-ms spectrum,
  // created from teh overlap of A and B
  LC_MS* lcmsMap;
  
  std::vector<int> idsToRemove;
  std::map< double, std::vector<feature*> > mzClusters;

  
  ////////////////////////////////////////////////
  // declaration of the public members:
  
public:
  
  static double INTENSITY_APEX_THRESHOLD;
  static double MS1_PEAK_AREA_TR_RESOLUTION;

  static double INITIAL_TR_TOLERANCE;
  static double MS1_FEATURE_MERGING_TR_TOLERANCE;
  static double PERCENTAGE_INTENSITY_ELUTION_BORDER_VARIATION;
  static double PPM_TOLERANCE_FOR_MZ_CLUSTERING;
  static bool MS1_FEATURE_CLUSTERING;
  
  // class destructor
  ~MS1_feature_merger();
  // class constructor
  MS1_feature_merger(LC_MS*);
  
  
  //////////////////////////////////////////////////
  // start the merging process
  void startFeatureMerging();
  // create a distribution of delta Tr for the splited features
  void createMZFeatureClusters();
  // process a vector of m/z features
  void processMZFeatureVector(std::vector<feature*>* );
  // find to this feature the features which should be merged
  std::vector<feature*>::iterator findFeaturesToMerge(feature* , std::vector<feature*>::iterator , std::vector<feature*>* );
  // compare if a feature belongs to another feature
  bool compareMZFeatureBeloning(feature* , feature* );
  // merge the target to the search feature
  void mergeFeatures( feature* , feature* );
  // copmute new parameters for the merged MS1 feature
  void computeNewMS1FeatureParameters( feature* );
  // computes the area of between 2 peaks:
  double computeDeltaArea(double, double, double, double);

  
  // this structure provides the function to compare
  // in the sorting algorithm:
  struct OPERATOR_FEATURE_TR{
    // provide the compare function for sort:
    bool operator()(const feature A,const feature B) const{
      // check if they have same mass
      return A.TR < B.TR;
    }
  };
  
  
  ///////////////////////////////
  // start here all the get / set
  // function to access the
  // variables of the class
 
};

} // ns

#endif
