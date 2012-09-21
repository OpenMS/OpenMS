// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Florian Zeller $
// $Authors: Lukas Mueller, Markus Mueller $
// --------------------------------------------------------------------------
//
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

#ifndef MS1_FEATURE_MERGER_H
#define MS1_FEATURE_MERGER_H

namespace OpenMS
{

class OPENMS_DLLAPI MS1FeatureMerger{

    
  ////////////////////////////////////////////////
  // declaration of the private members:
  
private:
  
  //////
  // the common lc-ms spectrum,
  // created from teh overlap of A and B
  LCMS* lcmsMap;
  
  std::vector<int> idsToRemove;
  std::map< double, std::vector<SHFeature*> > mzClusters;

  
  ////////////////////////////////////////////////
  // declaration of the public members:
  
public:
  
  /*
  static double INTENSITY_APEX_THRESHOLD;
  static double MS1_PEAK_AREA_TR_RESOLUTION;

  static double INITIAL_TR_TOLERANCE;
  static double MS1_FEATURE_MERGING_TR_TOLERANCE;
  static double PERCENTAGE_INTENSITY_ELUTION_BORDER_VARIATION;
  static double PPM_TOLERANCE_FOR_MZ_CLUSTERING;
  static bool MS1_FEATURE_CLUSTERING;
*/
  
  // class destructor
  ~MS1FeatureMerger();
  // class constructor
  MS1FeatureMerger(LCMS*);
  
  
  //////////////////////////////////////////////////
  // start the merging process
  void startFeatureMerging();
  // create a distribution of delta Tr for the splited features
  void createMZFeatureClusters();
  // process a vector of m/z features
  void processMZFeatureVector(std::vector<SHFeature*>* );
  // find to this feature the features which should be merged
  std::vector<SHFeature*>::iterator findFeaturesToMerge(SHFeature* , std::vector<SHFeature*>::iterator , std::vector<SHFeature*>* );
  // compare if a feature belongs to another feature
  bool compareMZFeatureBeloning(SHFeature* , SHFeature* );
  // merge the target to the search feature
  void mergeFeatures( SHFeature* , SHFeature* );
  // copmute new parameters for the merged MS1 feature
  void computeNewMS1FeatureParameters( SHFeature* );
  // computes the area of between 2 peaks:
  double computeDeltaArea(double, double, double, double);

  
  // this structure provides the function to compare
  // in the sorting algorithm:
  struct OPERATOR_FEATURE_TR{
    // provide the compare function for sort:
    bool operator()(const SHFeature A,const SHFeature B) const{
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
