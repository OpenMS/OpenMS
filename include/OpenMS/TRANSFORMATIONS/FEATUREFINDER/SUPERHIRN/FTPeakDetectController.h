// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_SUPERHIRN_FTPEAKDETECTCONTROLLER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_SUPERHIRN_FTPEAKDETECTCONTROLLER_H

#include <boost/shared_ptr.hpp>

namespace OpenMS
{

  class OPENMS_DLLAPI FTPeakDetectController
  {

    ////////////////////////////////////////////////
    // declaration of the private members:

private:

    ////////////////////////////////////////////////
    // declaration of the public members:

    // LCMS runs
    //LCMS* THIS_LCMS;
    LCMS * lcms_;
//  std::vector<SHFeature> fakeFeatureList_;
    std::vector<LCMS> lcmsRuns_;

    // paths:
    std::string targetMzXML;
    std::string SOURCE_DIR;
    std::string OUTPUT_DIR;

public:

    typedef std::pair<double, boost::shared_ptr<RawData> > Map;
    typedef std::vector<Map> Vec;

//  static bool CREATE_FEATURE_ELUTION_PROFILES;
//  static bool LCelutionPeakDebugging;
//  static double LCelutionPeakMassMin;
//  static double LCelutionPeakMassMax;

//  static MS2Feature* SearchedM2Feature;

//  static bool FEATURE_FAKE_INSERTION_BASED_ON_MS2_FEATURE;

    // class destructor
    ~FTPeakDetectController();

    // class constructor
    FTPeakDetectController();
    // class copy constructor
    FTPeakDetectController(const FTPeakDetectController &);

    /////////////////////////////////////////////////////////
    // function for batch processing of mzXML data
    // parses LC-MS from runs from a directory or file of raw mzXML data:
    void parseMzXMLData();

    //////////////////////////////////////////////////
    // mzXML parsing functions for a single MzXML file:
    // start the scan parsing of a mzXML file:
    void startScanParsing(Vec datavec);

    // **** for the MS1 level post processing:
    // process MS1 level data
    void process_MS1_level_data_structure(ProcessData *);
    // adds an elution peak to the LC/MS run:
    void add_raw_peak_to_LC_MS_run(LCElutionPeak *);
    // function to add the elution profile to the feature:
    void addLCelutionProfile(SHFeature *, LCElutionPeak *);

    /////////////////////////////////////////////////////////////
    // reads already paths of existing LC-MS runs in xml format into the
    // memory
    // for now, open file system for every check, but otherwise could eb done
    // in the constructor
    bool checkIfFeatureExtractionExists(std::string);

    // **** for the MS2 level post processing:
    // process MS2 level data
    void process_MS2_level_data_structure(ProcessData *);
    // processes the tracted signals on teh MS2 level
    void extract_MS2_elution_features();
    // combine the MS2 feature trace data to the MS1 features:
    void associateMS2FeatureToMS1Feature(MS2Feature *);
    // add an observed MS2 feature to the MS1 feature
    // if an observation is already there, then
    // construct a merged MS2 feature
    void addMS2FeatureToMS1Feature(MS2Feature *, SHFeature *);

    // construct here fake ms1 features based on a observed MS2 feature
    // which however could not be matched to a exiting ms1 feature
    void constructMS1FeatureFromMS2Feature(MS2Feature *);

    /////////////////////////////////////////////////////////
    // write a parsed LC/MS into directory:
    void write_out_parsed_LC_MS(LCMS *);
    // add fake MS/MS information for the MS1 feature:
    void addFakeMSMSToFeature(SHFeature *);


    //////////////////////////////////////////////////
    // overload operators:
    FTPeakDetectController & operator=(const FTPeakDetectController &);
    FTPeakDetectController & operator<=(const FTPeakDetectController &);
    FTPeakDetectController & operator>=(const FTPeakDetectController &);
    FTPeakDetectController & operator<(const FTPeakDetectController &);
    FTPeakDetectController & operator>(const FTPeakDetectController &);

    ///////////////////////////////
    // start here all the get / set
    // function to access the
    // variables of the class

    // target file:
    void set_target_file(std::string IN);
    std::string get_target_file();

// get the vector of LC/MS runs:
//std::vector<LCMS> getParsedData();
//bool getParsedDataEmpty();
//std::vector<LCMS>::iterator get_parsed_DATA_START();
//std::vector<LCMS>::iterator get_parsed_DATA_END();
    LCMS * getLCMS();
  };

  inline void FTPeakDetectController::set_target_file(std::string IN)
  {
    targetMzXML = IN;
  }

  inline std::string FTPeakDetectController::get_target_file()
  {
    return targetMzXML;
  }

  /*
  // get the vector of LC/MS runs:
  inline std::vector<LCMS> FTPeakDetectController::getParsedData()
  {
      return lcmsRuns_;
  }

  inline bool FTPeakDetectController::getParsedDataEmpty()
  {
      return LC_MS_RUNS.empty();
  }

  inline std::vector<LCMS>::iterator FTPeakDetectController::get_parsed_DATA_START()
  {
      return LC_MS_RUNS.begin();
  }

  inline std::vector<LCMS>::iterator FTPeakDetectController::get_parsed_DATA_END()
  {
      return LC_MS_RUNS.end();
  }
*/

  inline LCMS * FTPeakDetectController::getLCMS()
  {
    return lcms_;
  }

} // ns

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_SUPERHIRN_FTPEAKDETECTCONTROLLER_H
