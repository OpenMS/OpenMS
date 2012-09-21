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
//  PEAK DETECTION OF FOURIER TRANSFORME MS INSTRUMENT DATA
//
//  by Lukas Mueller, Lukas.Mueller@imsb.biol.ethz.ch
//  October 2005
//  
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
//
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
// 
//

#include <fstream>
#include <iostream>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/RawData.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/ms_peak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/CentroidPeak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/CentroidData.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/LC_elution_peak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/BackgroundIntensityBin.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/PeptideIsotopeDistribution.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/ExternalIsotopicDistribution.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/BackgroundControl.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/Deisotoper.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/LCMSCData.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/simple_math2.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/Process_Data.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/FT_PEAK_DETEC_mzXML_reader.h>

namespace OpenMS
{

// debugging classes:
int	FT_PEAK_DETEC_mzXML_reader::sfReportMonoPeaks = 0; // 1 if info about monoisotopic peaks should be written to mono_peaks.txt
std::string FT_PEAK_DETEC_mzXML_reader::sfDebugDirectory; // Directory where peak detection debug files are written
int	FT_PEAK_DETEC_mzXML_reader::sfReportScanNumber = -1; // if sfReportMonoPeaks is set to 1, details about this spectrum will be written to debug files 

std::vector<double> FT_PEAK_DETEC_mzXML_reader::FRAGMENT_MASS_SCAN_LEVELS;
std::vector<double> FT_PEAK_DETEC_mzXML_reader::PEAK_EXTRACTION_SCAN_LEVELS;
bool FT_PEAK_DETEC_mzXML_reader::MS2_PEAK_PROCESSING = false;
int FT_PEAK_DETEC_mzXML_reader::MS1_base_inter_scan_distance;
int FT_PEAK_DETEC_mzXML_reader::MS2_base_inter_scan_distance;
double FT_PEAK_DETEC_mzXML_reader::TR_MIN;
double FT_PEAK_DETEC_mzXML_reader::TR_MAX;

////////////////////////////////////////////////
// constructor for the object ana_summarizer:
FT_PEAK_DETEC_mzXML_reader::FT_PEAK_DETEC_mzXML_reader(){
  
  MS1_LC_MS_DATA_PROCESSOR = NULL;
  //MS2_LC_MS_DATA_PROCESSOR = NULL;
  
  // initialize the variables
  index_offset = 0;
  scan_index = NULL;
  total_scan = 0;
  nbMS2Scans = 0;
  
  // FLOFLO
  MS1_LC_MS_DATA_PROCESSOR = new Process_Data( );

}

//////////////////////////////////////////////////
// class desctructor
FT_PEAK_DETEC_mzXML_reader::~FT_PEAK_DETEC_mzXML_reader(){
  if(MS1_LC_MS_DATA_PROCESSOR != NULL ){
    delete MS1_LC_MS_DATA_PROCESSOR;
    MS1_LC_MS_DATA_PROCESSOR = NULL;
  }
}

//////////////////////////////////////////////////
// set indexes of the current mzXML file;
void FT_PEAK_DETEC_mzXML_reader::set_current_indexes(double pminrt, double pmaxrt){
  
    // get min / max retention time in mzXML:
    // read the mzXML header at this scan
    //struct ScanHeaderStruct scan_header;
    // start:
    //readHeader( get_Ramp_file_handler(), get_scan( 0 ), &scan_header);
    //minRT = scan_header.retentionTime / 60.00;
  minRT = pminrt;
    //readHeader( get_Ramp_file_handler(), get_scan( total_scan ), &scan_header);
    //maxRT = scan_header.retentionTime / 60.00;
  maxRT = pmaxrt;
  
  // FLOFLO
  std::cout << "minRT = " << minRT << ", maxRT = " << maxRT << "\n";

    ExternalIsotopicDistribution::initRetentionTimeSegments( minRT, maxRT );

}





//////////////////////////////////////////////////
// reads the ms data from a mzXML file opened by teh handler
void FT_PEAK_DETEC_mzXML_reader::read_mzXML_DATA(Vec datavec){

  unsigned int i;
  double minrt = datavec[0].begin()->first;
  double maxrt = datavec[datavec.size()-1].begin()->first;
  set_current_indexes(minrt, maxrt);
  
  std::cout << "Anzahl scans: " << datavec.size() << std::endl;
  
  std::map<double, RawData*>::const_iterator it;
  for (i=0; i<datavec.size(); i++)
  {
    
    Map::iterator it = datavec.at(i).begin();
    
    // debug code to compare inputs in CleanUpSH2 and Superhirn
    RawData* pRawData = it->second;
    std::vector<double> masses,intens;
    pRawData->get(masses,intens);
    get_MS_scan(i, it->first, it->second);
  }
  
  std::cout << "Number of scans: " << i << "\n";
    
  /*
  // check all the scans of the current region:
  off_t scan = 0;
  off_t MAX_SCAN = total_scan;
  
  while(scan <= MAX_SCAN){
    // extract ms info within a mass range of a scan
    get_MS_scan(scan);
    scan++;
  }
   */
}

//////////////////////////////////////////////////
// get a MS scan at a given scan number within
// a mass range
void FT_PEAK_DETEC_mzXML_reader::get_MS_scan(off_t IN, double TR, RawData* data){
  
    ///////////////////////
    // get retention time:
    //TR = scan_header.retentionTime / 60.00;
  
    /////////////////////////////////////
    // check if this scan (MSx) is in the TR range
    // defned by max / min TR:
    // should have some peaks!
    if( (TR >= TR_MIN) && ( TR <= TR_MAX )){  //  &&( scan_header.peaksCount > 1 ) 
      
      // print the scan header:
      // print_scan_header( &scan_header );
      
      // build up an index scan vs retention time:
      insert_into_scan_TR_index(IN, (float)TR);
      
      //////////////////////////////////////////////////////
      // check here the MS Precursor Mass Spectrum Level 
      //if( checkMSPrecursorMassScan( scan_header.msLevel ) ){
        
        // set the maximal inter-monoisotopic distance
        // for the same LC-elution peak
        // set the inter-monoistopic distance:
        //int max_scan = setInterMonoIsotopicLCDistance( IN , 1, FT_PEAK_DETEC_mzXML_reader::MS1_base_inter_scan_distance);
        int max_scan = 0;
        MS1_LC_MS_DATA_PROCESSOR->setMaxScanDistance( max_scan);
        
        // process the data:
        processMS1InputData(IN, (float)TR, data);
        
      //}
  }

}


///////////////////////////////////////////////////////////////////////////////////////
// process the MS1 level input data:
// - construct a RawData object with input peaks
// - centoid them
// - add them to Process_Data Structure:
void FT_PEAK_DETEC_mzXML_reader::processMS1InputData(int SCAN, float TR, RawData* data){
  
  //////////////////////////////////
  // bool debug = (FT_PEAK_DETEC_mzXML_reader::sfReportMonoPeaks ==1 && SCAN== FT_PEAK_DETEC_mzXML_reader::sfReportScanNumber);
//  IsotopicDist::setDebug(false);
  
  //////////////////////////////////
  // construct a raw_data
  //RawData data(ms_peaks, CentroidData::sfIntensityFloor);
  
  //////////////////////////////////
  // centroid it:
  //CentroidData cd(CentroidPeak::sfCentroidWindowWidth, data, Process_Data::CENTROID_DATA_MODUS);
  CentroidData cd(CentroidPeak::sfCentroidWindowWidth, *data, TR, Process_Data::CENTROID_DATA_MODUS);
  
  //////////////////////////////////
  //  store it:
  MS1_LC_MS_DATA_PROCESSOR->add_scan_raw_data(SCAN, TR, &cd);
  
}

/////////////////////////////////////////////////////////////////////////
// set the maximal inter-monoisotopic distance
// for the same LC-elution peak
//int FT_PEAK_DETEC_mzXML_reader::setInterMonoIsotopicLCDistance(int my_scan, int level, int max_inter_scan_distance){
  
  
  
  // TODO: Rewrite, but setInterMonoIsotopicLCDistance in Superhirn returns always 0, too
  //printf("setInterMonoIsotopicLCDistance is FAKE\n");
//  return 0;
  
  
  
  /*
  if( my_scan == 0 ){
    return 0;
  }
  
  // check how many ms scan of not this level have been 
  // performed since the last MS scan on the input level:
  struct ScanHeaderStruct scan_header;
  int tmp_scan = my_scan - 1;
  int nonMSLevelCount = 0;
  int MSLevelCount = 0;
  readHeader( get_Ramp_file_handler(), get_scan( tmp_scan ), &scan_header);
  while( MSLevelCount < max_inter_scan_distance){
    
    // count the MS scans of this level:
    nonMSLevelCount++;

    // count the MS scans of this level:
    if( scan_header.msLevel == level){
      MSLevelCount++;
    }     
    tmp_scan--;
    
    // dont go below 1
    if( tmp_scan < 1){
      break;
    }
    
    readHeader( get_Ramp_file_handler(), get_scan(tmp_scan), &scan_header);
  }
  
  return nonMSLevelCount;
   */
//}

}
