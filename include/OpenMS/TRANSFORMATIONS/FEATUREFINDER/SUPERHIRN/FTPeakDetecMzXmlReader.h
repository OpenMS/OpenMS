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

namespace OpenMS
{

// file structure remapping
// for the mzXML parser ramp:
// - for old ramp use FILE*
// - for new cramp use RAMPFILE*
//typedef FILE* MZXML_FILE;
//typedef RAMPFILE RAMP_FILE;



class OPENMS_DLLAPI FTPeakDetecMzXmlReader{

    
    ////////////////////////////////////////////////
    // declaration of the private members:
    
private:
  
  ProcessData* dataProcessor_;

  
  ////////////////////////////////////////////////
  // declaration of the public members:
  
public:
  
  
  typedef std::map<double, RawData*> Map;
  typedef std::vector<Map> Vec;

  // class destructor
  ~FTPeakDetecMzXmlReader();
  // class constructor
  FTPeakDetecMzXmlReader();
    
  // reads the ms data 
  void readData(Vec datavec);

  // get a MS scan at a given scan number within
  // a mass range and stores it into the ms region structure
  void getMsScan(off_t IN, double TR, RawData* data);

  
  ///////////////////////////////////////////////////////
  // MS1 level processing routine functions:
  // process the MS1 level input data:
  // - construct a RawData object with input peaks
  // - centoid them
  // - add them to Process_Data Structure:
  void processMS1InputData(int , float, RawData* data);
     

  ///////////////////////////////////////////////////////////////////////////////
  // check if the scan number should be processed by MS Precursor Mass extraction
  bool checkMSPrecursorMassScan( int );
  // check if the scan number should be processed as MSn FragmentMass spectrum
  bool checkMSFragmentMassScan( int );

  
  // get the MS1 processed data:
  ProcessData* getData();

};

	inline ProcessData* FTPeakDetecMzXmlReader::getData() 
	{
		return dataProcessor_;
	}

} // ns

#endif

