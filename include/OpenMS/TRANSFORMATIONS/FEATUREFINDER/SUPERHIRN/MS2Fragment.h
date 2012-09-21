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

#include <OpenMS/CONCEPT/Types.h>

namespace OpenMS
{

class OPENMS_DLLAPI MS2Fragment{

    
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

} // ns

#endif

    
