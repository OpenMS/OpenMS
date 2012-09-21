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

#ifndef USE_BACKGROUND_INTENSITY_BIN_H
#define USE_BACKGROUND_INTENSITY_BIN_H

namespace OpenMS
{

class OPENMS_DLLAPI BackgroundIntensityBin {

private:

  // mz and tr coordinates of the bin:
  double mzCoord;
  double trCoord;
  double zCoord;  
  
  std::vector<double> IntensityMap;  
  std::map<double, double> IntensityHist;

  double mean;
  double median;
  void computeIntensityHist();

public:
  
  static double TR_BINS;
  static double MZ_BINS;
  static double INTENS_BINS;
  static int MIN_BIN_COUNT;
  
  ~BackgroundIntensityBin();

  BackgroundIntensityBin(double, double);
  
  // check if a peak belongs to this intenisty bin
  bool checkBelonging(ms_peak*);
  // add intensity to BackgroundIntensityBin
  void addIntensity( double );
  // add peak to BackgroundIntensityBin
  void addMSPeak( ms_peak* );
  // process collected intensities in the map
  void processIntensities( );

  std::vector<double>* getIntensityMap(){ return &IntensityMap;};
  std::map<double, double>* getIntensityHist(){ return &IntensityHist;};
  double getMean(){return mean;};
};

} // ns

#endif

    
