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
// $Authors: Florian Zeller $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerSH.h>

#include <vector>

namespace OpenMS
{
  PeakPickerSH::PeakPickerSH()
    : DefaultParamHandler("PeakPickerSH"),
      ProgressLogger()
  {
    defaultsToParam_();
  }
  
  PeakPickerSH::~PeakPickerSH()
  {
    // FLO: Do not care at the moment
  }
  
  template <typename PeakType>
  void PeakPickerSH::pick(const MSSpectrum<PeakType>& input, MSSpectrum<PeakType>& output, float fWindowWidth)
  {
    int i, hw, j;
    double cm, toti, min_dh;
    
    // Hack: Prepare data structures for Lukas' algorithm
    std::vector<double> masses,intens;
    // TODO: Probably we could save some time when we resize the vectors...
    //masses.resize(input.size());
    //intens.resize(input.size());
    for (Size k = 0; k < input.size()-1; ++k)
    {
      // Lukas requires a minimum of intensity (=50). His vectors do not contain
      // other data, so I strip the low ones out right here.
      // TODO: Read 50.0 from parameters
      if (input[k].getIntensity() >= 50.0)
      {
        masses.push_back(input[k].getMZ());
        intens.push_back(input[k].getIntensity());
      }
    }
    
    min_dh = 50.0;				// min height
    hw = fWindowWidth/2;
    
    for (i=2;i<(int)masses.size()-2;i++) { 
      
      // Peak must be concave in the interval [i-2 .. i+2]
      if (intens[i]>min_dh && intens[i]>intens[i-1]+min_dh && intens[i]>=intens[i+1] && intens[i-1]>intens[i-2]+min_dh && intens[i+1]>=intens[i+2]) {
        
        cm = 0.0;					// centroid mass:
        toti = 0.0;				// total intensity:
        
        for (j= -hw;j<=hw;j++) {
          double inte = intens[i-j];
          double mz = masses[i-j];
          
          cm += inte*mz;
          toti += (double) intens[i-j];
        }
        cm = cm/toti;			// Centre of gravity = centroid
        
        PeakType peak;
        peak.setMZ(cm);
        peak.setIntensity(intens[i]);
        output.push_back(peak);
      }
    }
  }
  
  void PeakPickerSH::pickExperiment(const MSExperiment<>& input, MSExperiment<>& output)
  {
    // make sure that output is clear
    output.clear(true);
    
    // copy experimental settings
    static_cast<ExperimentalSettings&>(output) = input;
    
    // resize output with respect to input
    output.resize(input.size());
    
    std::cout << "Before loop, input size = " << input.size() << std::endl;
    Size progress = 0;
    for (Size scan_idx = 0; scan_idx != input.size(); ++scan_idx)
    {	
      output[scan_idx].clear(true);
      output[scan_idx].SpectrumSettings::operator=(input[scan_idx]);
      output[scan_idx].MetaInfoInterface::operator=(input[scan_idx]);
      output[scan_idx].setRT(input[scan_idx].getRT());
      output[scan_idx].setMSLevel(input[scan_idx].getMSLevel());
      output[scan_idx].setName(input[scan_idx].getName());
      output[scan_idx].setType(SpectrumSettings::PEAKS);
      
      if (input[scan_idx].getMSLevel() != 1)
      {
        // When not considering MS2 data (MS2 fragment mass tracing=0), Lukas leaves out
        // the entire scan (instead of just copying it to the output as seen in 
        // another plugin).
        // pick(input[scan_idx], output[scan_idx], 4.0);
      }
      else
      {
        // TODO: Read value 4.0 from parameters
        pick(input[scan_idx], output[scan_idx], 5.0);
      }
      setProgress(++progress);
    }
    std::cout << "After loop" << std::endl;
    
    endProgress();
  }
}
