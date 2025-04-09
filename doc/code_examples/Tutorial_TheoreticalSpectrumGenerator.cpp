// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

//! [doxygen_snippet_TSG]

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

int main()
{
  // initialize a TheoreticalSpectrumGenerator
  TheoreticalSpectrumGenerator tsg;

  // get current parameters
  // in this case default parameters, since we have not changed any yet
  Param tsg_settings = tsg.getParameters();
    
  // with default parameters, only b- and y-ions are generated,
  // so we will add a-ions
  tsg_settings.setValue("add_a_ions", "true");
    
  // store ion types for each peak
  tsg_settings.setValue("add_metainfo", "true");
    
  // set the changed parameters for the TSG
  tsg.setParameters(tsg_settings);
                     

  PeakSpectrum theoretical_spectrum;

  // initialize peptide to be fragmented
  AASequence peptide = AASequence::fromString("DEFIANGER");
  
  // generate a-, b- and y- ion spectrum of the peptide
  // with all fragment charges from 1 to 2
  tsg.getSpectrum(theoretical_spectrum, peptide, 1, 2);
 
  // output of masses and meta information (ion-types) of some peaks
  const PeakSpectrum::StringDataArray& ion_types = theoretical_spectrum.getStringDataArrays().at(0);
  cout << "Mass of second peak: " << theoretical_spectrum[1].getMZ()
       << " | Ion type of second peak: " << ion_types[1] << endl;

  cout << "Mass of tenth peak: " << theoretical_spectrum[9].getMZ()
       << " | Ion type of tenth peak: " << ion_types[9] << endl;
  
  return 0;
} //end of main

//! [doxygen_snippet_TSG]
