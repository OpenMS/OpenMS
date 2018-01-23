// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

//! [TSG]

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

//! [TSG]
