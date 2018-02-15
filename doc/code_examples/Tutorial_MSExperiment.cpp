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

//! [MSExperiment]

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

int main()
{

  // create a peak map containing 4 dummy spectra and peaks
  MSExperiment exp;

  // The following examples creates a MSExperiment containing four MSSpectrum instances.
  for (Size i = 0; i < 4; ++i)
  {
    MSSpectrum spectrum;
    spectrum.setRT(i);
    spectrum.setMSLevel(1);
    for (float mz = 500.0; mz <= 900; mz += 100.0)
    {
      Peak1D peak;
      peak.setMZ(mz + i);
      spectrum.push_back(peak);
    }
    
    exp.addSpectrum(spectrum);
  }

  // Iteration over the RT range (2,3) and the m/z range (603,802) and print the peak positions.
  for (auto it = exp.areaBegin(2.0, 3.0, 603.0, 802.0); it != exp.areaEnd(); ++it)
  {
    cout << it.getRT() << " - " << it->getMZ() << endl;
  }

  // Iteration over all peaks in the experiment. 
  // Output: RT, m/z, and intensity
  // Note that the retention time is stored in the spectrum (not in the peak object)
  for (auto s_it = exp.begin(); s_it != exp.end(); ++s_it)
  {
    for (auto p_it = s_it->begin(); p_it != s_it->end(); ++p_it)
    {
      cout << s_it->getRT() << " - " << p_it->getMZ() << " " << p_it->getIntensity() << endl;
    }
  }

  // We could store the spectra to a mzML file with:
  // MzMLFile mzml;
  // mzml.store(filename, exp);
  
  // And load it with
  // mzml.load(filename, exp);
  // If we wanted to load only the MS2 spectra we could speed up reading by setting:
  // mzml.getOptions().addMSLevel(2);
  // before executing: mzml.load(filename, exp);

  return 0;
} //end of main

//! [MSExperiment]
