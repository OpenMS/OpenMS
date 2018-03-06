//---------------------------------------------------------------------------
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

//! [Precursor]

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

int main(int argc, const char** argv)
{
  
  if (argc < 2) return 1;
  
  // the path to the data should be given on the command line
  String tutorial_data_path(argv[1]);
  
  MSExperiment spectra;
  MzMLFile f;

  // load mzML from code examples folder
  f.load(tutorial_data_path + "/data/Tutorial_GaussFilter.mzML", spectra);

  // iterate over map and output MS2 precursor information
  for (auto s_it = spectra.begin(); s_it != spectra.end(); ++s_it)
  {
    // we are only interested in MS2 spectra so we skip all other levels
    if (s_it->getMSLevel() != 2) continue;

    // get a reference to the precursor information
    const MSSpectrum& spectrum = *s_it;
    const vector<Precursor>& precursors = spectrum.getPrecursors();

    // size check & throw exception if needed 
    if (precursors.empty()) throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, precursors.size());

    // get m/z and intensity of precursor
    double precursor_mz = precursors[0].getMZ();
    float precursor_int = precursors[0].getIntensity();
  
    // retrieve the precursor spectrum (the most recent MS1 spectrum)
    PeakMap::ConstIterator precursor_spectrum = spectra.getPrecursorSpectrum(s_it);
    double precursor_rt = precursor_spectrum->getRT();
  
    // output precursor information
    std::cout << " precusor m/z: " << precursor_mz
              << " intensity: " << precursor_int
              << " retention time (sec.): " << precursor_rt 
              << std::endl;
   }
                                                            
  return 0;
} // end of main

//! [Precursor]
