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

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

class TICWritingConsumer : public MSDataWritingConsumer 
{
  // Inheriting from MSDataWritingConsumer allows to change the data before
  // they are written to disk (to "filename") using the processSpectrum_ and
  // processChromatogram_ functions.
public:
  double TIC;
  int nr_spectra;

  // Create new consumer, set TIC to zero
  TICWritingConsumer(String filename) : MSDataWritingConsumer(filename) 
    { TIC = 0.0; nr_spectra = 0;}

  // Add a data processing step for spectra before they are written to disk
  void processSpectrum_(MSDataWritingConsumer::SpectrumType & s) override
  {
    for (Size i = 0; i < s.size(); i++) { TIC += s[i].getIntensity(); }
    nr_spectra++;
  }
  // Empty chromatogram data processing
  void processChromatogram_(MSDataWritingConsumer::ChromatogramType& /* c */) override {}
};

int main(int argc, const char** argv)
{
  if (argc < 2) return 1;
  // the path to the data should be given on the command line
  String tutorial_data_path(argv[1]);
  
  // Create the consumer, set output file name, transform
  TICWritingConsumer * consumer = new TICWritingConsumer("Tutorial_FileIO_output.mzML");
  MzMLFile().transform(tutorial_data_path + "/data/Tutorial_FileIO_indexed.mzML", consumer);

  std::cout << "There are " << consumer->nr_spectra << " spectra in the input file." << std::endl;
  std::cout << "The total ion current is " << consumer->TIC << std::endl;
  delete consumer;

  return 0;
} //end of main
