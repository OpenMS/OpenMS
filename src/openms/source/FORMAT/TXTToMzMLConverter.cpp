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
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/TXTToMzMLConverter.h>

namespace OpenMS
{
  TXTToMzMLConverter::TXTToMzMLConverter() {}

  TXTToMzMLConverter::~TXTToMzMLConverter() {}

  MSExperiment TXTToMzMLConverter::loadInputFile(const String& filename) const
  {
    std::ifstream ifs(filename, std::ifstream::in);
    if (!ifs.is_open())
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
    const Size BUFSIZE = 65536;
    char line[BUFSIZE];
    bool header_found = false;
    String header("Time (min)\tStep (s)\tValue (mAU)");
    while (!header_found && !ifs.eof())
    {
      ifs.getline(line, BUFSIZE);
      String s(line);
      if (std::strncmp(header.c_str(), line, header.size()) == 0)
      {
        header_found = true;
      }
    }
    MSChromatogram chromatogram;
    while (!ifs.eof())
    {
      ifs.getline(line, BUFSIZE);
      double rt, intensity;
      int ret = std::sscanf(line, "%lf\t%*f\t%lf", &rt, &intensity);
      if (ret == 2)
      {
        chromatogram.push_back(ChromatogramPeak(rt, intensity));
      }
      else if (!strcmp(line, "\r") || !strcmp(line, ""))
      {
        break;
      }
      else
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, std::string(line), "Couldn't parse the raw data.");
      }
    }
    MSExperiment experiment;
    experiment.addChromatogram(chromatogram);
    return experiment;
  }

  void TXTToMzMLConverter::storeMzMLFile(const String& filename, const MSExperiment& experiment) const
  {
    MzMLFile mzml;
    mzml.store(filename, experiment);
  }
}
