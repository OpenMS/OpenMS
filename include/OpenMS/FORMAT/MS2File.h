// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_MS2FILE_H
#define OPENMS_FORMAT_MS2FILE_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/METADATA/DocumentIdentifier.h>

#include <vector>
#include <fstream>

namespace OpenMS
{
  /**
      @brief MS2 input file adapter.

      For the format description take a look at:
      Rapid Commun Mass Spectrom. 2004;18(18):2162-8.

      MS1, MS2, and SQT-three unified, compact, and easily parsed file formats for the
      storage of shotgun proteomic spectra and identifications.

      McDonald WH, Tabb DL, Sadygov RG, MacCoss MJ, Venable J, Graumann J, Johnson JR,
      Cociorva D, Yates JR 3rd.

      PMID: 15317041

  @ingroup FileIO
  */
  class OPENMS_DLLAPI MS2File :
    public ProgressLogger
  {
public:

    /// constructor
    MS2File();

    /// constructor
    virtual ~MS2File();

    template <typename MapType>
    void load(const String & filename, MapType & exp)
    {
      //startProgress(0,0,"loading DTA2D file");

      if (!File::exists(filename))
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
      }
      if (!File::readable(filename))
      {
        throw Exception::FileNotReadable(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
      }

      exp.reset();

      //set DocumentIdentifier
      exp.setLoadedFileType(filename);
      exp.setLoadedFilePath(filename);

      std::ifstream in(filename.c_str());

      UInt spectrum_number = 0;
      typename MapType::SpectrumType spec;
      typename MapType::SpectrumType::PeakType p;

      String line;
      bool first_spec(true);

      // line number counter
      Size line_number = 0;

      while (getline(in, line, '\n'))
      {
        ++line_number;

        line.trim();
        if (line.empty()) continue;

        // header
        if (line[0] == 'H')
        {
          continue;
        }

        // scan
        if (line[0] == 'S')
        {
          if (!first_spec)
          {
            spec.setMSLevel(2);
            spec.setNativeID(String("index=") + (spectrum_number++));
            exp.addSpectrum(spec);
          }
          else
          {
            first_spec = false;
          }
          spec.clear(true);
          line.simplify();
          std::vector<String> split;
          line.split(' ', split);
          if (split.size() != 4)
          {
            throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "line (" + String(line_number) + ") '" + line  + "' should contain four values, got " + String(split.size()) + "!", "");
          }
          spec.getPrecursors().resize(1);
          spec.getPrecursors()[0].setMZ(split[3].toDouble());
          continue;
        }

        // charge-independent analysis
        if (line[0] == 'I')
        {
          continue;
        }

        // charge specification
        if (line[0] == 'Z')
        {
          continue;
        }

        // charge-dependent analysis
        if (line[0] == 'D')
        {
          continue;
        }

        // yet another peak, hopefully
        line.simplify();
        std::vector<String> split;
        line.split(' ', split);
        if (split.size() != 2)
        {
          throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "line (" + String(line_number) + ") '" + line  + "' should contain two values, got " + String(split.size()) + "!", "");
        }

        try
        {
          p.setPosition(split[0].toDouble());
          p.setIntensity(split[1].toFloat());
        }
        catch (Exception::ConversionError /*&e*/)
        {
          throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "ConversionError: line (" + String(line_number) + ") '" + line  + "' does not contain two numbers!", "");
        }
        spec.push_back(p);
      }

      if (!first_spec)
      {
        spec.setMSLevel(2);
        spec.setNativeID(String("index=") + (spectrum_number++));
        exp.addSpectrum(spec);
      }
    }

    /*
    template <typename MapType> void store(const String& filename, MapType& map)
    {

    }
    */

protected:

  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_MS2FILE_H
