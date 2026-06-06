// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#pragma once

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
      Rapid Communications in Mass Spectrometry. 2004;18(18):2162-8.

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
    ~MS2File() override;

    template <typename MapType>
    void load(const std::string & filename, MapType & exp)
    {
      //startProgress(0,0,"loading DTA2D file");

      if (!File::exists(filename))
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
      if (!File::readable(filename))
      {
        throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }

      exp.reset();

      //set DocumentIdentifier
      exp.setLoadedFileType(filename);
      exp.setLoadedFilePath(filename);

      std::ifstream in(filename.c_str());

      UInt spectrum_number = 0;
      typename MapType::SpectrumType spec;
      typename MapType::SpectrumType::PeakType p;

      std::string line;
      bool first_spec(true);

      // line number counter
      Size line_number = 0;

      while (getline(in, line, '\n'))
      {
        ++line_number;

        StringUtils::trim(line);
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
            spec.setNativeID(StringUtils::toStr("index=") + (spectrum_number++));
            exp.addSpectrum(spec);
          }
          else
          {
            first_spec = false;
          }
          spec.clear(true);
          StringUtils::simplify(line);
          std::vector<std::string> split;
          StringUtils::split(line, ' ', split);
          if (split.size() != 4)
          {
            throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "line (" + StringUtils::toStr(line_number) + ") '" + line  + "' should contain four values, got " + StringUtils::toStr(split.size()) + "!", "");
          }
          spec.getPrecursors().resize(1);
          spec.getPrecursors()[0].setMZ(StringUtils::toDouble(split[3]));
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
        StringUtils::simplify(line);
        std::vector<std::string> split;
        StringUtils::split(line, ' ', split);
        if (split.size() != 2)
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "line (" + StringUtils::toStr(line_number) + ") '" + line  + "' should contain two values, got " + StringUtils::toStr(split.size()) + "!", "");
        }

        try
        {
          p.setPosition(StringUtils::toDouble(split[0]));
          p.setIntensity(StringUtils::toFloat(split[1]));
        }
        catch ( Exception::ConversionError& )
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "ConversionError: line (" + StringUtils::toStr(line_number) + ") '" + line  + "' does not contain two numbers!", "");
        }
        spec.push_back(p);
      }

      if (!first_spec)
      {
        spec.setMSLevel(2);
        spec.setNativeID(StringUtils::toStr("index=") + (spectrum_number++));
        exp.addSpectrum(spec);
      }
      exp.updateRanges();
    }

protected:

  };

} // namespace OpenMS

