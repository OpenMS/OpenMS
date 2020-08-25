// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/FORMAT/TextFile.h>

#include <fstream>
#include <vector>

namespace OpenMS
{
  /**
      @brief File adapter for MsInspect files.

  Lines with "#" are comments and are ignored.

  The first non-comment line is the header and contains the column names:<br>
  scan	time	mz	accurateMZ	mass	intensity	charge	chargeStates	kl	background	median	peaks	scanFirst	scanLast	scanCount	totalIntensity	sumSquaresDist	description

  Every subsequent line is a feature.

  @ingroup FileIO
*/
  class OPENMS_DLLAPI MsInspectFile
  {
public:
    /// Default constructor
    MsInspectFile();
    /// Destructor
    virtual ~MsInspectFile();

    /**
              @brief Loads a MsInspect file into a featureXML.

              The content of the file is stored in @p features.

              @exception Exception::FileNotFound is thrown if the file could not be opened
              @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    template <typename FeatureMapType>
    void load(const String& filename, FeatureMapType& feature_map)
    {
      // load input
      TextFile input(filename);

      // reset map
      FeatureMapType fmap;
      feature_map = fmap;

      bool first_line = true;
      for (TextFile::ConstIterator it = input.begin(); it != input.end(); ++it)
      {
        String line = *it;

        //ignore comment lines
        if (line.empty() || line[0] == '#') continue;

        //skip leader line
        if (first_line)
        {
          first_line = false;
          continue;
        }

        //split lines: scan\ttime\tmz\taccurateMZ\tmass\tintensity\tcharge\tchargeStates\tkl\tbackground\tmedian\tpeaks\tscanFirst\tscanLast\tscanCount\ttotalIntensity\tsumSquaresDist\tdescription
        std::vector<String> parts;
        line.split('\t', parts);

        if (parts.size() < 18)
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "", String("Failed to convert line ")  + String((it - input.begin()) + 1) + ". Not enough columns (expected 18 or more, got " + String(parts.size()) + ")");
        }

        //create feature
        Feature f;
        Size column_to_convert = 0;
        try
        {
          column_to_convert = 1;
          f.setRT(parts[1].toDouble());
          column_to_convert = 2;
          f.setMZ(parts[2].toDouble());
          column_to_convert = 5;
          f.setIntensity(parts[5].toDouble());
          column_to_convert = 6;
          f.setCharge(parts[6].toInt());
          column_to_convert = 8;
          f.setOverallQuality(parts[8].toDouble());

          column_to_convert = 3;
          f.setMetaValue("accurateMZ", parts[3]);
          column_to_convert = 4;
          f.setMetaValue("mass", parts[4].toDouble());
          column_to_convert = 7;
          f.setMetaValue("chargeStates", parts[7].toInt());
          column_to_convert = 9;
          f.setMetaValue("background", parts[9].toDouble());
          column_to_convert = 10;
          f.setMetaValue("median", parts[10].toDouble());
          column_to_convert = 11;
          f.setMetaValue("peaks", parts[11].toInt());
          column_to_convert = 12;
          f.setMetaValue("scanFirst", parts[12].toInt());
          column_to_convert = 13;
          f.setMetaValue("scanLast", parts[13].toInt());
          column_to_convert = 14;
          f.setMetaValue("scanCount", parts[14].toInt());
          column_to_convert = 15;
          f.setMetaValue("totalIntensity", parts[15].toDouble());
          column_to_convert = 16;
          f.setMetaValue("sumSquaresDist", parts[16].toDouble());
        }
        catch ( Exception::BaseException& )
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "", String("Failed to convert value in column ") + String(column_to_convert + 1) + " into a number (line '" + String((it - input.begin()) + 1) + ")");
        }
        f.setMetaValue("description", parts[17]);
        feature_map.push_back(f);
      }

    }

    /**
      @brief Stores a featureXML as a MsInspect file.

      NOT IMPLEMENTED

              @exception Exception::UnableToCreateFile is thrown if the file could not be created
    */
    template <typename SpectrumType>
    void store(const String& filename, const SpectrumType& spectrum) const
    {
      std::cerr << "Store() for MsInspectFile not implemented. Filename was: " << filename << ", spec of size " << spectrum.size() << "\n";
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

  };
} // namespace OpenMS

