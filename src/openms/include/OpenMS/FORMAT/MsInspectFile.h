// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Exception.h>
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
    void load(const std::string& filename, FeatureMapType& feature_map)
    {
      // load input
      TextFile input(filename);

      // reset map
      FeatureMapType fmap;
      feature_map = fmap;

      bool first_line = true;
      for (TextFile::ConstIterator it = input.begin(); it != input.end(); ++it)
      {
        std::string line = *it;

        //ignore comment lines
        if (line.empty() || line[0] == '#') continue;

        //skip leader line
        if (first_line)
        {
          first_line = false;
          continue;
        }

        //split lines: scan\ttime\tmz\taccurateMZ\tmass\tintensity\tcharge\tchargeStates\tkl\tbackground\tmedian\tpeaks\tscanFirst\tscanLast\tscanCount\ttotalIntensity\tsumSquaresDist\tdescription
        std::vector<std::string> parts;
        StringUtils::split(line, '\t', parts);

        if (parts.size() < 18)
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "",StringUtils::toStr("Failed to convert line ")  + StringUtils::toStr((it - input.begin()) + 1) + ". Not enough columns (expected 18 or more, got " + StringUtils::toStr(parts.size()) + ")");
        }

        //create feature
        Feature f;
        Size column_to_convert = 0;
        try
        {
          column_to_convert = 1;
          f.setRT(StringUtils::toDouble(parts[1]));
          column_to_convert = 2;
          f.setMZ(StringUtils::toDouble(parts[2]));
          column_to_convert = 5;
          f.setIntensity(StringUtils::toDouble(parts[5]));
          column_to_convert = 6;
          f.setCharge(StringUtils::toInt32(parts[6]));
          column_to_convert = 8;
          f.setOverallQuality(StringUtils::toDouble(parts[8]));

          column_to_convert = 3;
          f.setMetaValue("accurateMZ", parts[3]);
          column_to_convert = 4;
          f.setMetaValue("mass", StringUtils::toDouble(parts[4]));
          column_to_convert = 7;
          f.setMetaValue("chargeStates", StringUtils::toInt32(parts[7]));
          column_to_convert = 9;
          f.setMetaValue("background", StringUtils::toDouble(parts[9]));
          column_to_convert = 10;
          f.setMetaValue("median", StringUtils::toDouble(parts[10]));
          column_to_convert = 11;
          f.setMetaValue("peaks", StringUtils::toInt32(parts[11]));
          column_to_convert = 12;
          f.setMetaValue("scanFirst", StringUtils::toInt32(parts[12]));
          column_to_convert = 13;
          f.setMetaValue("scanLast", StringUtils::toInt32(parts[13]));
          column_to_convert = 14;
          f.setMetaValue("scanCount", StringUtils::toInt32(parts[14]));
          column_to_convert = 15;
          f.setMetaValue("totalIntensity", StringUtils::toDouble(parts[15]));
          column_to_convert = 16;
          f.setMetaValue("sumSquaresDist", StringUtils::toDouble(parts[16]));
        }
        catch ( Exception::BaseException& )
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "",StringUtils::toStr("Failed to convert value in column ") + StringUtils::toStr(column_to_convert + 1) + " into a number (line '" + StringUtils::toStr((it - input.begin()) + 1) + ")");
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
    void store(const std::string& filename, const SpectrumType& spectrum) const
    {
      std::cerr << "Store() for MsInspectFile not implemented. Filename was: " << filename << ", spec of size " << spectrum.size() << "\n";
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

  };
} // namespace OpenMS

