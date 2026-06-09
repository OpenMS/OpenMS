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
    @brief File adapter for SpecArray (.pepList) files.

    The first line is the header and contains the column names:<br>
     m/z	     rt(min)	       snr	      charge	   intensity

    Every subsequent line is a feature.
    Entries are separated by Tab (\\t).


    @ingroup FileIO
  */
  class OPENMS_DLLAPI SpecArrayFile
  {
public:
    /// Default constructor
    SpecArrayFile();
    /// Destructor
    virtual ~SpecArrayFile();

    /**
      @brief Loads a SpecArray file into a featureXML.

      The content of the file is stored in @p features.

      @exception Exception::FileNotFound is thrown if the file could not be opened
      @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    template <typename FeatureMapType>
    void load(const std::string& filename, FeatureMapType& feature_map)
    {
      // load input
      TextFile input(filename, false);

      // reset map
      FeatureMapType fmap;
      feature_map = fmap;

      TextFile::ConstIterator it = input.begin();
      if (it == input.end()) return; // no data to load

      // skip header line
      ++it;
      // process content
      for (; it != input.end(); ++it)
      {
        std::string line = *it;

        std::vector<std::string> parts;
        StringUtils::split(line, '\t', parts);

        if (parts.size() < 5)
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "",std::string("Failed to convert line")  + StringUtils::toStr((it - input.begin()) + 1) + "not enough columns (expected 5 or more, got " + StringUtils::toStr(parts.size()) + ")");
        }

        Feature f;
        try
        {
          f.setMZ(StringUtils::toDouble(parts[0]));
          f.setRT(StringUtils::toDouble(parts[1]) * 60.0);
          f.setMetaValue("s/n", StringUtils::toDouble(parts[2]));
          f.setCharge(StringUtils::toInt32(parts[3]));
          f.setIntensity(StringUtils::toDouble(parts[4]));
        }
        catch ( Exception::BaseException& )
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "",std::string("Failed to convert value into a number (line '") + StringUtils::toStr((it - input.begin()) + 1) + ")");
        }
        feature_map.push_back(f);
      }
    }

    /**
      @brief Stores a featureXML as a SpecArray file.

      NOT IMPLEMENTED

              @exception Exception::UnableToCreateFile is thrown if the file could not be created
    */
    template <typename SpectrumType>
    void store(const std::string& filename, const SpectrumType& spectrum) const
    {
      std::cerr << "Store() for SpecArrayFile not implemented. Filename was: " << filename << ", spec of size " << spectrum.size() << "\n";
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

  };
} // namespace OpenMS

