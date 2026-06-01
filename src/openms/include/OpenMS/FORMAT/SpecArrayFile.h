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
    void load(const String& filename, FeatureMapType& feature_map)
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
        String line = *it;

        std::vector<String> parts;
        line.split('\t', parts);

        if (parts.size() < 5)
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "", String("Failed to convert line")  + String((it - input.begin()) + 1) + "not enough columns (expected 5 or more, got " + String(parts.size()) + ")");
        }

        Feature f;
        try
        {
          f.setMZ(parts[0].toDouble());
          f.setRT(parts[1].toDouble() * 60.0);
          f.setMetaValue("s/n", parts[2].toDouble());
          f.setCharge(parts[3].toInt());
          f.setIntensity(parts[4].toDouble());
        }
        catch ( Exception::BaseException& )
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "", String("Failed to convert value into a number (line '") + String((it - input.begin()) + 1) + ")");
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
    void store(const String& filename, const SpectrumType& spectrum) const
    {
      std::cerr << "Store() for SpecArrayFile not implemented. Filename was: " << filename << ", spec of size " << spectrum.size() << "\n";
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

  };
} // namespace OpenMS

