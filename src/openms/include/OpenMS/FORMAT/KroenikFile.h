// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Exception.h>

#include <OpenMS/KERNEL/FeatureMap.h>

namespace OpenMS
{
  class String;

  /**
    @brief File adapter for Kroenik (HardKloer sibling) files.

    The first line is the header and contains the column names:<br>
    File,  First Scan,  Last Scan,  Num of Scans,  Charge,  Monoisotopic Mass,  Base Isotope Peak,  Best Intensity,  Summed Intensity,  First RTime,  Last RTime,  Best RTime,  Best Correlation,  Modifications

    Every subsequent line is a feature.

    All properties in the file are converted to Feature properties, whereas "First Scan", "Last Scan", "Num of Scans" and "Modifications" are stored as
    metavalues with the following names "FirstScan", "LastScan", "NumOfScans" and "AveragineModifications".

    The width in m/z of the overall convex hull of each feature is set to 3 Th in lack of a value provided by the Kroenik file.

    @note Kroenik files are Tab (\\t) separated files.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI KroenikFile
  {
public:
    /// Default constructor
    KroenikFile();
    /// Destructor
    virtual ~KroenikFile();

    /**
      @brief Loads a Kroenik file into a featureXML.

      The content of the file is stored in @p features.

      @exception Exception::FileNotFound is thrown if the file could not be opened
      @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void load(const String& filename, FeatureMap& feature_map);

    /**
      @brief Stores a featureXML as a Kroenik file.

      NOT IMPLEMENTED

      @exception Exception::UnableToCreateFile is thrown if the file could not be created
    */
    template <typename SpectrumType>
    void store(const String& filename, const SpectrumType& spectrum) const
    {
      std::cerr << "Store() for KroenikFile not implemented. Filename was: " << filename << ", spec of size " << spectrum.size() << "\n";
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

  };
} // namespace OpenMS

