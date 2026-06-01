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
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/FORMAT/TextFile.h>

#include <fstream>
#include <vector>

namespace OpenMS
{
  /**
      @brief File adapter for Enhanced DTA files.

  Input text file containing tab, space or comma separated columns.
  The separator between columns is checked in the first line in this order.

  It supports three variants of this format.

  - Columns are: RT, MZ, Intensity
    A header is optional.

  - Columns are: RT, MZ, Intensity, Charge, <Meta-Data> columns{0,}
    A header is mandatory.

    Example:
    @code
    RT m/z Intensity charge mymeta1 mymeta2
    321 405.233 24543534 2 lala  lili
    321 406.207 4343344  2 blubb blabb
    @endcode

  - Columns are: (RT, MZ, Intensity, Charge){1,}, <Meta-Data> columns{0,}
    Header is mandatory.
    First quadruplet is the consensus. All following quadruplets describe the sub-features.
    This variant is discerned from variant #2 by the name of the fifth column, which is required to be RT1 (or rt1).
    All other column names for sub-features are faithfully ignored.

    Example:
    @code
    RT MZ INT CHARGE RT1 MZ1 INT1 CHARGE1 RT2 MZ2 INT2 CHARGE2
    321 405 100 2 321 405 100 2 321 406 50 2
    323 406 200 2 323 406 200 2 323 407 100 2 323 407 50 2
    @endcode

  @ingroup FileIO
*/
  class OPENMS_DLLAPI EDTAFile
  {
public:
    /// Default constructor
    EDTAFile();
    /// Destructor
    virtual ~EDTAFile();

private:
    /**
     * Check if column exists and convert String into double.
     */
    double checkedToDouble_(const std::vector<String> & parts, Size index, double def = -1);

    /**
     * Check if column exists and convert String into Int.
     */
    Int checkedToInt_(const std::vector<String> & parts, Size index, Int def = -1);

public:
    /**
      @brief Loads a EDTA file into a consensusXML.

              The content of the file is stored in @p features.

              @exception Exception::FileNotFound is thrown if the file could not be opened
              @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void load(const String & filename, ConsensusMap & consensus_map);

    /**
      @brief Stores a ConsensusMap as an enhanced DTA file.

      NOT IMPLEMENTED

              @exception Exception::UnableToCreateFile is thrown if the file could not be created
    */
    void store(const String & filename, const ConsensusMap & map) const;


    /**
      @brief Stores a FeatureMap as an enhanced DTA file.

      Creating columns: RT, m/z, intensity, charge

              @exception Exception::UnableToCreateFile is thrown if the file could not be created
    */
    void store(const String & filename, const FeatureMap & map) const;
  };
} // namespace OpenMS

