// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h> // OPENMS_DLLAPI
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/KERNEL/MSExperiment.h>

namespace OpenMS
{
  /**
    @brief Load Chromeleon HPLC text file and save it into a `MSExperiment`.

    An example of the expected format:
    > Raw Data:
    > Time (min)	Step (s)	Value (mAU)
    > 0.003333	0.200	-0.002496
    > 0.006667	0.200	-0.017589
    > ...
  */
  class OPENMS_DLLAPI ChromeleonFile
  {
public:
    /// Constructor
    ChromeleonFile() = default;
    /// Destructor
    virtual ~ChromeleonFile() = default;

    /**
      @brief Load the file's data and metadata, and save it into a `MSExperiment`.

      @param[in] filename Path to the Chromeleon input file
      @param[out] experiment The variable into which the extracted information will be saved
    */
    void load(const String& filename, MSExperiment& experiment) const;

    /**
      @brief Remove commas from the string (used as thousands separators) and
      parse its value

      @param[in] number A string representing a floating-point number
      @return The value converted to `double`
    */
    double removeCommasAndParseDouble(String& number) const;
  };
}

