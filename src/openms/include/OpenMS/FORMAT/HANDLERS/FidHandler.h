// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Guillaume Belz $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
#include <fstream>

namespace OpenMS
{
  namespace Internal
  {
    /**
      @brief Read-only fid File handler for XMass Analysis.

      fid File contains intensity array. Intensity for each point are coded in 4 bytes integer.

      @note Do not use this class directly. It is only needed for XMassFile.
    */
    class OPENMS_DLLAPI FidHandler :
      public std::ifstream
    {
public:
      /**
        @brief Constructor with filename.

        Open fid File as stream and initialize index.

        @param filename to fid File.
      */
      explicit FidHandler(const String & filename);

      /// Destructor
      ~FidHandler() override;

      /// Get index of current position (without position moving).
      Size getIndex() const;

      /// Get intensity of current position and move to next position.
      Size getIntensity();

private:
      /// Private default constructor
      FidHandler();

      /// Index of position
      Size index_;
    };

  }   // namespace Internal
} // namespace OpenMS

