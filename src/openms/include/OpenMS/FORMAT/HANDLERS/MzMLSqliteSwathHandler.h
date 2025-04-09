// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>

// forward declarations
struct sqlite3;

namespace OpenMS
{

  namespace Internal
  {

    /**
        @brief Sqlite handler for SWATH data sets

        This class represents a single sqMass file acquired in SWATH / DIA mode
        and provides some useful access to the indices of the individual SWATH
        windows.

    */
    class OPENMS_DLLAPI MzMLSqliteSwathHandler
    {

public:

      /**
          @brief Constructor

          @param filename The sqMass filename
      */
      MzMLSqliteSwathHandler(const String& filename) :
        filename_(filename)
      {}

      /**
          @brief Read SWATH windows boundaries from file

          @return A vector populated with SwathMap, with the following attributes initialized: center, lower and upper

      */
      std::vector<OpenSwath::SwathMap> readSwathWindows();

      /**
          @brief Read indices of MS1 spectra from file

          @return A list of spectral indices for the MS1 spectra
      */
      std::vector<int> readMS1Spectra();

      /**
          @brief Read indices of spectra belonging to specified SWATH window from file

          @param swath_map Contains the upper/lower boundaries of the SWATH window
          @return A list of spectral indices for the provided SWATH window

      */
      std::vector<int> readSpectraForWindow(const OpenSwath::SwathMap & swath_map);

protected:

      String filename_;

      /*
       * These are spectra and chromatogram ids that are global for a specific
       * database file. Keeping track of them allows us to append spectra and
       * chromatograms multiple times to a database.
      */
      Int spec_id_;
      Int chrom_id_;

    };


  }   // namespace Internal
} // namespace OpenMS

