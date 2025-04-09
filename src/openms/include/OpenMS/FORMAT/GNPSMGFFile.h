// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Dorrestein Lab - University of California San Diego - https://dorresteinlab.ucsd.edu/$
// $Authors: Abinesh Sarvepalli and Louis Felix Nothias$
// $Contributors: Fabian Aicheler and Oliver Alka from Oliver Kohlbacher's group at Tubingen University$

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/BinnedSpectrum.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

namespace OpenMS
{
  class OPENMS_DLLAPI GNPSMGFFile : 
    public DefaultParamHandler,
    public ProgressLogger
  {
    public:
      // default c'tor
      GNPSMGFFile();

      // see GNPSExport tool documentation
      /**
      * @brief Create file for GNPS molecular networking.
      * @param consensus_file_path path to consensusXML with spectrum references
      * @param mzml_file_paths path to mzML files referenced in consensusXML. Used to extract spectra as MGF.
      * @param out MGF file with MS2 peak data for molecular networking.
      */
      void store(const String& consensus_file_path, const StringList& mzml_file_paths, const String& out) const;

    private:
      static constexpr double DEF_COSINE_SIMILARITY = 0.9;
      static constexpr double DEF_MERGE_BIN_SIZE = static_cast<double>(BinnedSpectrum::DEFAULT_BIN_WIDTH_HIRES);

//      static constexpr double DEF_PREC_MASS_TOL = 0.5;
//      static constexpr bool DEF_PREC_MASS_TOL_ISPPM = false;

      static constexpr int DEF_PEPT_CUTOFF = 5;
      static constexpr int DEF_MSMAP_CACHE = 50;
  };
}
