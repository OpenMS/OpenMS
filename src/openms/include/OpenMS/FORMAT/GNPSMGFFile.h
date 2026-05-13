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
  /**
    @brief Writer for the consensus-spectrum MGF file consumed by GNPS Feature-Based / Ion-Identity Molecular Networking.

    For each consensus feature in a @ref ConsensusMap, the GNPS workflow expects one
    representative MS2 spectrum. This class reads the spectrum references stored on
    each consensus element, resolves them to the original mzML files, optionally merges
    multiple MS2 scans per feature into a single binned consensus spectrum, and writes
    the result as an MGF file ready to be uploaded to GNPS.

    Configurable behaviour is exposed through @ref DefaultParamHandler — see the defaults
    installed by the constructor for the supported keys (@c output_type — @c "most_intense"
    or @c "merged_spectra"; @c peptide_cutoff — how many top-intensity peptides to consider
    per consensus element; @c ms2_bin_size — bin width when merging fragment ions;
    @c merged_spectra:cos_similarity — cosine threshold for merging).

    Used by the @ref TOPP_GNPSExport tool. See its page for the recommended OpenMS
    metabolomics workflow producing the consensusXML input.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI GNPSMGFFile :
    public DefaultParamHandler,
    public ProgressLogger
  {
    public:
      /// Default constructor; installs the export parameters (see class docs)
      GNPSMGFFile();

      /**
        @brief Write the consensus-spectrum MGF for GNPS molecular networking.

        The consensusXML at @p consensus_file_path must reference MS2 spectra in
        @p mzml_file_paths (typically the same mzMLs that fed the feature finder).
        Per consensus feature, the referenced MS2 scans are optionally merged into a
        single representative spectrum (per the cosine-similarity / bin-width / peak-cutoff
        parameters) and appended to the output MGF.

        @param[in] consensus_file_path Path to a consensusXML file with spectrum references
                                       on each consensus element.
        @param[in] mzml_file_paths     Paths to the mzML files referenced by the consensusXML.
        @param[in] out                 Output MGF path; one spectrum block per consensus feature.
      */
      void store(const String& consensus_file_path, const StringList& mzml_file_paths, const String& out) const;

    private:
      /// Default cosine-similarity threshold used when merging multiple MS2 scans into one representative spectrum
      static constexpr double DEF_COSINE_SIMILARITY = 0.9;
      /// Default m/z bin width for the binned-spectrum representation used during merging (high-res default)
      static constexpr double DEF_MERGE_BIN_SIZE = static_cast<double>(BinnedSpectrum::DEFAULT_BIN_WIDTH_HIRES);

//      static constexpr double DEF_PREC_MASS_TOL = 0.5;
//      static constexpr bool DEF_PREC_MASS_TOL_ISPPM = false;

      /// Default value of the @c peptide_cutoff parameter — how many top-intensity peptides per consensus element to consider when @c output_type is @c merged_spectra (-1 = all)
      static constexpr int DEF_PEPT_CUTOFF = 5;
      /// Declared default but currently not wired to any parameter — kept for source compatibility with earlier versions of the writer
      static constexpr int DEF_MSMAP_CACHE = 50;
  };
}
