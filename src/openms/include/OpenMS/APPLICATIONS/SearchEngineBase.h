// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>


namespace OpenMS
{
  
  /**
    @brief Base class for TOPP search-engine adapters.

    Sits between @ref TOPPBase and concrete search-engine adapters
    (CometAdapter, MSGFPlusAdapter, ...) and bundles a few conventions
    that every adapter shares:
      - common @c -in (spectra) and @c -database (FASTA) parameters with
        helpers (@ref getRawfileName, @ref getDBFilename) that read
        them and resolve them against the search-engine constraints,
      - optional post-search reindexing via @ref PeptideIndexing with a
        registered @c -reindex switch (@ref registerPeptideIndexingParameter_,
        @ref reindex_).

    @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI SearchEngineBase : public TOPPBase
  {
  public:
    /// Default construction is disabled; a tool name and description are required.
    SearchEngineBase() = delete;

    /// Copy construction is disabled.
    SearchEngineBase(const SearchEngineBase&) = delete;

    /**
      @brief Construct the adapter (signature mirrors @ref TOPPBase to forward arguments verbatim).

      @param[in] name             Tool name.
      @param[in] description      One-line tool description.
      @param[in] official         If @c true, the tool is treated as an
                                  official TOPP tool: the name is
                                  cross-checked against the TOPP-tool
                                  list and a warning is printed if it
                                  is missing.
      @param[in] citations        Citations associated with this TOPP
                                  tool; printed during @c --help.
      @param[in] toolhandler_test Whether to check that the tool is
                                  registered with the @ref ToolHandler.
                                  Disable for unit tests only.
    */
    SearchEngineBase(const String& name, const String& description, bool official = true, const std::vector<Citation>& citations = {}, bool toolhandler_test = true);

    /// Destructor.
    ~SearchEngineBase() override;

    /**
      @brief Resolve the @c -in argument and verify the spectra match a search engine's expectations.

      Reads the @c -in parameter, identifies the file format and -- if
      it is mzML -- inspects the per-MS-level centroid metadata to make
      sure spectra of MS level @p ms_level are centroided. mzML
      profile spectra trigger an @c IllegalArgument exception unless
      the @c -force flag is set (in which case the call only logs a
      warning and proceeds). Format-specific shortcuts: MGF and Bruker
      TDF data are accepted without further inspection; for any other
      format only a general warning is issued.

      @param[in] ms_level MS level to check for their peak type
                          (centroided/profile).
      @return @c -in path (relative or absolute as supplied).

      @throws OpenMS::Exception::FileEmpty       when the mzML file
                                                 has no spectra of the
                                                 requested level.
      @throws OpenMS::Exception::IllegalArgument when the mzML file
                                                 contains profile
                                                 spectra at this level
                                                 (or no centroided
                                                 spectra) and the
                                                 @c -force flag is not
                                                 set.
    */
    String getRawfileName(int ms_level = 2) const;


    /**
      @brief Resolve the search-database path.

      Returns @p db if non-empty, otherwise the value of the
      @c -database parameter. If the resulting path is not directly
      readable, the OpenMS database search paths (see
      @ref OpenMS::File::findDatabase) are scanned.

      @param[in] db Optional explicit database name; used in place of
                    the @c -database parameter when non-empty
                    (occasionally required for special database
                    formats, e.g. OMSSA).
      @return Resolved database path (relative or absolute).

      @throws OpenMS::Exception::FileNotFound when the name cannot be
                                              resolved against the
                                              OpenMS database paths.
    */
    String getDBFilename(const String& db = "") const;


    /**
      @brief Register the @c -reindex switch and a hidden @c PeptideIndexing parameter sub-tree.

      Adds a @c -reindex string option (@c "true"/@c "false", default
      @c "true") plus the @c PeptideIndexing parameter tree under
      a @c PeptideIndexing: prefix. The supplied
      @p peptide_indexing_parameter is patched in two ways before
      being registered: @c missing_decoy_action is forced to @c "warn",
      and a fixed list of advanced options
      (@c decoy_string, @c decoy_string_position,
      @c missing_decoy_action, @c enzyme:name, @c enzyme:specificity,
      @c write_protein_sequence, @c write_protein_description,
      @c keep_unreferenced_proteins, @c unmatched_action, @c aaa_max,
      @c mismatches_max, @c IL_equivalent) are tagged @c "advanced".

      @param[in] peptide_indexing_parameter Defaults for the indexer;
                                            search-engine adapters can
                                            customise these (e.g.
                                            non-tryptic digestion)
                                            before passing them in.
    */
    virtual void registerPeptideIndexingParameter_(Param peptide_indexing_parameter);

    /**
      @brief Optionally re-run @ref PeptideIndexing against the parsed identifications.

      A no-op when the @c -reindex parameter is anything other than
      @c "true". When enabled, runs @ref PeptideIndexing with the
      adapter's @c PeptideIndexing: sub-tree merged onto the indexer's
      defaults, against the database resolved via @ref getDBFilename.

      @param[in,out] protein_identifications  Protein hits; updated in place.
      @param[in,out] peptide_identifications  Peptide hits; updated in place.
      @return @ref EXECUTION_OK on success or when @c -reindex is
              disabled; @ref INPUT_FILE_EMPTY when the indexer reports
              @c DATABASE_EMPTY; @ref UNEXPECTED_RESULT when the
              indexer reports @c UNEXPECTED_RESULT; @ref UNKNOWN_ERROR
              for any other non-OK indexer exit code (excluding
              @c PEPTIDE_IDS_EMPTY, which is treated as success).
    */
    virtual SearchEngineBase::ExitCodes reindex_(
      std::vector<ProteinIdentification>& protein_identifications,
      PeptideIdentificationList& peptide_identifications) const;
  }; // end SearchEngineBase

}   // end NS OpenMS
