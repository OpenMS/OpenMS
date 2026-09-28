// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/DATASTRUCTURES/TypeAliases.h>
#include <OpenMS/config.h>

#include <map>
#include <vector>

namespace OpenMS
{
  // Forward declarations
  class ProteinIdentification;
  class PeptideIdentification;

  /**
    @brief Two-way mapping from ms-run-path to protID|pepID-identifier.

    This class encapsulates the mapping between protein identification identifiers
    and their associated MS run paths. It is used to resolve the correct source file
    for peptide identifications, especially in merged identification results.

    Usage example:
    @code
    // Create mapping from protein identifications
    IdentifierMSRunMapper mapping(protein_ids);

    // Get the source file for a peptide identification
    std::string source_file = mapping.getPrimaryMSRunPath(pep_id);

    // Build a USI using the mapping
    USI usi = pep_id.buildUSI(mapping, "PXD000561");
    @endcode

    @note The internal data structures are private to avoid exposing nested containers
          that are difficult to wrap in Python bindings.

    @see ProteinIdentification, PeptideIdentification

    @ingroup Metadata
  */
  class OPENMS_DLLAPI IdentifierMSRunMapper
  {
  public:
    IdentifierMSRunMapper() = default;

    /// Construct mapping from a vector of ProteinIdentifications
    explicit IdentifierMSRunMapper(const std::vector<ProteinIdentification>& prot_ids);

    /// Create/update mapping from a vector of ProteinIdentifications
    void create(const std::vector<ProteinIdentification>& prot_ids);

    /// Get the primary MS run path for a PeptideIdentification (using id_merge_index metadata)
    std::string getPrimaryMSRunPath(const PeptideIdentification& pepid) const;

    /**
      @brief Refuse a PSM of a merged run that cannot be attributed to an origin file

      A run with several @c spectra_data entries resolves per PSM only through
      @c id_merge_index. Without a usable one, getPrimaryMSRunPath() falls back to the run's
      @e first file, silently labelling every PSM of the run with it. Callers that use the
      resolved path as a key -- QPX @c run_file_name is a primary-key component of the psm,
      feature and pg views -- must refuse such input rather than publish it, so run this as a
      preflight before any output is opened.

      Runs with 0 or 1 path are exempt: unmerged input has nothing to disambiguate, and
      resolution cannot be wrong.

      @param[in] pep_id The PSM to check
      @param[in] psm_index Only used to make the diagnostic locatable
      @throw Exception::MissingInformation if the run is merged and @c id_merge_index is
             missing, is not an integer, or is out of range
    */
    void validateMergeIndex(const PeptideIdentification& pep_id, size_t psm_index) const;

    /// Check if the mapping contains an entry for the given identifier
    bool hasIdentifier(const std::string& identifier) const;

    /// Get the identifier for a given MS run path list (throws if not found)
    const std::string& getIdentifier(const StringList& ms_run_paths) const;

    /// Check if the mapping is empty
    bool empty() const;

    /// Get the number of identifier mappings
    Size size() const;

    /// Get the MS run paths for a given identifier (returns empty list if not found)
    const StringList& getMSRunPaths(const std::string& identifier) const;

    /// Get all identifiers in this mapping
    std::vector<std::string> getIdentifiers() const;

    /// Check if the mapping contains an entry for the given MS run paths
    bool hasRunPath(const StringList& ms_run_paths) const;

    /// Try to get identifier for a given MS run path list (returns false if not found)
    bool tryGetIdentifier(const StringList& ms_run_paths, std::string& identifier) const;

  private:
    static const StringList empty_stringlist_; ///< Empty list returned by getMSRunPaths when identifier not found
    std::map<std::string, StringList> identifier_to_msrunpath_;
    std::map<StringList, std::string> runpath_to_identifier_;
  };

} // namespace OpenMS
