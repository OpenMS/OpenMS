// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
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
    String source_file = mapping.getPrimaryMSRunPath(pep_id);

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
    String getPrimaryMSRunPath(const PeptideIdentification& pepid) const;

    /// Check if the mapping contains an entry for the given identifier
    bool hasIdentifier(const String& identifier) const;

    /// Get the identifier for a given MS run path list (throws if not found)
    const String& getIdentifier(const StringList& ms_run_paths) const;

    /// Check if the mapping is empty
    bool empty() const;

    /// Get the number of identifier mappings
    Size size() const;

    /// Get the MS run paths for a given identifier (returns empty list if not found)
    const StringList& getMSRunPaths(const String& identifier) const;

    /// Get all identifiers in this mapping
    std::vector<String> getIdentifiers() const;

    /// Check if the mapping contains an entry for the given MS run paths
    bool hasRunPath(const StringList& ms_run_paths) const;

    /// Try to get identifier for a given MS run path list (returns false if not found)
    bool tryGetIdentifier(const StringList& ms_run_paths, String& identifier) const;

  private:
    static const StringList empty_stringlist_; ///< Empty list returned by getMSRunPaths when identifier not found
    std::map<String, StringList> identifier_to_msrunpath_;
    std::map<StringList, String> runpath_to_identifier_;
  };

} // namespace OpenMS
