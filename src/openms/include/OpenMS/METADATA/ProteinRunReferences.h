// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <string>
#include <vector>

namespace OpenMS
{
  class PeptideIdentification;
  class PeptideIdentificationList;
  class ProteinIdentification;

  /**
    @brief Checks that every peptide identification references an existing protein identification run.

    Every PeptideIdentification needs the ProteinIdentification run (the search run) whose identifier it names
    (PeptideIdentification::getIdentifier() == ProteinIdentification::getIdentifier()). The run may have no protein
    hits and no database (e.g. peptidomics or de novo results). File writers and loaders enforce this contract with
    the checks below, so peptide identifications are never silently dropped or re-assigned.

    @ingroup Metadata
  */
  class OPENMS_DLLAPI ProteinRunReferences
  {
  public:
    /// The error message for a peptide identification whose @p identifier names no protein identification run
    static std::string missingRunMessage(const std::string& identifier);

    /// Throws Exception::InvalidParameter (with missingRunMessage()) if the @p identifier names none of the @p runs
    static void check(const std::vector<ProteinIdentification>& runs, const std::string& identifier);

    /**
      @brief Throws Exception::InvalidParameter if a peptide identification names no protein identification run

      The message (see missingRunMessage()) names the identifier of the first such peptide identification.
    */
    static void check(const std::vector<ProteinIdentification>& runs, const PeptideIdentificationList& peptides);

    /// Same as above, for peptide identifications given by pointers
    static void check(const std::vector<ProteinIdentification>& runs, const std::vector<const PeptideIdentification*>& peptides);

    /**
      @brief Throws Exception::ElementNotFound if a peptide evidence names a protein that is no protein hit of its run

      featureXML and consensusXML reference the proteins of a peptide hit's evidences among the protein hits of the
      peptide identification's run (empty accessions are not written), e.g. for the peptide identifications of a map
      (see IdentificationDataConverter::exportMapIdentifications()). Includes check(): a peptide identification
      without its run throws Exception::InvalidParameter (with missingRunMessage()).
    */
    static void checkProteinAccessions(const std::vector<ProteinIdentification>& runs, const std::vector<const PeptideIdentification*>& peptides);
  };
} // namespace OpenMS
