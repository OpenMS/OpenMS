// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Nico Pfeifer $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief Used to load Mascot PepXML files

        A schema for this format can be found at http://www.matrixscience.com/xmlns/schema/pepXML_v18/pepXML_v18.xsd.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI PepXMLFileMascot :
    protected Internal::XMLHandler,
    public Internal::XMLFile
  {
public:

    /// Constructor
    PepXMLFileMascot();

    /**
        @brief Loads peptide sequences with modifications out of a PepXML file

        @exception Exception::FileNotFound is thrown if the file could not be opened
        @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void load(const String & filename, std::map<String, std::vector<AASequence> > & peptides);

protected:

    // Docu in base class
    void endElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname) override;

    // Docu in base class
    void startElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname, const xercesc::Attributes & attributes) override;

    void matchModification_(double mass, String & modification_description);

    /// @name members for loading data
    //@{
    /// Pointer to fill in protein identifications

    /// The title of the actual spectrum
    String actual_title_;

    /// The sequence of the actual peptide hit
    String actual_sequence_;

    /// The modifications of the actual peptide hit (position is 1-based)
    std::vector<std::pair<String, UInt> > actual_modifications_;

    /// The peptides together with the spectrum title
    std::map<String, std::vector<AASequence> > * peptides_;

    /// stores the actual peptide sequences
    std::vector<AASequence> actual_aa_sequences_;

    /// stores the fixed residue modifications
    std::vector<String> fixed_modifications_;

    /// stores the variable residue modifications
    std::vector<std::pair<String, double> > variable_modifications_;
    //@}
  };

} // namespace OpenMS

