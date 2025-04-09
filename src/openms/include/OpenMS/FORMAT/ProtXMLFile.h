// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/XMLFile.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief Used to load (storing not supported, yet) ProtXML files

    This class is used to load (storing not supported, yet) documents that implement
    the schema of ProtXML files.

    A documented schema for this format comes with the TPP and can also be
    found at https://github.com/OpenMS/OpenMS/tree/develop/share/OpenMS/SCHEMAS

    OpenMS can only read parts of the protein_summary subtree to extract
    protein-peptide associations. All other parts are silently ignored.

    For protein groups, only the "group leader" (which is annotated with a
    probability and coverage) receives these attributes. All indistinguishable
    proteins of the same group only have an accession and score of -1.

    @todo Document which metavalues of Protein/PeptideHit are filled when reading ProtXML (Chris)
    @todo Writing of protXML is currently not supported

    @ingroup FileIO
  */
  class OPENMS_DLLAPI ProtXMLFile :
    protected Internal::XMLHandler,
    public Internal::XMLFile
  {
public:

    /// A protein group (set of indices into ProteinIdentification)
    typedef ProteinIdentification::ProteinGroup ProteinGroup;

    /// Constructor
    ProtXMLFile();

    /**
        @brief Loads the identifications of an ProtXML file without identifier

        The information is read in and the information is stored in the
        corresponding variables

        @exception Exception::FileNotFound is thrown if the file could not be opened
        @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void load(const String & filename, ProteinIdentification & protein_ids, PeptideIdentification & peptide_ids);

    /**
        @brief [not implemented yet!] Stores the data in an ProtXML file

        [not implemented yet!]
        The data is stored in the file 'filename'.

        @exception Exception::UnableToCreateFile is thrown if the file could not be created
    */
    void store(const String & filename, const ProteinIdentification & protein_ids, const PeptideIdentification & peptide_ids, const String & document_id = "");

protected:

    /// reset members after reading/writing
    void resetMembers_();

    /// Docu in base class
    void endElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname) override;

    /// Docu in base class
    void startElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname, const xercesc::Attributes & attributes) override;

    /// Creates a new protein entry (if not already present) and appends it to the current group
    void registerProtein_(const String & protein_name);

    /**
        @brief find modification name given a modified AA mass

        Matches a mass of a modified AA to a mod in our modification db
        For ambiguous mods, the first (arbitrary) is returned
        If no mod is found an error is issued and the return string is empty
        @note A duplicate of this function is also used in PepXMLFile

        @param mass Modified AA's mass
        @param origin AA one letter code
        @param modification_description [out] Name of the modification, e.g. 'Carboxymethyl (C)'
    */
    void matchModification_(const double mass, const String & origin, String & modification_description);

    /// @name members for loading data
    //@{
    /// Pointer to protein identification
    ProteinIdentification * prot_id_;
    /// Pointer to peptide identification
    PeptideIdentification * pep_id_;
    /// Temporary peptide hit
    PeptideHit * pep_hit_;
    /// protein group
    ProteinGroup protein_group_;


    //@}
  };

} // namespace OpenMS

