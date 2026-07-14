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

    For protein groups, only the "group leader" carries its own probability and
    coverage in the protXML. Indistinguishable siblings (\<indistinguishable_protein\>)
    carry only a protein_name in the file (no probability or coverage of their own);
    on read they inherit the group leader's score (so that score-based filtering does
    not tear groups apart) and are left without a coverage value.

    @note ProteinProphet assigns probability=0 to "unneeded" (subsumable) proteins
    whose peptides are fully explained by higher-probability proteins in the same
    protein_group. These are distinct from truly indistinguishable proteins (which
    share 100% of their peptides and are represented as \<indistinguishable_protein\>
    elements in the protXML). During parsing, probability=0 \<protein\> elements are
    filtered out entirely: no ProteinHit, group entry, or peptide evidence is created
    for them. This prevents downstream tools like ProteinQuantifier from treating
    subsumable proteins as separate indistinguishable groups, which would incorrectly
    cause shared peptides to be discarded during quantification.

    Data filled when reading a protXML:
    - ProteinHit: accession (from @p protein_name), score (ProteinProphet
      probability; indistinguishable siblings inherit the group leader's score),
      and coverage (from @p percent_coverage, group leader only). No meta values
      are set on ProteinHit.
    - PeptideHit: sequence (from @p peptide_sequence, including any modifications
      derived from \<mod_aminoacid_mass\>), score (from @p nsp_adjusted_probability),
      charge, and one PeptideEvidence per protein of the enclosing indistinguishable
      group. Meta values: @p is_unique (1/0 from @p is_nondegenerate_evidence) and
      @p is_contributing (1/0 from @p is_contributing_evidence).
    - The score type of both the ProteinIdentification and PeptideIdentification is
      set to "ProteinProphet probability" (higher is better).

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
    void load(const std::string & filename, ProteinIdentification & protein_ids, PeptideIdentification & peptide_ids);

    /**
        @brief [not implemented yet!] Stores the data in an ProtXML file

        [not implemented yet!]
        The data is stored in the file 'filename'.

        @exception Exception::UnableToCreateFile is thrown if the file could not be created
    */
    void store(const std::string & filename, const ProteinIdentification & protein_ids, const PeptideIdentification & peptide_ids, const std::string & document_id = "");

protected:

    /// reset members after reading/writing
    void resetMembers_();

    /// Docu in base class
    void endElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname) override;

    /// Docu in base class
    void startElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname, const xercesc::Attributes & attributes) override;

    /// Creates a new protein entry (if not already present) and appends it to the current group
    void registerProtein_(const std::string & protein_name);

    /**
        @brief find modification name given a modified AA mass

        Matches a mass of a modified AA to a mod in our modification db
        For ambiguous mods, the first (arbitrary) is returned
        If no mod is found an error is issued and the return string is empty
        @note A duplicate of this function is also used in PepXMLFile

        @param[in,out] mass Modified AA's mass
        @param[in] origin AA one letter code
        @param[in] modification_description [out] Name of the modification, e.g. 'Carboxymethyl (C)'
    */
    void matchModification_(const double mass, const std::string & origin, std::string & modification_description);

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
    /// true when inside a \<protein\> element with probability=0 (subsumable/"unneeded" protein);
    /// all child elements (peptides, indistinguishable_protein) are skipped during parsing
    bool skip_protein_;

    //@}
  };

} // namespace OpenMS

