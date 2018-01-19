// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_PROTXMLFILE_H
#define OPENMS_FORMAT_PROTXMLFILE_H

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

        A documented schema for this format comes with the TPP.

        OpenMS can only read parts of the protein_summary subtree to extract protein-peptide associations. All other parts are silently ignored.

        For protein groups, only the "group leader" (which is annotated with a probability and coverage)
        receives these attributes. All indistinguishable proteins of the same group only have
        an accession and score of -1.

    @note This format will eventually be replaced by the HUPO-PSI (mzIdentML and mzQuantML) AnalysisXML formats!

    @todo Document which metavalues of Protein/PeptideHit are filled when reading ProtXML (Chris)

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

#endif // OPENMS_FORMAT_PROTXMLFILE_H
