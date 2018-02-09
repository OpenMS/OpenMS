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
// $Maintainer: Timo Sachsenberg $
// $Authors: Immanuel Luhn $
// --------------------------------------------------------------------------
#ifndef OPENMS_ANALYSIS_ID_IDRIPPER_H
#define OPENMS_ANALYSIS_ID_IDRIPPER_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>


namespace OpenMS
{
  /**
    @brief Ripping protein/peptide identification according their file origin.

    Helper class, which is used by @ref TOPP_ProteinQuantifier. See there for further documentation.

    @htmlinclude OpenMS_IDRipper.parameters

    @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI IDRipper :
    public DefaultParamHandler
  {
public:

    /// Default constructor
    IDRipper();

    /// Destructor
    ~IDRipper() override;

    /**
      @brief Ripping protein/peptide identification according their file origin

      Iteration over all @p PeptideIdentification. For each annotated file origin create a map entry and store the
      respective @p PeptideIdentification and @p ProteinIdentification.

      @param ripped Contains the protein identification and peptide identification for each file origin annotated in proteins and peptides
      @param proteins Peptide identification annotated with file origin
      @param peptides Protein identification
    */
    void rip(std::map<String, std::pair<std::vector<ProteinIdentification>, std::vector<PeptideIdentification> > > & ripped, std::vector<ProteinIdentification> & proteins, std::vector<PeptideIdentification> & peptides);

private:

    //Not implemented
    /// Copy constructor
    IDRipper(const IDRipper & rhs);

    // Not implemented
    /// Assignment
    IDRipper & operator=(const IDRipper & rhs);

    /// helper function, extracts all protein hits that match the protein accession
    void getProteinHits_(std::vector<ProteinHit> & result, const std::vector<ProteinHit> & protein_hits, const std::vector<String> & protein_accessions);
    /// helper function, returns the string representation of the peptide hit accession
    void getProteinAccessions_(std::vector<String> & result, const std::vector<PeptideHit> & peptide_hits);
    /// helper function, returns the protein identification for the given peptide identification based on the same identifier
    void getProteinIdentification_(ProteinIdentification & result, PeptideIdentification pep_ident, std::vector<ProteinIdentification> & prot_idents);
  };

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_ID_IDRIPPER_H
