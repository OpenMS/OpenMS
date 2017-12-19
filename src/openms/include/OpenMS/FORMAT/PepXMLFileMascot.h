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
// $Authors: Nico Pfeifer $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_PEPXMLFILEMASCOT_H
#define OPENMS_FORMAT_PEPXMLFILEMASCOT_H

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

#endif // OPENMS_FORMAT_PEPXMLFILEMASCOT_H
