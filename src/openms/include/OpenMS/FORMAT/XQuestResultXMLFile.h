// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
// $Maintainer: Lukas Zimmermann $
// $Authors: Lukas Zimmermann $
// --------------------------------------------------------------------------
#ifndef OPENMS_FORMAT_XQUESTRESULTXMLFILE_H
#define OPENMS_FORMAT_XQUESTRESULTXMLFILE_H

#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

namespace OpenMS
{

  class OPENMS_DLLAPI XQuestResultXMLFile :
    public Internal::XMLFile
  {
public:
    XQuestResultXMLFile();
    ~XQuestResultXMLFile();

    /**
     * @brief Load the content of the xquest.xml file into the provided data structures.
     * @param filename Filename of the file which is to be loaded.
     * @param csms Where the spectra with identifications of the input file will be loaded to.
     * @param prot_ids Where the protein identification of the input file will be loaded to.
     * @param min_n_hits_per_spectrum How many minimum hits a spectrum must contain to be loaded to @p csms.
     * @param load_to_peptideHit Whether the data will be loaded as meta values also into the peptide hits, instead just into the PeptideIdentification
     */
    void load(const String & filename,
              std::vector< std::vector< PeptideIdentification > > & csms,
              std::vector< ProteinIdentification > & prot_ids,
              Size min_n_hits_per_spectrum = 0,
              bool load_to_peptideHit = false);

    // Currently not implemented
    //void store(const String &, std::vector< std::vector< PeptideIdentification > > & );

    /**
     * @brief Returns the total number of hits in the file
     * @return total number of hits in the file
     */
    int getNumberOfHits() const;

    /**
     * @brief Returns minimum score among the hits in the file.
     * @return Minimum score among the hits in the file.
     */
    double getMinScore() const;

    /**
     * @brief Returns maximum score among the hits in the file.
     * @return Maximum score among the hits in the file.
     */
    double getMaxScore() const;


private:
    int n_hits_; // Total number of hits within the result file

    double min_score_; // Minimum score encountered in file
    double max_score_; // Maximum score encountered in file
  };
} // namespace OpenMS
#endif // OPENMS_FORMAT_XQUESTRESULTXMLFILE_H
