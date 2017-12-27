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
// $Authors: Immanuel Luhn$
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_QUANTITATION_QUANTITATIVEEXPERIMENTALDESIGN_H
#define OPENMS_ANALYSIS_QUANTITATION_QUANTITATIVEEXPERIMENTALDESIGN_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/ANALYSIS/QUANTITATION/PeptideAndProteinQuant.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ProteinResolver.h>

//#include <vector>

namespace OpenMS
{
  /**
    @brief Merge files according to experimental design

    This class

        @htmlinclude OpenMS_QuantitativeExperimentalDesign.parameters

        @ingroup Analysis_QUANTITATION
  */
  class OPENMS_DLLAPI QuantitativeExperimentalDesign :
    public DefaultParamHandler
  {
public:

    ///Default constructor
    QuantitativeExperimentalDesign();

    //Destructor
    ~QuantitativeExperimentalDesign() override;

    /**
        @brief xxxxxxx

        @note xxxxxx
    */
    void applyDesign2Resolver(ProteinResolver & resolver, TextFile & file, StringList & fileNames);

    /**
        @brief xxxxxxx

        @note xxxxxx
    */
    void applyDesign2Quantifier(PeptideAndProteinQuant & quantifier, TextFile & file, StringList & fileNames);
    //std::vector< std::pair<PeptideAndProteinQuant::PeptideData,PeptideAndProteinQuant::ProteinQuant> >& result);

private:
    ///Not implemented
    QuantitativeExperimentalDesign(const QuantitativeExperimentalDesign &);

    ///Not implemented
    QuantitativeExperimentalDesign & operator=(const QuantitativeExperimentalDesign &);

    void mergeFeatureMaps_(FeatureMap & map, const String & experiment, StringList & file_paths);

    void mergeConsensusMaps_(ConsensusMap & map, const String & experiment, StringList & file_paths);

    void mergeIDFiles_(std::vector<ProteinIdentification> & proteins, std::vector<PeptideIdentification> & peptides, const String & experiment, StringList & file_paths);

    void findRelevantFilePaths_(std::map<String, StringList> & design2FileBaseName, std::map<String, StringList> & design2FilePath, StringList & filePaths);

    /// Analyze header, which column represents filename and experimental setting
    void analyzeHeader_(UInt & expCol, UInt & fileCol, StringList & header);

    /// Helper class to get separator of the experimental design text file
    void getSeparator_(String & separator);

    /// helper class. which file name belongs to which experimental setting
    void mapFiles2Design_(std::map<String, StringList> & experiments, TextFile & file);


  };

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_QUANTITATION_QUANTITATIVEEXPERIMENTALDESIGN_H
