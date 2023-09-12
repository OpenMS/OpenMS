// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Immanuel Luhn$
// --------------------------------------------------------------------------

#pragma once

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

private:
    ///Not implemented
    QuantitativeExperimentalDesign(const QuantitativeExperimentalDesign &);

    ///Not implemented
    QuantitativeExperimentalDesign & operator=(const QuantitativeExperimentalDesign &);

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

