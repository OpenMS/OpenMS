// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Immanuel Luhn $
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
  class OPENMS_DLLAPI QuantitativeExperimentalDesign
    : public DefaultParamHandler
  {
    public:

        ///Default constructor
        QuantitativeExperimentalDesign();

        //Destructor
        virtual ~QuantitativeExperimentalDesign();

        /**
            @brief xxxxxxx

            @note xxxxxx
        */
        void applyDesign2Resolver(ProteinResolver& resolver, TextFile& file, StringList& fileNames);

        /**
            @brief xxxxxxx

            @note xxxxxx
        */
        void applyDesign2Quantifier(PeptideAndProteinQuant& quantifier, TextFile& file, StringList& fileNames);
                                    //std::vector< std::pair<PeptideAndProteinQuant::PeptideData,PeptideAndProteinQuant::ProteinQuant> >& result);

    private:
      ///Not implemented
      QuantitativeExperimentalDesign(const QuantitativeExperimentalDesign&);

      ///Not implemented
      QuantitativeExperimentalDesign& operator = (const QuantitativeExperimentalDesign&);

      void mergeFeatureMaps_(FeatureMap<>& map, const String& experiment, StringList& file_paths);

      void mergeConsensusMaps_(ConsensusMap& map, const String& experiment, StringList& file_paths);

      void mergeIDFiles_(std::vector<ProteinIdentification>& proteins, std::vector<PeptideIdentification>& peptides, const String& experiment, StringList& file_paths);

      void findRelevantFilePaths_( std::map<String, StringList>& design2FileBaseName, std::map<String, StringList>& design2FilePath, StringList& filePaths);

      /// Analyze header, which column represents filename and experimental setting
      void analyzeHeader_(UInt& expCol, UInt& fileCol, StringList& header);

      /// Helper class to get separator of the experimental design text file
      void getSeparator_(String& separator);

      /// helper class. which file name belongs to which experimental setting
      void mapFiles2Design_(std::map<String, StringList>& experiments, TextFile& file);


  };

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_QUANTITATION_QUANTITATIVEEXPERIMENTALDESIGN_H
