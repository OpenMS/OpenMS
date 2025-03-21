// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>

// boost shared_ptr
#include <boost/shared_ptr.hpp>

namespace OpenMS
{

  // forward declaration
  class ConsensusMap;
  class String;
  class IsobaricQuantitationMethod;
  class AASequence;
  
  /**
    @brief Implements the export of consensusmaps into the IBSpectra format
           used by isobar to load quantification results.
  */
  class OPENMS_DLLAPI IBSpectraFile
  {
public:
    /**
      @brief Constructor.
    */
    IBSpectraFile();

    /**
      @brief Copy constructor.
    */
    IBSpectraFile(const IBSpectraFile& other);

    /**
      @brief Assignment operator.
    */
    IBSpectraFile& operator=(const IBSpectraFile& rhs);

    /**
      @brief Writes the contents of the ConsensusMap cm into the file named by filename.

      @param filename The name of the file where the contents of cm should be stored.
      @param cm The ConsensusMap that should be exported to filename.

      @throws Exception::InvalidParameter if the ConsensusMap does not hold the result of an isobaric quantification experiment (e.g., itraq).
    */
    void store(const String& filename, const ConsensusMap& cm);
private:

    /**
      @brief Guesses the type of isobaric quantitation performed on the experiment.

      @throws Exception::InvalidParameter if the ConsensusMap does not hold the result of an isobaric quantification experiment (e.g., itraq).
    */
    boost::shared_ptr<IsobaricQuantitationMethod> guessExperimentType_(const ConsensusMap& cm);

    /**
      @brief Constructs the matching file header for the given quantitation method.
      @param quantMethod The used quantitation method.
      @return The header of the IBSpectra file for the given quantitation method.
    */
    StringList constructHeader_(const IsobaricQuantitationMethod& quantMethod);
    
    /**
      @brief Generates the modification string for the given AASequence.
     
      @param sequence The sequence for which the modification string should be generated.
      @return The modification string for the given sequence.
    */
    String getModifString_(const AASequence& sequence);
  };

}

