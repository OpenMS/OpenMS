// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

