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
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger, Hannes Roest $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_OPENSWATH_TRANSITIONPQPFILE_H
#define OPENMS_ANALYSIS_OPENSWATH_TRANSITIONPQPFILE_H

#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>

#include <boost/range/algorithm.hpp>
#include <boost/range/algorithm_ext/erase.hpp>
#include <sqlite3.h>
#include <iostream>

namespace OpenMS
{

  /**
      @brief This class can convert TraML and PQP files into each other

      The PQP files are SQLite databases consisting of several tables representing the data contained in TraML files.

  @htmlinclude OpenMS_TransitionPQPFile.parameters

  */
  class OPENMS_DLLAPI TransitionPQPFile :
    public TransitionTSVFile
  {

private:
    /** @brief Read PQP SQLite file
     *
     * @param filename The input file
     * @param transition_list The output list of transitions
     * @param legacy_traml_id Should legacy TraML IDs be used (boolean)?
     *
    */
    void readPQPInput_(const char* filename, std::vector<TSVTransition>& transition_list, bool legacy_traml_id = false);

    /** @brief Write a TargetedExperiment to a file
     *
     * @param filename Name of the output file
     * @param targeted_exp The data structure to be written to the file
    */
    void writePQPOutput_(const char* filename, OpenMS::TargetedExperiment& targeted_exp);

public:

    //@{
    /// Constructor
    TransitionPQPFile();

    /// Destructor
    ~TransitionPQPFile() override;
    //@}

    /** @brief Write out a targeted experiment (TraML structure) into a PQP file
     *
     * @param filename The output file
     * @param targeted_exp The targeted experiment
     *
    */
    void convertTargetedExperimentToPQP(const char* filename, OpenMS::TargetedExperiment& targeted_exp);

    /** @brief Read in a PQP file and construct a targeted experiment (TraML structure)
     *
     * @param filename The input file
     * @param targeted_exp The output targeted experiment
     * @param legacy_traml_id Should legacy TraML IDs be used (boolean)?
     *
    */
    void convertPQPToTargetedExperiment(const char* filename, OpenMS::TargetedExperiment& targeted_exp, bool legacy_traml_id = false);

    /** @brief Read in a PQP file and construct a targeted experiment (Light transition structure)
     *
     * @param filename The input file
     * @param targeted_exp The output targeted experiment
     * @param legacy_traml_id Should legacy TraML IDs be used (boolean)?
     *
    */
    void convertPQPToTargetedExperiment(const char* filename, OpenSwath::LightTargetedExperiment& targeted_exp, bool legacy_traml_id = false);

  };
}

#endif // OPENMS_ANALYSIS_OPENSWATH_TRANSITIONPQPREADER_H

