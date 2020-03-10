// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Svetlana Kutuzova, Douglas McCloskey $
// $Authors: Svetlana Kutuzova, Douglas McCloskey $
// --------------------------------------------------------------------------
 
#pragma once

#include <OpenMS/FORMAT/CsvFile.h>
#include <map>

namespace OpenMS
{
class MzMLFile;
/*
    @brief Scheduler for FIA-MS data batches. Works with FIAMSDataProcessor.

    The class is initialised with the path to the csv file that must contain the following columns:
    - *filename* - the mzML filename for the sample. Must not contains the extension and path. The results filenames follow the pattern @filename_@n_seconds.
    - *dir_input* - the location of the directory with the mzML file
    - *dir_output* - the location of the directory where the results will be stored
    - *resolution* - the resolution of the instrument
    - *polarity* - the charge of the instrument, accepts "positive" or "negative"
    - *db:mapping* - database input file containing three tab-separated columns of mass, formula, identifier for the accurate mass search
    - *db:struct* - database input file containing four tab-separated columns of identifier, name, SMILES, INCHI for the accurate mass search
    - *positive_adducts* - file containing the list of potential positive adducts for the accurate mass search
    - *negative_adducts* - file containing the list of potential negative adducts for the accurate mass search
    - *time* - ";"-separated numbers of seconds to process, f.e. "30;60;90;180"
*/
class OPENMS_DLLAPI FIAMSScheduler 
  {
public:
    FIAMSScheduler() = default;

    /// Default constructor
    FIAMSScheduler(
      String filename,
      String base_dir = "/",
      bool load_cached_ = true
    );

    /// Default desctructor
    ~FIAMSScheduler() = default;

    /// Copy constructor
    FIAMSScheduler(const FIAMSScheduler& cp) = default;

    /// Assignment
    FIAMSScheduler& operator=(const FIAMSScheduler& fdp) = default;

    /**
      @brief Run the FIA-MS data analysis for the batch defined in the @filename_
    */
    void run();

    /**
      @brief Get the batch
    */
    const std::vector<std::map<String, String>>& getSamples();

    /**
      @brief Get the base directory for the relevant paths from the csv file
    */
    const String& getBaseDir();

private:
    /**
      @brief Load the batch from the csv file and store as the vector of maps
    */
    void loadSamples_();

    String filename_;
    String base_dir_;
    bool load_cached_;
    std::vector<std::map<String, String>> samples_;
  };
} // namespace OpenMS
