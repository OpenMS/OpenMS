// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_PERCOLATOROUTFILE_H
#define OPENMS_FORMAT_PERCOLATOROUTFILE_H

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <boost/regex.hpp>

#include <map>
#include <vector>

namespace OpenMS
{
  /**
     @brief Class for reading Percolator tab-delimited output files.

     For PSM-level output, the file extension should be ".psms".
  */
  class OPENMS_DLLAPI PercolatorOutfile
  {

  public:
    
    /// Types of Percolator scores
    enum ScoreType { QVALUE, POSTERRPROB, SCORE, SIZE_OF_SCORETYPE };

    /// Names of Percolator scores (to match ScoreType)
    static const std::string score_type_names[];

    /// Return a score type given its name
    static enum ScoreType getScoreType(String score_type_name);
    
    /// Constructor
    PercolatorOutfile();

    /// Loads a Percolator output file
    void load(const String& filename, ProteinIdentification& proteins, 
              std::vector<PeptideIdentification>& peptides, 
              enum ScoreType output_score = QVALUE, 
              const String& psm_regex = "", bool count_from_zero = false,
              const MSExperiment<>* experiment_p = 0);

  private:

    /// Information about a single fragment spectrum
    struct ScanInfo
    {
      Int charge;
      double rt, mz;
    };

    /// Mapping: spectrum index -> fragment spectrum details
    typedef std::map<Size, struct ScanInfo> ScanInfoMap;

    /// Description of how to extract information from PSM IDs
    struct PSMInfoExtractor
    {
      boost::regex re; // regular expression for meta data extraction
      bool count_from_zero; // start counting scans at 0 or 1?
    };

    /// List of data extractors to try by default
    std::vector<struct PSMInfoExtractor> extractors_;

    /**
       @brief Extract information from a Percolator PSM ID

       @return Returns whether extraction was successful
    */
    bool getPSMInfo_(const String& PSM_ID,
                     const std::vector<struct PSMInfoExtractor>& extractors,
                     Int& scan_number, Int& charge, double& rt, double& mz);

    /// Converts the peptide string to an 'AASequence' instance
    void getPeptideSequence_(String peptide, AASequence& seq) const;

    /// Extracts information from the raw data
    void preprocessExperiment_(const MSExperiment<>& experiment,
                               ScanInfoMap& scan_map);
  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_PERCOLATOROUTFILE_H
