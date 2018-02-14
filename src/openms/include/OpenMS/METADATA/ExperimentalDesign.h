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
// $Maintainer:	Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_EXPERIMENTALDESIGN_H
#define OPENMS_KERNEL_EXPERIMENTALDESIGN_H

#include <OpenMS/DATASTRUCTURES/String.h>

#include <vector>
#include <string>
#include <map>
#include <set>
#include <algorithm>

namespace OpenMS
{
  /**
  @brief A TSV and user friendly representation of the experimental design.
  Used for loading and storing the experimental design in OpenMS.

  @ingroup Metadata
  */
  class OPENMS_DLLAPI ExperimentalDesign
  {
  public:
    /// 1) Mandatory section with run-level information of the experimental design.
    ///    Required to process fractionated data and technical replicates.
/*
Run	Fraction	Path (Spectra File)	Assay (Quantification Channel)	Sample (Condition)
1	1	SPECTRAFILE_F1_TR1.mzML	iTRAQ reagent 114	1
1	2	SPECTRAFILE_F2_TR1.mzML	iTRAQ reagent 114	1
1	3	SPECTRAFILE_F3_TR1.mzML	iTRAQ reagent 114	1
1	1	SPECTRAFILE_F1_TR1.mzML	iTRAQ reagent 115	2
1	2	SPECTRAFILE_F2_TR1.mzML	iTRAQ reagent 115	2
1	3	SPECTRAFILE_F3_TR1.mzML	iTRAQ reagent 115	2
1	1	SPECTRAFILE_F1_TR1.mzML	iTRAQ reagent 116	3
1	2	SPECTRAFILE_F2_TR1.mzML	iTRAQ reagent 116	3
1	3	SPECTRAFILE_F3_TR1.mzML	iTRAQ reagent 116	3
1	1	SPECTRAFILE_F1_TR1.mzML	iTRAQ reagent 117	4
1	2	SPECTRAFILE_F2_TR1.mzML	iTRAQ reagent 117	4
1	3	SPECTRAFILE_F3_TR1.mzML	iTRAQ reagent 117	4
2	1	SPECTRAFILE_F1_TR2.mzML	iTRAQ reagent 114	5
2	2	SPECTRAFILE_F2_TR2.mzML	iTRAQ reagent 114	5
2	3	SPECTRAFILE_F3_TR2.mzML	iTRAQ reagent 114	5
2	1	SPECTRAFILE_F1_TR2.mzML	iTRAQ reagent 115	6
2	2	SPECTRAFILE_F2_TR2.mzML	iTRAQ reagent 115	6
2	3	SPECTRAFILE_F3_TR2.mzML	iTRAQ reagent 115	6
2	1	SPECTRAFILE_F1_TR2.mzML	iTRAQ reagent 116	7
2	2	SPECTRAFILE_F2_TR2.mzML	iTRAQ reagent 116	7
2	3	SPECTRAFILE_F3_TR2.mzML	iTRAQ reagent 116	7
2	1	SPECTRAFILE_F1_TR2.mzML	iTRAQ reagent 117	8
2	2	SPECTRAFILE_F2_TR2.mzML	iTRAQ reagent 117	8
2	3	SPECTRAFILE_F3_TR2.mzML	iTRAQ reagent 117	8
*/
    class OPENMS_DLLAPI Row
    {
    public:
      Row() = default;
      unsigned run = 1; ///< run index (before prefractionation)
      unsigned fraction = 1; ///< fraction 1..m, mandatory, 1 if not set
      std::string path = "NA"; ///< file name, mandatory
      std::string assay = "label-free";  ///< needed to determine if and how many multiplexed channels are in a file
      unsigned sample = 1;  ///< allows grouping by sample
    };

    // the rows of the experimental design
    std::vector<Row> rows;

    /// return fraction index to file paths
    std::map<unsigned int, std::set<String> > getFractionToMSFilesMapping() const;

    // @return the number of samples measured (= highest sample index)
    unsigned getNumberOfSamples() const 
    {
      return std::max_element(rows.begin(), rows.end(), 
        [](const Row& f1, const Row& f2) 
        {
          return f1.sample < f2.sample;
        })->sample;
    }

    // @return the number of fractions (= highest fraction index)
    unsigned getNumberOfFractions() const 
    {
      return std::max_element(rows.begin(), rows.end(), 
        [](const Row& f1, const Row& f2) 
        {
          return f1.fraction < f2.fraction;
        })->fraction;
    }

    // @return the number of MS files (= fractions * runs)
    unsigned getNumberOfMSFiles() const
    {
      std::set<std::string> unique_paths;
      for (auto const & r : rows) { unique_paths.insert(r.path); }
      return unique_paths.size();
    }

    // @return the number of runs (before fractionation)
    unsigned getNumberOfRuns() const
    {
      // TODO: assert fractions * runs == max_element ...
      return std::max_element(rows.begin(), rows.end(), 
        [](const Row& f1, const Row& f2) 
        {
          return f1.fraction < f2.fraction;
        })->run;
    }

    /// return if each fraction number is associated with the same number of runs 
    bool sameNrOfMSFilesPerFraction() const;

    /// Loads an experimental design from a tabular separated file
    void load(const String & tsv_file, ExperimentalDesign & design) const;
  };
}


#endif // header guard
