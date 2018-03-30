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
// $Maintainer:	Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_EXPERIMENTALDESIGN_H
#define OPENMS_KERNEL_EXPERIMENTALDESIGN_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/SYSTEM/File.h>

#include <vector>
#include <map>
#include <set>

namespace OpenMS
{
  /**
  @brief Representation of the Experimental Design in OpenMS. Instances are loaded via
   the ExperimentalDesignIO class.

  @ingroup Format
  */
  class OPENMS_DLLAPI ExperimentalDesign
  {
  public:
  ExperimentalDesign() = default;

    /// 1) Mandatory section with run-level information of the experimental design.
    ///    Required to process fractionated data.
/*
 * Run Section Format:
   Format: Single header line
         Run:                         Run index (prior fractionation) used to group fractions and source files.
                                      Note: For label-free this has same cardinality as sample.
                                      For multiplexed experiments, these might differ as multiple samples can be measured in single files
         Fraction:                    1st, 2nd, .., fraction. Note: All runs must have the same number of fractions.
         Path(Spectra File):          Path to mzML files
         Channel:                     Channel in MS file:
                                      label-free: always 1
                                      TMT6Plex: 1..6
                                      SILAC with light and heavy: 1..2
         Sample:                      Index of sample measured in the specified channel X, in fraction Y of run Z

	Run	Fraction	Path(Spectra File)	Channel		Sample
	1	1		SPECTRAFILE_F1_TR1.mzML	1		1
	1	2		SPECTRAFILE_F2_TR1.mzML	1		1
	1	3		SPECTRAFILE_F3_TR1.mzML	1		1
	1	1		SPECTRAFILE_F1_TR1.mzML	2		2
	1	2		SPECTRAFILE_F2_TR1.mzML	2		2
	1	3		SPECTRAFILE_F3_TR1.mzML	2		2
	1	1		SPECTRAFILE_F1_TR1.mzML	3		3
	1	2		SPECTRAFILE_F2_TR1.mzML	3		3
	1	3		SPECTRAFILE_F3_TR1.mzML	3		3
	1	1		SPECTRAFILE_F1_TR1.mzML	4		4
	1	2		SPECTRAFILE_F2_TR1.mzML	4		4
	1	3		SPECTRAFILE_F3_TR1.mzML	4		4
	2	1		SPECTRAFILE_F1_TR2.mzML	1		5
	2	2		SPECTRAFILE_F2_TR2.mzML	1		5
	2	3		SPECTRAFILE_F3_TR2.mzML	1		5
	2	1		SPECTRAFILE_F1_TR2.mzML	2		6
	2	2		SPECTRAFILE_F2_TR2.mzML	2		6
	2	3		SPECTRAFILE_F3_TR2.mzML	2		6
	2	1		SPECTRAFILE_F1_TR2.mzML	3		7
	2	2		SPECTRAFILE_F2_TR2.mzML	3		7
	2	3		SPECTRAFILE_F3_TR2.mzML	3		7
	2	1		SPECTRAFILE_F1_TR2.mzML	4		8
	2	2		SPECTRAFILE_F2_TR2.mzML	4		8
	2	3		SPECTRAFILE_F3_TR2.mzML	4		8

  /// 2) Mandatory section with sample information of the experimental design.
  ///    Required to process fractionated data. One Column must be 'Sample', other columns
  ///    are unspecified and can contain arbitrary factors

 Sample	Some_Condition	Technical_Replicate
  1     1               1
  2	    2	              1
  3	    3	              1
  4	    4	              1
  5	    1	              2
  6	    2	              2
  7	    3	              2
  8	    4	              2

*/
    class OPENMS_DLLAPI RunRow
    {
    public:
      RunRow() = default;
      unsigned run = 1; ///< run index (before prefractionation)
      unsigned fraction = 1; ///< fraction 1..m, mandatory, 1 if not set
      std::string path = "UNKNOWN_FILE"; ///< file name, mandatory
      unsigned channel = 1;  ///< if and how many multiplexed channels are in a file
      unsigned sample = 1;  ///< allows grouping by sample
    };

    class OPENMS_DLLAPI SampleSection
    {
    public:
      SampleSection() = default;

      // Get set of all samples that are present in the sample section
      std::set< unsigned > getSamples() const;

      // Get set of all factors (column names) that were defined for the sample section
      std::set< String > getFactors() const;

      // Checks whether sample section has row for a sample number
      bool hasSample(unsigned sample) const;

      // Checks whether Sample Section has a specific factor (i.e. column name)
      bool hasFactor(const String &factor) const;

      // Returns value of factor for given sample and factor name
      String getFactorValue(unsigned sample, const String &factor);

    private:

      // The entries of the Sample Section, filled while parsing
      // the Experimental Design File
      std::vector< std::vector < String > > content_;

      // Maps the Sample Entry to the row where the sample
      // appears in the Sample section
      std::map< unsigned, Size > sample_to_rowindex_;

      // Maps the column name of the SampleSection to the
      // Index of the column
      std::map< String, Size > columnname_to_columnindex_;

      friend class ExperimentalDesignIO;
    };

    using RunRows = std::vector<RunRow>;

    const RunRows& getRunSection() const;

    void setRunSection(const RunRows& run_section);

    // Returns the Sample Section of the experimental design file
    const ExperimentalDesign::SampleSection& getSampleSection() const;

    // Gets vector of Filenames that appears in the run section, optionally trims to basename
    std::vector< String > getFileNames(bool basename) const;

    // Returns vector of channels of the run section
    std::vector<unsigned> getChannels() const;

    std::vector<unsigned> getFractions() const;

    /// return fraction index to file paths (ordered by run id)
    std::map<unsigned int, std::vector<String> > getFractionToMSFilesMapping() const;

   /*
    *   The (Path, Channel) tuples in the experimental design have to be unique, so we can map them
    *   uniquely to the sample number, fraction number, and run number
    */
    /// return <file_path, channel> to sample mapping
    std::map< std::pair< String, unsigned >, unsigned> getPathChannelToSampleMapping(bool) const;

    /// return <file_path, channel> to fraction mapping
    std::map< std::pair< String, unsigned >, unsigned> getPathChannelToFractionMapping(bool) const;

    /// return <file_path, channel> to run mapping
    std::map< std::pair< String, unsigned >, unsigned> getPathChannelToRunMapping(bool) const;

    // @return the number of samples measured (= highest sample index)
    unsigned getNumberOfSamples() const;

    // @return the number of fractions (= highest fraction index)
    unsigned getNumberOfFractions() const;

    // @return the number of channels per file
    unsigned getNumberOfChannels() const;

    // @return the number of MS files (= fractions * runs)
    unsigned getNumberOfMSFiles() const;

    // @return the number of runs (before fractionation)
    // Allows to group fraction ids and source files
    unsigned getNumberOfPrefractionationRuns() const;

    // @return sample index (depends on run and channel)
    unsigned getSample(unsigned run, unsigned channel = 1);

    /// @return whether at least one run in this experimental design is fractionated
    bool isFractionated() const;

    /// return if each fraction number is associated with the same number of runs
    bool sameNrOfMSFilesPerFraction() const;

    friend class ExperimentalDesignIO;

    private:

      /// Generic Mapper (Path, Channel) -> f(row)
      std::map< std::pair< String, unsigned >, unsigned> pathChannelMapper(
          bool,
          unsigned (*f)(const ExperimentalDesign::RunRow&)) const;

      // sort to obtain the default order
      void sort_();

      template<typename T>
      static void errorIfAlreadyExists(std::set<T> &container, T &item, const String &message);

      void checkValidRunSection_();

      RunRows run_section_;
      SampleSection sample_section_;
  };
}
#endif // header guard
