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

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <vector>
#include <map>
#include <set>
#include <algorithm>

namespace OpenMS
{
  /**

  @brief Representation of the Experimental Design in OpenMS. Instances can be loaded via
         the ExperimentalDesignFile class.

  The experimental design in OpenMS consists of two sections:
    1. FileSection captures the mapping of quantitative values to files and, optionally, to samples.
    2. SampleSection captures the experimental factors and conditions associated with a sample.

  The FileSection is mandatory while the SampleSection is optional and only required for downstream analysis.

  Details on the FileSection:

  To unambiguously map a single quantitative value we need to define:
    a. the label (e.g., label = 2 in the case of a heavy peptide in a light/heavy experiment)
    b. which spectra file did produce the result (e.g., path = "/data/SILAC_file.mzML")  
    c. which fraction the file corresponds to (e.g., fraction = 1)
    d. a fraction group identifier that groups fractions (e.g., fraction_group = 1)
       Note: in the case of label-free data, the fraction group identifier has 
             the same cardinality as the sample identifier.
   (e.) optionally, the sample that has been measured (e.g., sample = 1)

   This information is defined in the FileSection and enables fraction aware data processing.

  Details on the SampleSection:
  
  To map a sample to conditions / factors we need to define:
    a. the sample ( e.g., sample = 1)
    b. multiple columns containing conditions / factors

  @ingroup Metadata

  **/

  class OPENMS_DLLAPI ExperimentalDesign
  {

  public:
    /// MSFileSectionEntry links single quant. values back the MS file
    /// It supports:
    ///  - multiplexed/labeled data via specification of the quantified label
    ///  - multiple fractions via specification of the:
    ///    - fraction index (e.g., 1..10 if ten fractions were measured)
    ///    - fraction_group to trace which fractions belong together
    class OPENMS_DLLAPI MSFileSectionEntry
    {
    public:
      MSFileSectionEntry() = default;
      unsigned fraction_group = 1; ///< fraction group id
      unsigned fraction = 1; ///< fraction 1..m, mandatory, 1 if not set
      std::string path = "UNKNOWN_FILE"; ///< file name, mandatory
      unsigned label = 1;  ///< the label (e.g.,: 1 for label-free, 1..8 for TMT8plex)
      unsigned sample = 1;  ///< allows grouping by sample
    };

    class OPENMS_DLLAPI SampleSection
    {
    public:

      SampleSection() = default;

      SampleSection(
        const std::vector< std::vector < String > >& content,
        const std::map< unsigned, Size >& sample_to_rowindex,
        const std::map< String, Size >& columnname_to_columnindex
      );

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
    };

    using MSFileSection = std::vector<MSFileSectionEntry>;

    // Experimental Design c'tors
    ExperimentalDesign() = default;

    ExperimentalDesign(const MSFileSection& msfile_section, const SampleSection& sample_section);

    const MSFileSection& getMSFileSection() const;

    void setMSFileSection(const MSFileSection& msfile_section);

    // Returns the Sample Section of the experimental design file
    const ExperimentalDesign::SampleSection& getSampleSection() const;

    void setSampleSection(const SampleSection& sample_section);

    /// return fraction index to file paths (ordered by fraction_group)
    std::map<unsigned int, std::vector<String> > getFractionToMSFilesMapping() const;

   /*
    *   The (Path, Label) tuples in the experimental design have to be unique, so we can map them
    *   uniquely to the sample number, fraction number, and fraction_group number
    */
    /// return <file_path, label> to sample mapping
    std::map< std::pair< String, unsigned >, unsigned> getPathLabelToSampleMapping(bool) const;

    /// return <file_path, label> to fraction mapping
    std::map< std::pair< String, unsigned >, unsigned> getPathLabelToFractionMapping(bool) const;

    /// return <file_path, label> to fraction_group mapping
    std::map< std::pair< String, unsigned >, unsigned> getPathLabelToFractionGroupMapping(bool) const;

    // @return the number of samples measured (= highest sample index)
    unsigned getNumberOfSamples() const;

    // @return the number of fractions (= highest fraction index)
    unsigned getNumberOfFractions() const;

    // @return the number of labels per file
    unsigned getNumberOfLabels() const;

    // @return the number of MS files (= fractions * fraction groups)
    unsigned getNumberOfMSFiles() const;

    // @return the number of fraction_groups
    // Allows to group fraction ids and source files
    unsigned getNumberOfFractionGroups() const;

    // @return sample index (depends on fraction_group and label)
    unsigned getSample(unsigned fraction_group, unsigned label = 1);

    /// @return whether we have a fractionated design 
    // This is the case if we have at least one fraction group with >= 2 fractions
    bool isFractionated() const;

    /// @returns whether all fraction groups have the same number of fractions
    bool sameNrOfMSFilesPerFraction() const;

    /// Extract experimental design from consensus map
    static ExperimentalDesign fromConsensusMap(const ConsensusMap& c);

    /// Extract experimental design from feature map
    static ExperimentalDesign fromFeatureMap(const FeatureMap& f);

    /// Extract experimental design from identifications
    static ExperimentalDesign fromIdentifications(const std::vector<ProteinIdentification> & proteins);

    private:
    // MS filename column, optionally trims to basename
    std::vector< String > getFileNames_(bool basename) const;

    // returns label column
    std::vector<unsigned> getLabels_() const;

    // returns fraction column
    std::vector<unsigned> getFractions_() const;

    /// Generic Mapper (Path, Label) -> f(row)
    std::map< std::pair< String, unsigned >, unsigned> pathLabelMapper_(
        bool,
        unsigned (*f)(const ExperimentalDesign::MSFileSectionEntry&)) const;

    // sort to obtain the default order
    void sort_();

    template<typename T>
    static void errorIfAlreadyExists(std::set<T> &container, T &item, const String &message);

    // basic consistency checks
    void isValid_();

    MSFileSection msfile_section_;
    SampleSection sample_section_;
  };
}

