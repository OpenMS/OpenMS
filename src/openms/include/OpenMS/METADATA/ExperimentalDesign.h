// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:	Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <vector>
#include <map>
#include <set>

namespace OpenMS
{
  class ConsensusMap;
  class FeatureMap;

  /**

  @brief Representation of an experimental design in OpenMS. Instances can be loaded with
         the ExperimentalDesignFile class.

    Experimental designs can be provided in two formats: the one-table format and the two-table format.

    The one-table format is simpler but slightly more redundant.

    The one-table format consists of mandatory (file columns) and optional sample metadata (sample columns).

    The mandatory file columns are Fraction_Group, Fraction, Spectra_Filepath and Label.
    These columns capture the mapping of quantitative values to files for label-free and multiplexed experiments
    and enables fraction-aware data processing.

    - Fraction_Group: a numeric identifier that indicates which fractions are grouped together. Please do NOT
        reuse the same identifiers across samples! Assign identifiers continuously.

    - Fraction: a numeric identifier that indicates which fraction was measured in this file. 
              In the case of unfractionated data, the fraction identifier is 1 for all samples.
              Make sure the same identifiers are used across different Fraction_Groups, as this
              determines which fractions correspond to each other.

    - Label: a numeric identifier for the label. 1 for label-free, 1 and 2 for SILAC light/heavy, e.g., 1-10 for TMT10Plex

    - Spectra_Filepath: a filename or path as string representation (e.g., SILAC_file.mzML)

    For processing with MSstats, the optional sample columns are typically MSstats_Condition and MSstats_BioReplicate with an additional MSstats_Mixture
    column in the case of TMT labeling.
    They capture the experimental factors and conditions associated with a sample.

    - MSstats_Condition: a string that indicates the condition (e.g., control or 1000 mMol). Will be forwarded to MSstats and
                       can then be used to specify test contrasts.

    - MSstats_BioReplicate: a string identifier to indicate biological replication of a sample.
                        Entries with the same Sample/Condition/BioReplicate but different Filepath (and therefore FractionGroup number)
                        will be treated as technical replicates.
                          
    - MSstats_Mixture: (for TMT labeling only): a string identifier to indicate the mixture of samples labeled with different TMT reagents, which can be analyzed in
                                             a single mass spectrometry experiment. E.g., same samples labeled with different TMT reagents have a different mixture identifier. 
                                             Technical replicates need to have the same mixture identifier.

    For details on the MSstats columns please refer to the MSstats manual for details
    (https://www.bioconductor.org/packages/release/bioc/vignettes/MSstats/inst/doc/MSstats.html).

    <table>
    <tr>
        <th>Fraction_Group</th>
        <th>Fraction</th>
        <th>Spectra_Filepath</th>
        <th>Label</th>
        <th>MSstats_Condition</th>
        <th>MSstats_BioReplicate</th>
    </tr>
    <tr>
        <td>1</td>
        <td>1</td>
        <td>UPS1_12500amol_R1.mzML</td>
        <td>1</td>
        <td>12500 amol</td>
        <td>1</td>
    </tr>
    <tr>
        <td>2</td>
        <td>1</td>
        <td>UPS1_12500amol_R2.mzML</td>
        <td>1</td>
        <td>12500 amol</td>
        <td>2</td>
    </tr>
    <tr>
        <td>3</td>
        <td>1</td>
        <td>UPS1_12500amol_R3.mzML</td>
        <td>1</td>
        <td>12500 amol</td>
        <td>3</td>
    </tr>
    <tr>
        <td>...</td>
        <td>...<br></td>
        <td>...</td>
        <td>...<br></td>
        <td>...<br></td>
        <td>...<br></td>
    </tr>
    <tr>
        <td>22</td>
        <td>1</td>
        <td>UPS1_500amol_R1.mzML</td>
        <td>1</td>
        <td>500 amol</td>
        <td>1</td>
    </tr>
    <tr>
        <td>23</td>
        <td>1</td>
        <td>UPS1_500amol_R2.mzML</td>
        <td>1</td>
        <td>500 amol</td>
        <td>2</td>
    </tr>
    <tr>
        <td>24</td>
        <td>1</td>
        <td>UPS1_500amol_R3.mzML</td>
        <td>1</td>
        <td>500 amol</td>
        <td>3</td>
    </tr>
    </table>

    Alternatively, the experimental design can be specified with a file consisting of two tables whose headers are separated
    by a blank line. The two tables are:

    - The file section table and the sample section table.
    - The file section consists of columns Fraction_Group, Fraction, Spectra_Filepath, Label and Sample

    The sample section consists of columns Sample, MSstats_Condition and MSstats_BioReplicate.

    The content is the same as described for the one table format, except that the additional numeric sample column
    allows referencing between file and sample section.

    <table>
    <tr>
        <th>Fraction_Group</th>
        <th>Fraction</th>
        <th>Spectra_Filepath</th>
        <th>Label</th>
        <th>Sample</th>
    </tr>
    <tr>
        <td>1</td>
        <td>1</td>
        <td>UPS1_12500amol_R1.mzML</td>
        <td>1</td>
        <td>1</td>
    </tr>
    <tr>
        <td>2</td>
        <td>1</td>
        <td>UPS1_12500amol_R2.mzML</td>
        <td>1</td>
        <td>2</td>
    </tr>
    <tr>
        <td>...</td>
        <td>...<br></td>
        <td>...</td>
        <td>...<br></td>
        <td>...<br></td>
    </tr>
    <tr>
        <td>22</td>
        <td>1</td>
        <td>UPS1_500amol_R1.mzML</td>
        <td>1</td>
        <td>22</td>
    </tr>
    </table>

    <table>
    <tr>
        <th>Sample</th>
        <th>MSstats_Condition</th>
        <th>MSstats_BioReplicate</th>        
    </tr>
    <tr>
        <td>1</td>
        <td>12500 amol</td>
        <td>1</td>
    </tr>
    <tr>
        <td>2</td>
        <td>12500 amol</td>
        <td>2</td>
    </tr>
    <tr>
        <td>...</td>
        <td>...<br></td>
        <td>...<br></td>
    </tr>
    <tr>
        <td>22</td>
        <td>500 amol</td>
        <td>3</td>
    </tr>
    </table>


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
      unsigned sample = 0;  ///< zero-based index for the sample section contents; allows grouping by sample
      String sample_name = "0"; ///< sample name for access of the sample row by string, not index
    };

    class OPENMS_DLLAPI SampleSection
    {
    public:

      SampleSection() = default;

      SampleSection(
        const std::vector< std::vector < String > >& content,
        const std::map< String, Size >& sample_to_rowindex,
        const std::map< String, Size >& columnname_to_columnindex
      );

      // Get set of all samples that are present in the sample section
      std::set< String > getSamples() const;

      // Add a sample as the last row
      void addSample(const String& sample, const std::vector<String>& content = {});

      // TODO should it include the Sample ID column or not??
      // Get set of all factors (column names) that were defined for the sample section
      std::set< String > getFactors() const;

      // Checks whether sample section has row for a sample number
      bool hasSample(const String& sample) const;

      // Checks whether Sample Section has a specific factor (i.e. column name)
      bool hasFactor(const String &factor) const;

      // Returns value of factor for given sample and factor name
      String getFactorValue(const String& sample_name, const String &factor) const;

      // Returns value of factor for given sample index and factor name
      String getFactorValue(unsigned sample_idx, const String &factor) const;

      // Returns column index of factor
      Size getFactorColIdx(const String &factor) const;

      // Returns the name/ID of the sample. Not necessarily the row index
      String getSampleName(unsigned sample_row) const;

      // Returns the row index in the sample section for a sample name/ID
      unsigned getSampleRow(const String& sample) const;

      /// returns the number of entries in content_ member
      Size getContentSize() const;

    private:

      // The entries of the Sample Section, filled while parsing
      // the Experimental Design File
      std::vector< std::vector < String > > content_;

      // Maps the Sample Entry name to the row where the sample
      // appears in the Sample section, its sample index
      std::map< String, Size > sample_to_rowindex_;

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

    /// returns a map from a sample section row to sample id for clustering
    /// duplicate sample rows (e.g. to find all fractions of the same "sample")
    std::map<std::vector<String>, std::set<String>> getUniqueSampleRowToSampleMapping() const;

    /// uses getUniqueSampleRowToSampleMapping to get the reversed map
    /// mapping sample ID to a real unique sample
    std::map<String, unsigned> getSampleToPrefractionationMapping() const;

    /// return fraction index to file paths (ordered by fraction_group)
    //TODO this probably needs a basename parameter to be fully compatible with the other mappings!! Implicit full path.
    std::map<unsigned int, std::vector<String> > getFractionToMSFilesMapping() const;

    /// return vector of filepath/label combinations that share the same conditions after removing
    /// replicate columns in the sample section (e.g. for merging across replicates)
    //TODO this probably needs a basename parameter to be fully compatible with the other mappings!! Implicit full path.
    std::vector<std::vector<std::pair<String, unsigned>>> getConditionToPathLabelVector() const;

    /// return a condition (unique combination of sample section values except replicate) to Sample index mapping
    std::map<std::vector<String>, std::set<unsigned>> getConditionToSampleMapping() const;

   /*
    *   The (Path, Label) tuples in the experimental design have to be unique, so we can map them
    *   uniquely to the sample number, fraction number, and fraction_group number
    */

    /// return <file_path, label> to prefractionation mapping (a prefractionation group is a unique combination of
    /// all columns in the sample section, except for replicates.
    std::map< std::pair< String, unsigned >, unsigned> getPathLabelToPrefractionationMapping(bool use_basename_only) const;

    /// return <file_path, label> to condition mapping (a condition is a unique combination of all columns in the
    /// sample section, except for replicates.
    std::map< std::pair< String, unsigned >, unsigned> getPathLabelToConditionMapping(bool use_basename_only) const;

    /// return Sample name to condition mapping (a condition is a unique combination of all columns in the
    /// sample section, except for replicates. Numbering of conditions is alphabetical due to map.
    std::map<String, unsigned> getSampleToConditionMapping() const;

    /// return <file_path, label> to sample index mapping
    std::map< std::pair< String, unsigned >, unsigned> getPathLabelToSampleMapping(bool use_basename_only) const;

    /// return <file_path, label> to fraction mapping
    std::map< std::pair< String, unsigned >, unsigned> getPathLabelToFractionMapping(bool use_basename_only) const;

    /// return <file_path, label> to fraction_group mapping
    std::map< std::pair< String, unsigned >, unsigned> getPathLabelToFractionGroupMapping(bool use_basename_only) const;

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

    /// filters the MSFileSection to only include a given subset of files whose basenames
    /// are given with @p bns
    /// @return number of files that have been filtered
    Size filterByBasenames(const std::set<String>& bns);

    /// @returns whether all fraction groups have the same number of fractions
    bool sameNrOfMSFilesPerFraction() const;

    /// Extract experimental design from consensus map
    static ExperimentalDesign fromConsensusMap(const ConsensusMap& c);

    /// Extract experimental design from feature map
    static ExperimentalDesign fromFeatureMap(const FeatureMap& f);

    /// Extract experimental design from identifications
    static ExperimentalDesign fromIdentifications(const std::vector<ProteinIdentification>& proteins);
    //TODO create another overload here, that takes two enums outerVec and innerVec with entries Replicate, Fraction, Sample

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

