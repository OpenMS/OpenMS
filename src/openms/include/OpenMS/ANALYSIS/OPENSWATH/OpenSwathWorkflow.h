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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_OPENSWATH_OPENSWATHWORKFLOW_H
#define OPENMS_ANALYSIS_OPENSWATH_OPENSWATHWORKFLOW_H

// Interfaces
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/DataStructures.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>

#include <OpenMS/FORMAT/MzMLFile.h> // debug file store only

// Kernel and implementations
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessTransforming.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMSInMemory.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/SwathMap.h>

// Helpers
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>
// #include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathTSVWriter.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathOSWWriter.h>

// Algorithms
#include <OpenMS/ANALYSIS/OPENSWATH/MRMRTNormalizer.h>
#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SwathMapMassCorrection.h>
#include <OpenMS/FILTERING/TRANSFORMERS/LinearResamplerAlign.h>

#include <assert.h>
#include <limits>

// #define OPENSWATH_WORKFLOW_DEBUG

// The workflow class
namespace OpenMS
{

  /** @brief ChromatogramExtractor parameters
   *
   * A small helper struct to pass the parameters for the chromatogram
   * extraction through to the actual algorithm.
   *
  */
  struct ChromExtractParams
  {
    /// Whether to not extract anything closer than this (in Da) from the upper edge
    double min_upper_edge_dist;
    /// Extraction window in Da or ppm (e.g. 50ppm means extraction +/- 25ppm)
    double mz_extraction_window;
    /// Whether the extraction window is given in ppm or Da
    bool ppm;
    /// The extraction function in mass space
    String extraction_function;
    /// The retention time extraction window
    double rt_extraction_window;
    /// Whether to extract some extra in the retention time (can be useful if one wants to look at the chromatogram outside the window)
    double extra_rt_extract;
  };

  /** @brief Simple OpenSwathWorkflow to perform RT and m/z correction based on a set of known peptides
   *
  */
  class OPENMS_DLLAPI OpenSwathRetentionTimeNormalization :
    public ProgressLogger
  {
  public:

    /** @brief Perform RT and m/z correction of the input data using RT-normalization peptides.
     *
     * This function extracts the RT normalization chromatograms and then uses
     * the chromatograms to find features (in RTNormalization).  If desired,
     * also m/z correction is performed using the lock masses of the given
     * peptides. Therefore, swath_maps may be changed in this function.
     *
     * @param irt_transitions A set of transitions used for the RT normalization peptides
     * @param swath_maps The raw data (swath maps)
     * @param min_rsq Minimal R^2 value that is expected for the RT regression
     * @param min_coverage Minimal coverage of the chromatographic space that needs to be achieved
     * @param feature_finder_param Parameter set for the feature finding in chromatographic dimension
     * @param cp_irt Parameter set for the chromatogram extraction
     * @param irt_detection_param Parameter set for the detection of the iRTs (outlier detection, peptides per bin etc)
     * @param mz_correction_function If correction in m/z is desired, which function should be used
     * @param debug_level Debug level (writes out the RT normalization chromatograms if larger than 1)
     *
    */
    TransformationDescription performRTNormalization(const OpenMS::TargetedExperiment & irt_transitions,
      std::vector< OpenSwath::SwathMap > & swath_maps,
      double min_rsq,
      double min_coverage,
      const Param & feature_finder_param,
      const ChromExtractParams & cp_irt,
      const Param & irt_detection_param,
      const String & mz_correction_function,
      Size debug_level,
      bool sonar = false,
      bool load_into_memory = false);

  private:

    /** @brief Perform RT and m/z correction using the MRMFeatureFinderScoring
     *
     * @param transition_exp_ The transitions for the normalization peptides
     * @param chromatograms The extracted chromatograms
     * @param min_rsq Minimal R^2 value that is expected for the RT regression
     * @param min_coverage Minimal coverage of the chromatographic space that needs to be achieved
     * @param feature_finder_param Parameter set for the feature finding in chromatographic dimension
     * @param irt_detection_param Parameter set for the detection of the iRTs (outlier detection, peptides per bin etc)
     * @param swath_maps The raw data for the m/z correction
     * @param mz_correction_function If correction in m/z is desired, which function should be used
     *
     * @note: feature_finder_param are copied because they are changed here.
     * @note: This function is based on the algorithm inside the OpenSwathRTNormalizer tool
     *
    */
    TransformationDescription RTNormalization(const TargetedExperiment& transition_exp_,
      const std::vector< OpenMS::MSChromatogram >& chromatograms,
      double min_rsq,
      double min_coverage,
      const Param& default_ffparam,
      const Param& irt_detection_param,
      std::vector< OpenSwath::SwathMap > & swath_maps,
      const String & mz_correction_function,
      double mz_extraction_window,
      bool ppm);

    /// Simple method to extract chromatograms (for the RT-normalization peptides)
    void simpleExtractChromatograms(const std::vector< OpenSwath::SwathMap > & swath_maps,
                                    const OpenMS::TargetedExperiment & irt_transitions,
                                    std::vector< OpenMS::MSChromatogram > & chromatograms,
                                    const ChromExtractParams & cp, bool sonar, bool load_into_memory);

    static void addChromatograms(MSChromatogram& base_chrom, const MSChromatogram& newchrom);
  };

  /**
   * @brief Class to execute an OpenSwath Workflow
   *
   * The workflow will perform a complete OpenSWATH analysis. Optionally, an RT
   * transformation (mapping peptides to normalized space) can be obtained
   * beforehand using the OpenSwathRetentionTimeNormalization class.
   *
   */
  class OPENMS_DLLAPI OpenSwathWorkflow :
    public ProgressLogger
  {
    typedef OpenSwath::LightTransition TransitionType;
    typedef MRMTransitionGroup< MSChromatogram, TransitionType> MRMTransitionGroupType;

  public:

    explicit OpenSwathWorkflow(bool use_ms1_traces) :
      use_ms1_traces_(use_ms1_traces)
    {
    }

    /** @brief Execute the OpenSWATH workflow on a set of SwathMaps and transitions.
     *
     * Executes the following operations on the given input:
     *
     * 1. Selecting the appropriate transitions for each SWATH window (using OpenSwathHelper::selectSwathTransitions)
     * 2. Extract the chromatograms from the SWATH maps (MS1 and MS2) using (ChromatogramExtractor)
     * 3. Pick peaks in the chromatograms and perform peak scoring (inside scoreAllChromatograms function)
     * 4. Write out chromatograms and found features
     *
     * @param swath_maps The raw data (swath maps)
     * @param trafo Transformation description (translating this runs' RT to normalized RT space)
     * @param cp Parameter set for the chromatogram extraction
     * @param feature_finder_param Parameter set for the feature finding in chromatographic dimension
     * @param transition_exp The set of assays to be extracted and scored
     * @param out_featureFile Output feature map to store identified features
     * @param store_features Whether features should be appended to the output feature map
     * @param tsv_writer TSV Writer object to store identified features in csv format
     * @param osw_writer OSW Writer object to store identified features in SQLite format
     * @param chromConsumer Chromatogram consumer object to store the extracted chromatograms
     * @param batchSize Size of the batches which should be extracted and scored
     * @param load_into_memory Whether to cache the current SWATH map in memory
     *
    */
    void performExtraction(const std::vector< OpenSwath::SwathMap > & swath_maps,
                           const TransformationDescription trafo,
                           const ChromExtractParams & cp,
                           const Param & feature_finder_param,
                           const OpenSwath::LightTargetedExperiment& transition_exp,
                           FeatureMap& out_featureFile,
                           bool store_features,
                           OpenSwathTSVWriter & tsv_writer,
                           OpenSwathOSWWriter & osw_writer,
                           Interfaces::IMSDataConsumer * chromConsumer,
                           int batchSize,
                           bool load_into_memory);

  protected:


    /** @brief Write output features and chromatograms to disk
     *
    */
    void writeOutFeaturesAndChroms_(std::vector< OpenMS::MSChromatogram > & chromatograms,
                                    const FeatureMap & featureFile,
                                    FeatureMap& out_featureFile,
                                    bool store_features,
                                    Interfaces::IMSDataConsumer * chromConsumer);

    /** @brief Perform MS1 extraction and store result in ms1_chromatograms
     *
    */
    void MS1Extraction_(const std::vector< OpenSwath::SwathMap > & swath_maps,
                        std::map< std::string, OpenSwath::ChromatogramPtr >& ms1_chromatograms,
                        Interfaces::IMSDataConsumer * chromConsumer,
                        const ChromExtractParams & cp,
                        const OpenSwath::LightTargetedExperiment& transition_exp,
                        const TransformationDescription& trafo_inverse,
                        bool load_into_memory,
                        bool ms1only = false);

    /** @brief Perform scoring on a set of chromatograms
     *
     *  Will iterate over all assays contained in transition_exp and for each
     *  assay fetch the corresponding chromatograms and find peak groups.
     *
     * @param input Input chromatograms (MS2 level)
     * @param ms1_chromatograms Input chromatograms for MS1-level
     * @param swath_maps Set of swath map(s) for the current swath window (for SONAR multiple maps are provided)
     * @param transition_exp The transition experiment (assay library)
     * @param feature_finder_param Parameters for the MRMFeatureFinderScoring
     * @param trafo RT Transformation function
     * @param rt_extraction_window RT extraction window
     * @param output Output map
     * @param tsv_writer TSV writer for storing output (on the fly)
     * @param osw_writer OSW Writer object to store identified features in SQLite format
     *
    */
    void scoreAllChromatograms(
        const OpenSwath::SpectrumAccessPtr input,
        const std::map< std::string, OpenSwath::ChromatogramPtr > & ms1_chromatograms,
        const std::vector< OpenSwath::SwathMap > swath_maps,
        OpenSwath::LightTargetedExperiment& transition_exp,
        const Param& feature_finder_param,
        TransformationDescription trafo,
        const double rt_extraction_window,
        FeatureMap& output,
        OpenSwathTSVWriter & tsv_writer,
        OpenSwathOSWWriter & osw_writer,
        bool ms1only = false);

    /** @brief Select which compounds to analyze in the next batch (and copy to output)
     *
     * This function will select which compounds to analyze in the next batch j
     * and will copy the corresponding compounds and transitions into the
     * output structure. The output will contain batch_size compounds.
     *
     * @param transition_exp_used_all input (all transitions for this swath)
     * @param transition_exp_used output (will contain only transitions for the next batch)
     * @param batch_size how many compounds per batch
     * @param j batch number (peptides from j*batch_size to j*batch_size+batch_size will be copied)
     *
     * @note The proteins will be copied completely without checking for a match
     *
    */
    void selectCompoundsForBatch_(const OpenSwath::LightTargetedExperiment& transition_exp_used_all,
      OpenSwath::LightTargetedExperiment& transition_exp_used, int batch_size, size_t j);

    /** @brief Copy the required transitions to output
     *
     * Copy all transitions matching to one of the compounds in the selected
     * peptide vector from all_transitions to the output.
     *
     * @param used_compounds Which peptides or metabolites to be used
     * @param all_transitions Transitions vector from which to select transitions
     * @param output Output vector containing matching transitions (taken from all_transitions)
     *
    */
    void copyBatchTransitions_(const std::vector<OpenSwath::LightCompound>& used_compounds,
      const std::vector<OpenSwath::LightTransition>& all_transitions,
      std::vector<OpenSwath::LightTransition>& output);

    /** @brief Function to prepare extraction coordinates that also correctly handles RT transformations
     *
     * Creates a set of (empty) chromatograms and extraction coordinates with
     * the correct ids, m/z and retention time start/end points to be extracted
     * by the ChromatogramExtractor.
     *
     * Handles rt extraction windows by calculating the correct transformation
     * for each coordinate.
     *
     * @param chrom_list Output of chromatograms (will be filled with empty chromatogram ptrs)
     * @param coordinates Output of extraction coordinates (will be filled with matching extraction coordinates)
     * @param transition_exp_used The transition experiment used to create the coordinates
     * @param ms1 Whether to perform MS1 (precursor ion) or MS2 (fragment ion) extraction
     * @param trafo_inverse Inverse transformation function
     * @param cp Parameter set for the chromatogram extraction
     *
    */
    void prepareExtractionCoordinates_(std::vector< OpenSwath::ChromatogramPtr > & chrom_list,
      std::vector< ChromatogramExtractorAlgorithm::ExtractionCoordinates > & coordinates,
      const OpenSwath::LightTargetedExperiment & transition_exp_used,
      const bool ms1, const TransformationDescription trafo_inverse,
      const ChromExtractParams & cp) const;

    /** @brief Simple function to prepare extraction coordinates
     *
     * This will take the targeted experiment and prepare extraction
     * coordinates (either MS1 or MS2) for extraction by the
     * ChromatogramExtractor.
     *
     * @param output_chromatograms Output of chromatograms (will be filled with empty chromatogram ptrs)
     * @param coordinates Output of extraction coordinates (will be filled with matching extraction coordinates)
     * @param transition_exp_used The transition experiment used to create the coordinates
     * @param rt_extraction_window Window for retention time extraction
     * @param ms1 Whether extraction coordinates should be created for MS1 (if false, it will be for MS2)
     *
    */
    void prepare_coordinates_sub(std::vector< OpenSwath::ChromatogramPtr > & output_chromatograms,
      std::vector< ChromatogramExtractorAlgorithm::ExtractionCoordinates > & coordinates,
      const OpenSwath::LightTargetedExperiment & transition_exp_used,
      const double rt_extraction_window, const bool ms1) const;

    /**
     * @brief Spectrum Access to the MS1 map (note that this is *not* threadsafe!)
     *
     * @note This pointer is not threadsafe, please use the lightClone() function to create a copy for each thread
     * @note This pointer may be NULL if use_ms1_traces_ is set to false
     *
     */
    OpenSwath::SpectrumAccessPtr ms1_map_;

    /// Whether to use the MS1 traces
    bool use_ms1_traces_;

  };

  /**
   * @brief Class to execute an OpenSwath Workflow for SONAR data
   *
   * The workflow will perform a complete OpenSWATH analysis, using scanning
   * SWATH data (SONAR data) instead of regular data. In this case, each
   * fragment ion may appear in multiple SWATH windows and thus needs to be
   * extracted from multiple maps.
   *
   */
  class OPENMS_DLLAPI OpenSwathWorkflowSonar :
    public OpenSwathWorkflow
  {

  public:
    explicit OpenSwathWorkflowSonar(bool use_ms1_traces) :
      OpenSwathWorkflow(use_ms1_traces)
    {}

    /** @brief Execute the OpenSWATH workflow on a set of SONAR SwathMaps and transitions.
     *
     * Executes the following operations on the given input:
     *
     * 1. Selecting the appropriate transitions for each SWATH window (using OpenSwathHelper::selectSwathTransitions)
     * 2. Extract the chromatograms from the SWATH maps (MS1 and MS2) using (ChromatogramExtractor)
     * 3. Pick peaks in the chromatograms and perform peak scoring (inside scoreAllChromatograms function)
     * 4. Write out chromatograms and found features
     *
     * Given that these are scanning SWATH maps, for each transition multiple
     * maps will be used for chromatogram extraction and scoring.
     *
     * @param swath_maps The raw data, expected to be scanning SWATH maps (SONAR)
     * @param trafo Transformation description (translating this runs' RT to normalized RT space)
     * @param cp Parameter set for the chromatogram extraction
     * @param feature_finder_param Parameter set for the feature finding in chromatographic dimension
     * @param transition_exp The set of assays to be extracted and scored
     * @param out_featureFile Output feature map to store identified features
     * @param store_features Whether features should be appended to the output feature map
     * @param tsv_writer TSV Writer object to store identified features in csv format
     * @param osw_writer OSW Writer object to store identified features in SQLite format
     * @param chromConsumer Chromatogram consumer object to store the extracted chromatograms
     * @param batchSize Size of the batches which should be extracted and scored
     * @param load_into_memory Whether to cache the current SONAR map(s) in memory
     *
    */
    void performExtractionSonar(const std::vector< OpenSwath::SwathMap > & swath_maps,
                                const TransformationDescription trafo,
                                const ChromExtractParams & cp,
                                const Param & feature_finder_param,
                                const OpenSwath::LightTargetedExperiment& transition_exp,
                                FeatureMap& out_featureFile,
                                bool store_features,
                                OpenSwathTSVWriter & tsv_writer,
                                OpenSwathOSWWriter & osw_writer,
                                Interfaces::IMSDataConsumer * chromConsumer,
                                int batchSize,
                                bool load_into_memory);

    /** @brief Compute start, end and total number of (virtual) SONAR windows
     *
    */
    void computeSonarWindows_(const std::vector< OpenSwath::SwathMap > & swath_maps,
                              double & sonar_winsize,
                              double & sonar_start,
                              double & sonar_end,
                              int & sonar_total_win);

    /** @brief Perform extraction from multiple SONAR windows
     *
    */
    void performSonarExtraction_(const std::vector< OpenSwath::SwathMap > & used_maps,
                                 const std::vector< ChromatogramExtractor::ExtractionCoordinates > & coordinates,
                                 std::vector< OpenSwath::ChromatogramPtr > & chrom_list,
                                 const ChromExtractParams & cp);

    /** @brief Add two chromatograms
     *
     * @param base_chrom The base chromatogram to which we will add intensity
     * @param newchrom The chromatogram to be added
     *
    */
    OpenSwath::ChromatogramPtr addChromatograms(OpenSwath::ChromatogramPtr base_chrom, OpenSwath::ChromatogramPtr newchrom);
  };

}

#endif // OPENMS_ANALYSIS_OPENSWATH_OPENSWATHWORKFLOW_H

