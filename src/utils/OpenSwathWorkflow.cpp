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

// Consumers
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataSqlConsumer.h>

// Files
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>
#include <OpenMS/FORMAT/SwathFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SwathWindowLoader.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathTSVWriter.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathOSWWriter.h>

// Kernel and implementations
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessTransforming.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMSInMemory.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/SwathMap.h>

// Helpers
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>

// Algorithms
#include <OpenMS/ANALYSIS/OPENSWATH/MRMRTNormalizer.h>
#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SwathMapMassCorrection.h>

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathWorkflow.h>

#include <cassert>
#include <limits>

// #define OPENSWATH_WORKFLOW_DEBUG

using namespace OpenMS;

static bool SortPairDoubleByFirst(const std::pair<double,double> & left, const std::pair<double,double> & right)
{
  return left.first < right.first;
}

// OpenMS base classes
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page UTILS_OpenSwathWorkflow OpenSwathWorkflow

  @brief Complete workflow to run OpenSWATH

  This implements the OpenSwath workflow as described in Rost and Rosenberger
  et al. (Nature Biotechnology, 2014) and provides a complete, integrated
  analysis tool without the need to run multiple tools consecutively.

  It executes the following steps in order:

  <ul>
    <li>Reading of input files, which can be provided as one single mzML or multiple "split" mzMLs (one per SWATH)</li>
    <li>Computing the retention time transformation using RT-normalization peptides</li>
    <li>Reading of the transition list</li>
    <li>Extracting the specified transitions</li>
    <li>Scoring the peak groups in the extracted ion chromatograms (XIC)</li>
    <li>Reporting the peak groups and the chromatograms</li>
  </ul>

  See below or have a look at the INI file (via "OpenSwathWorkflow -write_ini myini.ini") for available parameters and more functionality.

  <h3>Input: SWATH maps and transition list </h3>
  SWATH maps can be provided as mzML files, either as single file directly from
  the machine (this assumes that the SWATH method has 1 MS1 and then n MS2
  spectra which are ordered the same way for each cycle). E.g. a valid method
  would be MS1, MS2 [400-425], MS2 [425-450], MS1, MS2 [400-425], MS2 [425-450]
  while an invalid method would be MS1, MS2 [400-425], MS2 [425-450], MS1, MS2
  [425-450], MS2 [400-425] where MS2 [xx-yy] indicates an MS2 scan with an
  isolation window starting at xx and ending at yy. OpenSwathWorkflow will try
  to read the SWATH windows from the data, if this is not possible please
  provide a tab-separated list with the correct windows using the
  -swath_windows_file parameter (this is recommended). Note that the software
  expects extraction windows (e.g. which peptides to extract from
  which window) which cannot have overlaps, otherwise peptides will be
  extracted from two different windwos.

  Alternatively, a set of split files (n+1 mzML files) can be provided, each
  containing one SWATH map (or MS1 map).

  Since the file size can become rather large, it is recommended to not load the
  whole file into memory but rather cache it somewhere on the disk using a
  fast-access data format. This can be specified using the -readOptions cache
  parameter (this is recommended!).

  <h3>Parameters</h3>
  The current parameters are optimized for 2 hour gradients on SCIEX 5600 /
  6600 TripleTOF instruments with a peak width of around 30 seconds using iRT
  peptides.  If your chromatography differs, please consider adjusting
  -Scoring:TransitionGroupPicker:min_peak_width  to allow for smaller or larger
  peaks and adjust the -rt_extraction_window to use a different extraction
  window for the retention time. In m/z domain, it consider adjusting
  -mz_extraction_window to your instrument resolution, which can be in Th or
  ppm (using -ppm).

  Furthermore, if you wish to use MS1 information, use the -use_ms1_traces flag
  and provide an MS1 map in addition to the SWATH data.

  If you encounter issues with peak picking, try to disable peak filtering by
  setting -Scoring:TransitionGroupPicker:compute_peak_quality false which will
  disable the filtering of peaks by chromatographic quality. Furthermore, you
  can adjust the smoothing parameters for the peak picking, by adjusting
  -Scoring:TransitionGroupPicker:PeakPickerMRM:sgolay_frame_length or using a
  Gaussian smoothing based on your estimated peak width. Adjusting the signal
  to noise threshold will make the peaks wider or smaller.

  <h3>Output: Feature list and chromatograms </h3>
  The output of the OpenSwathWorkflow is a feature list, either as FeatureXML
  or as tsv (use -out_features or -out_tsv) while the latter is more memory
  friendly. If you analyze large dataset, it is recommended to only use
  -out_tsv and not -out_features. For downstream analysis (e.g. using mProphet or pyProphet)
  also the -out_tsv format is recommended.

  The feature list generated by -out_tsv is a tab-separated file. It can be
  used directly as input to the mProphet or pyProphet (a Python
  re-implementation of mProphet) software tool, see Reiter et al (2011, Nature
  Methods).

  In addition, the extracted chromatograms can be written out using the
  -out_chrom parameter.

  <h4> Feature list output format </h4>

  The tab-separated feature output contains the following information:

<CENTER>
  <table>
    <tr>
      <td ALIGN = "left" BGCOLOR="#EBEBEB"> Header row </td>
      <td ALIGN = "left" BGCOLOR="#EBEBEB"> Format </td>
      <td ALIGN = "left" BGCOLOR="#EBEBEB"> Description </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> transition_group_id </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> A unique id for the transition group (all chromatographic traces that are analyzed together)</td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> peptide_group_label </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> A unique id for the peptide group (will be the same for each charge state and heavy/light status) </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> run_id </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> An identifier for the run (currently always 0)</td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> filename </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> The input filename </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> RT </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Peak group retention time </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> id </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> A unique identifier for the peak group</td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Sequence </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Peptide sequence (no modifications) </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> FullPeptideName </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Full peptide sequence including modifications in Unimod format</td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Charge </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Int </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Assumed charge state</td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> m/z </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Precursor m/z</td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Intensity </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Peak group intensity (sum of all transitions)</td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> ProteinName </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Name of the associated protein</td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> decoy </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Whether the transition is decoy or not (0 = false, 1 = true) </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> assay_rt </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> The expected RT in seconds (based on normalized iRT value) </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> delta_rt </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> The difference between the expected RT and the peak group RT in seconds </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> leftWidth </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> The start of the peak group (left side) in seconds </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> main_var_xx_swath_prelim_score </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Initial score </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> norm_RT </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> The peak group retention time in normalized (iRT) space </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> nr_peaks </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Int </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> The number of transitions used </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> peak_apices_sum </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> The sum of all peak apices (may be used as alternative intensity) </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> potentialOutlier </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Potential outlier transitions (or "none" if none was detected)</td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> rightWidth </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> The end of the peak group (left side) in seconds </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> rt_score </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> The raw RT score (unnormalized) </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> sn_ratio </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> The raw S/N ratio </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> total_xic </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> The total XIC of the chromatogram </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> var_... </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> One of multiple sub-scores used by OpenSWATH to describe the peak group </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> aggr_prec_Peak_Area </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Intensity (peak area) of MS1 traces separated by semicolon </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> aggr_prec_Peak_Apex </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Intensity (peak apex) of MS1 traces separated by semicolon </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> aggr_prec_Fragment_Annotation </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Annotation of MS1 traces separated by semicolon </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> aggr_Peak_Area </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Intensity (peak area) of fragment ion traces separated by semicolon </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> aggr_Peak_Apex </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Intensity (peak apex) of fragment ion traces separated by semicolon </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> aggr_Fragment_Annotation </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=2> Annotation of fragment ion traces separated by semicolon </td>
    </tr>


  </table>
</CENTER>

  <B>The command line parameters of this tool are:</B>
  @verbinclude UTILS_OpenSwathWorkflow.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude UTILS_OpenSwathWorkflow.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
class TOPPOpenSwathWorkflow
  : public TOPPBase
{
public:

  TOPPOpenSwathWorkflow()
    : TOPPBase("OpenSwathWorkflow", "Complete workflow to run OpenSWATH", false)
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<files>", StringList(), "Input files separated by blank");
    setValidFormats_("in", ListUtils::create<String>("mzML,mzXML,sqMass"));

    registerInputFile_("tr", "<file>", "", "transition file ('TraML','tsv','pqp')");
    setValidFormats_("tr", ListUtils::create<String>("traML,tsv,pqp"));
    registerStringOption_("tr_type", "<type>", "", "input file type -- default: determined from file extension or content\n", false);
    setValidStrings_("tr_type", ListUtils::create<String>("traML,tsv,pqp"));

    // one of the following two needs to be set
    registerInputFile_("tr_irt", "<file>", "", "transition file ('TraML')", false);
    setValidFormats_("tr_irt", ListUtils::create<String>("traML"));

    registerInputFile_("rt_norm", "<file>", "", "RT normalization file (how to map the RTs of this run to the ones stored in the library). If set, tr_irt may be omitted.", false, true);
    setValidFormats_("rt_norm", ListUtils::create<String>("trafoXML"));

    registerInputFile_("swath_windows_file", "<file>", "", "Optional, tab separated file containing the SWATH windows for extraction: lower_offset upper_offset \\newline 400 425 \\newline ... Note that the first line is a header and will be skipped.", false, true);
    registerFlag_("sort_swath_maps", "Sort of input SWATH files when matching to SWATH windows from swath_windows_file", true);

    registerFlag_("use_ms1_traces", "Extract the precursor ion trace(s) and use for scoring", true);
    registerFlag_("enable_uis_scoring", "Enable additional scoring of identification assays", true);

    // one of the following two needs to be set
    registerOutputFile_("out_features", "<file>", "", "output file", false);
    setValidFormats_("out_features", ListUtils::create<String>("featureXML"));

    registerOutputFile_("out_tsv", "<file>", "", "TSV output file (mProphet compatible TSV file)", false);
    setValidFormats_("out_tsv", ListUtils::create<String>("tsv"));

    registerOutputFile_("out_osw", "<file>", "", "OSW output file (PyProphet compatible SQLite file)", false);
    setValidFormats_("out_osw", ListUtils::create<String>("osw"));

    registerOutputFile_("out_chrom", "<file>", "", "Also output all computed chromatograms output in mzML (chrom.mzML) or sqMass (SQLite format)", false, true);
    setValidFormats_("out_chrom", ListUtils::create<String>("mzML,sqMass"));

    registerDoubleOption_("min_upper_edge_dist", "<double>", 0.0, "Minimal distance to the edge to still consider a precursor, in Thomson", false, true);
    registerDoubleOption_("rt_extraction_window", "<double>", 600.0, "Only extract RT around this value (-1 means extract over the whole range, a value of 600 means to extract around +/- 300 s of the expected elution).", false);
    registerDoubleOption_("extra_rt_extraction_window", "<double>", 0.0, "Output an XIC with a RT-window that by this much larger (e.g. to visually inspect a larger area of the chromatogram)", false, true);
    registerDoubleOption_("mz_extraction_window", "<double>", 0.05, "Extraction window used (in Thomson, to use ppm see -ppm flag)", false);
    setMinFloat_("mz_extraction_window", 0.0);
    setMinFloat_("extra_rt_extraction_window", 0.0);
    registerFlag_("ppm", "m/z extraction_window is in ppm");
    registerFlag_("sonar", "data is scanning SWATH data");

    registerDoubleOption_("min_rsq", "<double>", 0.95, "Minimum r-squared of RT peptides regression", false, true);
    registerDoubleOption_("min_coverage", "<double>", 0.6, "Minimum relative amount of RT peptides to keep", false, true);

    registerFlag_("split_file_input", "The input files each contain one single SWATH (alternatively: all SWATH are in separate files)", true);
    registerFlag_("use_elution_model_score", "Turn on elution model score (EMG fit to peak)", true);

    registerStringOption_("readOptions", "<name>", "normal", "Whether to run OpenSWATH directly on the input data, cache data to disk first or to perform a datareduction step first. If you choose cache, make sure to also set tempDirectory", false, true);
    setValidStrings_("readOptions", ListUtils::create<String>("normal,cache,cacheWorkingInMemory,workingInMemory"));

    registerStringOption_("mz_correction_function", "<name>", "none", "Use the retention time normalization peptide MS2 masses to perform a mass correction (linear, weighted by intensity linear or quadratic) of all spectra.", false, true);
    setValidStrings_("mz_correction_function", ListUtils::create<String>("none,unweighted_regression,weighted_regression,quadratic_regression,weighted_quadratic_regression,weighted_quadratic_regression_delta_ppm,quadratic_regression_delta_ppm"));
    registerDoubleOption_("irt_mz_extraction_window", "<double>", 0.05, "Extraction window used for iRT and m/z correction (in Thomson, use ppm use -ppm flag)", false, true);
    registerFlag_("ppm_irtwindow", "iRT m/z extraction_window is in ppm", true);

    // TODO terminal slash !
    registerStringOption_("tempDirectory", "<tmp>", "/tmp/", "Temporary directory to store cached files for example", false, true);

    registerStringOption_("extraction_function", "<name>", "tophat", "Function used to extract the signal", false, true);
    setValidStrings_("extraction_function", ListUtils::create<String>("tophat,bartlett"));

    registerIntOption_("batchSize", "<number>", 0, "The batch size of chromatograms to process (0 means to only have one batch, sensible values are around 500-1000)", false, true);
    setMinInt_("batchSize", 0);

    registerSubsection_("Scoring", "Scoring parameters section");
    registerSubsection_("Library", "Library parameters section");

    registerSubsection_("RTNormalization", "Parameters for the RTNormalization for iRT petides. This specifies how the RT alignment is performed and how outlier detection is applied. Outlier detection can be done iteratively (by default) which removes one outlier per iteration or using the RANSAC algorithm.");
    registerSubsection_("Debugging", "Debugging");
  }

  Param getSubsectionDefaults_(const String& name) const override
  {
    if (name == "Scoring")
    {
      // set sensible default parameters
      Param feature_finder_param = MRMFeatureFinderScoring().getDefaults();
      feature_finder_param.remove("rt_extraction_window");
      feature_finder_param.setValue("rt_normalization_factor", 100.0); // for iRT peptides between 0 and 100 (more or less)

      feature_finder_param.setValue("TransitionGroupPicker:min_peak_width", 14.0);
      feature_finder_param.setValue("TransitionGroupPicker:recalculate_peaks", "true");
      feature_finder_param.setValue("TransitionGroupPicker:compute_peak_quality", "true");
      feature_finder_param.setValue("TransitionGroupPicker:minimal_quality", -1.5);
      feature_finder_param.setValue("TransitionGroupPicker:background_subtraction", "none");
      feature_finder_param.remove("TransitionGroupPicker:stop_after_intensity_ratio");

      // Peak Picker
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:use_gauss", "false");
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:sgolay_polynomial_order", 3);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:sgolay_frame_length", 11);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:peak_width", -1.0);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:remove_overlapping_peaks", "true");
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:write_sn_log_messages", "false"); // no log messages
      // TODO it seems that the legacy method produces slightly larger peaks, e.g. it will not cut off peaks too early
      // however the same can be achieved by using a relatively low SN cutoff in the -Scoring:TransitionGroupPicker:PeakPickerMRM:signal_to_noise 0.5
      feature_finder_param.setValue("TransitionGroupPicker:recalculate_peaks_max_z", 0.75);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:method", "corrected");
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:signal_to_noise", 0.1);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:gauss_width", 30.0);
      feature_finder_param.setValue("uis_threshold_sn",0);
      feature_finder_param.setValue("uis_threshold_peak_area",0);
      feature_finder_param.remove("TransitionGroupPicker:PeakPickerMRM:sn_win_len");
      feature_finder_param.remove("TransitionGroupPicker:PeakPickerMRM:sn_bin_count");
      feature_finder_param.remove("TransitionGroupPicker:PeakPickerMRM:stop_after_feature");

      // EMG Scoring - turn off by default since it is very CPU-intensive
      feature_finder_param.remove("Scores:use_elution_model_score");
      feature_finder_param.setValue("EMGScoring:max_iteration", 10);
      feature_finder_param.remove("EMGScoring:interpolation_step");
      feature_finder_param.remove("EMGScoring:tolerance_stdev_bounding_box");
      feature_finder_param.remove("EMGScoring:deltaAbsError");

      // remove these parameters
      feature_finder_param.remove("add_up_spectra");
      feature_finder_param.remove("spacing_for_spectra_resampling");
      feature_finder_param.remove("EMGScoring:statistics:mean");
      feature_finder_param.remove("EMGScoring:statistics:variance");
      return feature_finder_param;
    }
    else if (name == "RTNormalization")
    {
      Param p;

      p.setValue("alignmentMethod", "linear", "How to perform the alignment to the normalized RT space using anchor points. 'linear': perform linear regression (for few anchor points). 'interpolated': Interpolate between anchor points (for few, noise-free anchor points). 'lowess' Use local regression (for many, noisy anchor points). 'b_spline' use b splines for smoothing.");
      p.setValidStrings("alignmentMethod", ListUtils::create<String>("linear,interpolated,lowess,b_spline"));
      p.setValue("lowess:span", 2.0/3, "Span parameter for lowess");
      p.setMinFloat("lowess:span", 0.0);
      p.setMaxFloat("lowess:span", 1.0);
      p.setValue("b_spline:num_nodes", 5, "Number of nodes for b spline");
      p.setMinInt("b_spline:num_nodes", 0);

      p.setValue("outlierMethod", "iter_residual", "Which outlier detection method to use (valid: 'iter_residual', 'iter_jackknife', 'ransac', 'none'). Iterative methods remove one outlier at a time. Jackknife approach optimizes for maximum r-squared improvement while 'iter_residual' removes the datapoint with the largest residual error (removal by residual is computationally cheaper, use this with lots of peptides).");
      p.setValidStrings("outlierMethod", ListUtils::create<String>("iter_residual,iter_jackknife,ransac,none"));

      p.setValue("useIterativeChauvenet", "false", "Whether to use Chauvenet's criterion when using iterative methods. This should be used if the algorithm removes too many datapoints but it may lead to true outliers being retained.");
      p.setValidStrings("useIterativeChauvenet", ListUtils::create<String>("true,false"));

      p.setValue("RANSACMaxIterations", 1000, "Maximum iterations for the RANSAC outlier detection algorithm.");
      p.setValue("RANSACMaxPercentRTThreshold", 3, "Maximum threshold in RT dimension for the RANSAC outlier detection algorithm (in percent of the total gradient). Default is set to 3% which is around +/- 4 minutes on a 120 gradient.");
      p.setValue("RANSACSamplingSize", 10, "Sampling size of data points per iteration for the RANSAC outlier detection algorithm.");

      p.setValue("estimateBestPeptides", "false", "Whether the algorithms should try to choose the best peptides based on their peak shape for normalization. Use this option you do not expect all your peptides to be detected in a sample and too many 'bad' peptides enter the outlier removal step (e.g. due to them being endogenous peptides or using a less curated list of peptides).");
      p.setValidStrings("estimateBestPeptides", ListUtils::create<String>("true,false"));

      p.setValue("InitialQualityCutoff", 0.5, "The initial overall quality cutoff for a peak to be scored (range ca. -2 to 2)");
      p.setValue("OverallQualityCutoff", 5.5, "The overall quality cutoff for a peak to go into the retention time estimation (range ca. 0 to 10)");
      p.setValue("NrRTBins", 10, "Number of RT bins to use to compute coverage. This option should be used to ensure that there is a complete coverage of the RT space (this should detect cases where only a part of the RT gradient is actually covered by normalization peptides)");
      p.setValue("MinPeptidesPerBin", 1, "Minimal number of peptides that are required for a bin to counted as 'covered'");
      p.setValue("MinBinsFilled", 8, "Minimal number of bins required to be covered");
      return p;
    }
    else if (name == "Debugging")
    {
      Param p;
      p.setValue("irt_mzml", "", "Chromatogram mzML containing the iRT peptides");
      // p.setValidFormats_("irt_mzml", ListUtils::create<String>("mzML"));
      p.setValue("irt_trafo", "", "Transformation file for RT transform");
      // p.setValidFormats_("irt_trafo", ListUtils::create<String>("trafoXML"));
      return p;
    }
    else if (name == "Library")
    {
      return TransitionTSVFile().getDefaults();
    }
    else
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown subsection", name);
    }
  }

  void loadSwathFiles(StringList& file_list, bool split_file, String tmp, String readoptions,
    boost::shared_ptr<ExperimentalSettings > & exp_meta,
    std::vector< OpenSwath::SwathMap > & swath_maps)
  {
    SwathFile swath_file;
    swath_file.setLogType(log_type_);

    if (split_file || file_list.size() > 1)
    {
      // TODO cannot use data reduction here any more ...
      swath_maps = swath_file.loadSplit(file_list, tmp, exp_meta, readoptions);
    }
    else
    {
      FileTypes::Type in_file_type = FileTypes::nameToType(file_list[0]);
      if (in_file_type == FileTypes::MZML || file_list[0].suffix(4).toLower() == "mzml"
        || file_list[0].suffix(7).toLower() == "mzml.gz"  )
      {
        swath_maps = swath_file.loadMzML(file_list[0], tmp, exp_meta, readoptions);
      }
      else if (in_file_type == FileTypes::MZXML || file_list[0].suffix(5).toLower() == "mzxml"
        || file_list[0].suffix(8).toLower() == "mzxml.gz"  )
      {
        swath_maps = swath_file.loadMzXML(file_list[0], tmp, exp_meta, readoptions);
      }
      else if (in_file_type == FileTypes::SQMASS || file_list[0].suffix(6).toLower() == "sqmass")
      {
        swath_maps = swath_file.loadSqMass(file_list[0], exp_meta);
      }
      else
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "Input file needs to have ending mzML or mzXML");
      }
    }
  }

  /**
   * @brief Load the retention time transformation file
   *
   * This function will create the retention time transformation either by
   * loading a provided .trafoXML file or determine it from the data itself by
   * extracting the transitions specified in the irt_tr_file TraML file.
   *
   * @param trafo_in Input trafoXML file (if not empty, transformation will be
   *                 loaded from this file)
   * @param irt_tr_file  Input TraML file containing transitions (if trafo_in
   *                     is empty, this file will be loaded and transitions
   *                     will be extracted)
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
  TransformationDescription loadTrafoFile(String trafo_in, String irt_tr_file,
    std::vector< OpenSwath::SwathMap > & swath_maps, double min_rsq, double min_coverage,
    const Param& feature_finder_param, const ChromExtractParams& cp_irt,
    const Param& irt_detection_param, const String & mz_correction_function,
    Size debug_level, bool sonar, bool load_into_memory, const String& debug_output)
  {
    TransformationDescription trafo_rtnorm;
    if (!trafo_in.empty())
    {
      // get read RT normalization file
      TransformationXMLFile trafoxml;
      trafoxml.load(trafo_in, trafo_rtnorm, false);
      Param model_params = getParam_().copy("model:", true);
      model_params.setValue("symmetric_regression", "false");
      model_params.setValue("span", irt_detection_param.getValue("lowess:span"));
      model_params.setValue("num_nodes", irt_detection_param.getValue("b_spline:num_nodes"));
      String model_type = irt_detection_param.getValue("alignmentMethod");
      trafo_rtnorm.fitModel(model_type, model_params);
    }
    else if (!irt_tr_file.empty())
    {
      // Loading iRT file
      std::cout << "Will load iRT transitions and try to find iRT peptides" << std::endl;
      TraMLFile traml;
      OpenMS::TargetedExperiment irt_transitions;
      traml.load(irt_tr_file, irt_transitions);

      // perform extraction
      OpenSwathRetentionTimeNormalization wf;
      wf.setLogType(log_type_);
      trafo_rtnorm = wf.performRTNormalization(irt_transitions, swath_maps, min_rsq, min_coverage,
          feature_finder_param, cp_irt, irt_detection_param, mz_correction_function, debug_level, sonar, load_into_memory);

      if (!debug_output.empty())
      {
        TransformationXMLFile trafoxml;
        trafoxml.store(debug_output, trafo_rtnorm);
      }
    }
    return trafo_rtnorm;
  }

  ExitCodes main_(int, const char **) override
  {
    ///////////////////////////////////
    // Prepare Parameters
    ///////////////////////////////////
    StringList file_list = getStringList_("in");
    String tr_file = getStringOption_("tr");

    Param irt_detection_param = getParam_().copy("RTNormalization:", true);

    //tr_file input file type
    FileHandler fh_tr_type;
    FileTypes::Type tr_type = FileTypes::nameToType(getStringOption_("tr_type"));

    if (tr_type == FileTypes::UNKNOWN)
    {
      tr_type = fh_tr_type.getType(tr_file);
      writeDebug_(String("Input file type: ") + FileTypes::typeToName(tr_type), 2);
    }

    if (tr_type == FileTypes::UNKNOWN)
    {
      writeLog_("Error: Could not determine input file type!");
      return PARSE_ERROR;
    }

    String out = getStringOption_("out_features");
    String out_tsv = getStringOption_("out_tsv");
    String out_osw = getStringOption_("out_osw");

    String irt_tr_file = getStringOption_("tr_irt");
    String trafo_in = getStringOption_("rt_norm");

    String out_chrom = getStringOption_("out_chrom");
    bool ppm = getFlag_("ppm");
    bool irt_ppm = getFlag_("ppm_irtwindow");
    bool split_file = getFlag_("split_file_input");
    bool use_emg_score = getFlag_("use_elution_model_score");
    bool force = getFlag_("force");
    bool sonar = getFlag_("sonar");
    bool sort_swath_maps = getFlag_("sort_swath_maps");
    bool use_ms1_traces = getFlag_("use_ms1_traces");
    bool enable_uis_scoring = getFlag_("enable_uis_scoring");
    double min_upper_edge_dist = getDoubleOption_("min_upper_edge_dist");
    double mz_extraction_window = getDoubleOption_("mz_extraction_window");
    double irt_mz_extraction_window = getDoubleOption_("irt_mz_extraction_window");
    double rt_extraction_window = getDoubleOption_("rt_extraction_window");
    double extra_rt_extract = getDoubleOption_("extra_rt_extraction_window");
    String extraction_function = getStringOption_("extraction_function");
    String swath_windows_file = getStringOption_("swath_windows_file");
    int batchSize = (int)getIntOption_("batchSize");
    Size debug_level = (Size)getIntOption_("debug");

    double min_rsq = getDoubleOption_("min_rsq");
    double min_coverage = getDoubleOption_("min_coverage");

    Param debug_params = getParam_().copy("Debugging:", true);

    String readoptions = getStringOption_("readOptions");
    String mz_correction_function = getStringOption_("mz_correction_function");
    String tmp = getStringOption_("tempDirectory");

    ///////////////////////////////////
    // Parameter validation
    ///////////////////////////////////

    bool load_into_memory = false;
    bool is_sqmass_input  = (file_list[0].suffix(6).toLower() == "sqmass");
    if (readoptions == "cacheWorkingInMemory")
    {
      readoptions = "cache";
      load_into_memory = true;
    }
    else if (readoptions == "workingInMemory")
    {
      readoptions = "normal";
      load_into_memory = true;
    }

    if (is_sqmass_input && !load_into_memory)
    {
      std::cout << "When using sqMass input files, it is highly recommended to use the workingInMemory option as otherwise data access will be very slow." << std::endl;
    }

    if (trafo_in.empty() && irt_tr_file.empty())
    {
      std::cout << "Since neither rt_norm nor tr_irt is set, OpenSWATH will " <<
        "not use RT-transformation (rather a null transformation will be applied)" << std::endl;
    }
    if ( (out.empty() && out_tsv.empty() && out_osw.empty()) || (!out.empty() && !out_tsv.empty())  || (!out.empty() && !out_osw.empty())  || (!out_tsv.empty() && !out_osw.empty()) )
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Either out_features, out_tsv or out_osw needs to be set (but not two or three at the same time)");
    }
    if (!out_osw.empty() && tr_type != FileTypes::PQP)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "OSW output files can only be generated in combination with PQP input files (-tr).");
    }

    // Check swath window input
    if (!swath_windows_file.empty())
    {
      LOG_INFO << "Validate provided Swath windows file:" << std::endl;
      std::vector<double> swath_prec_lower;
      std::vector<double> swath_prec_upper;
      SwathWindowLoader::readSwathWindows(swath_windows_file, swath_prec_lower, swath_prec_upper);

      LOG_INFO << "Read Swath maps file with " << swath_prec_lower.size() << " windows." << std::endl;
      for (Size i = 0; i < swath_prec_lower.size(); i++)
      {
        LOG_DEBUG << "Read lower swath window " << swath_prec_lower[i] << " and upper window " << swath_prec_upper[i] << std::endl;
      }
    }

    ChromExtractParams cp;
    cp.min_upper_edge_dist   = min_upper_edge_dist;
    cp.mz_extraction_window  = mz_extraction_window;
    cp.ppm                   = ppm;
    cp.rt_extraction_window  = rt_extraction_window,
    cp.extraction_function   = extraction_function;
    cp.extra_rt_extract      = extra_rt_extract;

    ChromExtractParams cp_irt = cp;
    cp_irt.rt_extraction_window = -1; // extract the whole RT range
    cp_irt.mz_extraction_window = irt_mz_extraction_window;
    cp_irt.ppm                  = irt_ppm;

    Param feature_finder_param = getParam_().copy("Scoring:", true);
    Param tsv_reader_param = getParam_().copy("Library:", true);
    if (use_emg_score)
    {
      feature_finder_param.setValue("Scores:use_elution_model_score", "true");
    }
    else
    {
      feature_finder_param.setValue("Scores:use_elution_model_score", "false");
    }
    if (use_ms1_traces)
    {
      feature_finder_param.setValue("Scores:use_ms1_correlation", "true");
      feature_finder_param.setValue("Scores:use_ms1_fullscan", "true");
    }
    if (enable_uis_scoring)
    {
      feature_finder_param.setValue("Scores:use_uis_scores", "true");
    }

    ///////////////////////////////////
    // Load the transitions
    ///////////////////////////////////
    OpenSwath::LightTargetedExperiment transition_exp;
    ProgressLogger progresslogger;
    progresslogger.setLogType(log_type_);
    if (tr_type == FileTypes::TRAML || tr_file.suffix(5).toLower() == "traml"  )
    {
      progresslogger.startProgress(0, 1, "Load TraML file");
      TargetedExperiment targeted_exp;
      TraMLFile().load(tr_file, targeted_exp);
      OpenSwathDataAccessHelper::convertTargetedExp(targeted_exp, transition_exp);
      progresslogger.endProgress();
    }
    else if (tr_type == FileTypes::PQP || tr_file.suffix(3).toLower() == "pqp"  )
    {
      progresslogger.startProgress(0, 1, "Load PQP file");
      TransitionPQPFile().convertPQPToTargetedExperiment(tr_file.c_str(), transition_exp);
      progresslogger.endProgress();

      remove(out_osw.c_str());
      if (!out_osw.empty())
      {
        std::ifstream  src(tr_file.c_str(), std::ios::binary);
        std::ofstream  dst(out_osw.c_str(), std::ios::binary);

        dst << src.rdbuf();
      }
    }
    else if (tr_type == FileTypes::TSV || tr_file.suffix(3).toLower() == "tsv"  )
    {
      progresslogger.startProgress(0, 1, "Load TSV file");
      TransitionTSVFile tsv_reader;
      tsv_reader.setParameters(tsv_reader_param);
      tsv_reader.convertTSVToTargetedExperiment(tr_file.c_str(), tr_type, transition_exp);
      progresslogger.endProgress();
    }
    else
    {
      LOG_ERROR << "Provide valid TraML, TSV or PQP transition file." << std::endl;
      return PARSE_ERROR;
    }
    LOG_INFO << "Loaded " << transition_exp.getProteins().size() << " proteins, " << 
      transition_exp.getCompounds().size() << " compounds with " << transition_exp.getTransitions().size() << " transitions." << std::endl;


    ///////////////////////////////////
    // Load the SWATH files
    ///////////////////////////////////

    // (i) Load files
    boost::shared_ptr<ExperimentalSettings> exp_meta(new ExperimentalSettings);
    std::vector< OpenSwath::SwathMap > swath_maps;
    loadSwathFiles(file_list, split_file, tmp, readoptions, exp_meta, swath_maps);

    // (ii) Allow the user to specify the SWATH windows
    if (!swath_windows_file.empty())
    {
      SwathWindowLoader::annotateSwathMapsFromFile(swath_windows_file, swath_maps, sort_swath_maps);
    }

    for (Size i = 0; i < swath_maps.size(); i++)
    {
      LOG_DEBUG << "Found swath map " << i 
        << " with lower " << swath_maps[i].lower
        << " and upper " << swath_maps[i].upper 
        << " and " << swath_maps[i].sptr->getNrSpectra()
        << " spectra." << std::endl;
    }

    // (iii) Sanity check: there should be no overlap between the windows:
    std::vector<std::pair<double, double> > sw_windows;
    for (Size i = 0; i < swath_maps.size(); i++)
    {
      if (!swath_maps[i].ms1)
      {
        sw_windows.push_back(std::make_pair(swath_maps[i].lower, swath_maps[i].upper));
      }
    }
    std::sort(sw_windows.begin(), sw_windows.end(), SortPairDoubleByFirst);

    for (Size i = 1; i < sw_windows.size(); i++)
    {
      double lower_map_end = sw_windows[i-1].second - min_upper_edge_dist;
      double upper_map_start = sw_windows[i].first;
      LOG_DEBUG << "Extraction will go up to " << lower_map_end << " and continue at " << upper_map_start << std::endl;

      if (upper_map_start - lower_map_end > 0.01)
      {
        LOG_WARN << "Extraction will have a gap between " << lower_map_end << " and " << upper_map_start << std::endl;
        if (!force)
        {
          LOG_ERROR << "Extraction windows have a gap. Will abort (override with -force)" << std::endl;
          return PARSE_ERROR;
        }
      }

      if (sonar) {continue;} // skip next step as expect them to overlap ...

      if (lower_map_end - upper_map_start > 0.01)
      {
        LOG_WARN << "Extraction will overlap between " << lower_map_end << " and " << upper_map_start << std::endl;
        LOG_WARN << "This will lead to multiple extraction of the transitions in the overlapping region" <<
                    "which will lead to duplicated output. It is very unlikely that you want this." << std::endl;
        LOG_WARN << "Please fix this by providing an appropriate extraction file with -swath_windows_file" << std::endl;
        if (!force)
        {
          LOG_ERROR << "Extraction windows overlap. Will abort (override with -force)" << std::endl;
          return PARSE_ERROR;
        }
      }

    }

    ///////////////////////////////////
    // Get the transformation information (using iRT peptides)
    ///////////////////////////////////
    String debug_output = debug_params.getValue("irt_trafo");
    TransformationDescription trafo_rtnorm = loadTrafoFile(trafo_in,
        irt_tr_file, swath_maps, min_rsq, min_coverage, feature_finder_param,
        cp_irt, irt_detection_param, mz_correction_function, debug_level,
        sonar, load_into_memory, debug_output);

    ///////////////////////////////////
    // Set up chromatogram output
    // Either use chrom.mzML or sqlite DB
    ///////////////////////////////////
    Interfaces::IMSDataConsumer * chromatogramConsumer;
    if (!out_chrom.empty())
    {
      if (out_chrom.hasSuffix(".sqMass"))
      {
        bool full_meta = false; // can lead to very large files in memory
        bool lossy_compression = true;
        chromatogramConsumer = new MSDataSqlConsumer(out_chrom, 500, full_meta, lossy_compression);
      }
      else
      {
        PlainMSDataWritingConsumer * chromConsumer = new PlainMSDataWritingConsumer(out_chrom);
        int expected_chromatograms = transition_exp.transitions.size();
        chromConsumer->setExpectedSize(0, expected_chromatograms);
        chromConsumer->setExperimentalSettings(*exp_meta);
        chromConsumer->getOptions().setWriteIndex(true);  // ensure that we write the index
        chromConsumer->addDataProcessing(getProcessingInfo_(DataProcessing::SMOOTHING));

        // prepare data structures for lossy compression
        MSNumpressCoder::NumpressConfig npconfig_mz;
        MSNumpressCoder::NumpressConfig npconfig_int;
        npconfig_mz.estimate_fixed_point = true; // critical
        npconfig_int.estimate_fixed_point = true; // critical
        npconfig_mz.numpressErrorTolerance = -1.0; // skip check, faster
        npconfig_int.numpressErrorTolerance = -1.0; // skip check, faster
        npconfig_mz.setCompression("linear");
        npconfig_int.setCompression("slof");
        npconfig_mz.linear_fp_mass_acc = 0.05; // set the desired RT accuracy in seconds

        chromConsumer->getOptions().setNumpressConfigurationMassTime(npconfig_mz);
        chromConsumer->getOptions().setNumpressConfigurationIntensity(npconfig_int);
        chromConsumer->getOptions().setCompression(true);

        chromatogramConsumer = chromConsumer;
      }
    }
    else
    {
      chromatogramConsumer = new NoopMSDataWritingConsumer("");
    }

    ///////////////////////////////////
    // Extract and score
    ///////////////////////////////////
    FeatureMap out_featureFile;

    OpenSwathTSVWriter tsvwriter(out_tsv, file_list[0], use_ms1_traces, sonar, enable_uis_scoring); // only active if filename not empty
    OpenSwathOSWWriter oswwriter(out_osw, file_list[0], use_ms1_traces, sonar, enable_uis_scoring); // only active if filename not empty

    if (sonar)
    {
      OpenSwathWorkflowSonar wf(use_ms1_traces);
      wf.setLogType(log_type_);
      wf.performExtractionSonar(swath_maps, trafo_rtnorm, cp, feature_finder_param, transition_exp,
          out_featureFile, !out.empty(), tsvwriter, oswwriter, chromatogramConsumer, batchSize, load_into_memory);
    }
    else
    {
      OpenSwathWorkflow wf(use_ms1_traces);
      wf.setLogType(log_type_);
      wf.performExtraction(swath_maps, trafo_rtnorm, cp, feature_finder_param, transition_exp,
          out_featureFile, !out.empty(), tsvwriter, oswwriter, chromatogramConsumer, batchSize, load_into_memory);
    }

    if (!out.empty())
    {
      addDataProcessing_(out_featureFile, getProcessingInfo_(DataProcessing::QUANTITATION));
      out_featureFile.ensureUniqueId();
      FeatureXMLFile().store(out, out_featureFile);
    }

    delete chromatogramConsumer;

    return EXECUTION_OK;
  }

};

int main(int argc, const char ** argv)
{
  TOPPOpenSwathWorkflow tool;
  return tool.main(argc, argv);
}

/// @endcond
