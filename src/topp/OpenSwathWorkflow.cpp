// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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
#include <OpenMS/FORMAT/SwathFile.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataTransformingConsumer.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SwathWindowLoader.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SwathQC.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathTSVWriter.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathOSWWriter.h>
#include <OpenMS/SYSTEM/File.h>

// Kernel and implementations
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessTransforming.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMSInMemory.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>

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

// OpenMS base classes
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/APPLICATIONS/OpenSwathBase.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>


#include <QDir>

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_OpenSwathWorkflow OpenSwathWorkflow

@brief Complete workflow to run OpenSWATH

This implements the OpenSWATH workflow as described in Rost and Rosenberger
et al. (Nature Biotechnology, 2014) and provides a complete, integrated
analysis tool without the need to run multiple tools consecutively. See also
http://openswath.org/ for additional documentation.

It executes the following steps in order, which is implemented in @ref OpenMS::OpenSwathWorkflow "OpenSwathWorkflow":

<ul>
  <li>Reading of input files, which can be provided as one single mzML or multiple "split" mzMLs (one per SWATH)</li>
  <li>Computing the retention time transformation using RT-normalization peptides</li>
  <li>Reading of the transition list</li>
  <li>Extracting the specified transitions</li>
  <li>Scoring the peak groups in the extracted ion chromatograms (XIC)</li>
  <li>Reporting the peak groups and the chromatograms</li>
</ul>


See below or have a look at the INI file (via "OpenSwathWorkflow -write_ini myini.ini") for available parameters and more functionality.

<h3>Input: SWATH maps and assay library (transition list) </h3>
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
extracted from two different windows.

Alternatively, a set of split files (n+1 mzML files) can be provided, each
containing one SWATH map (or MS1 map).

Since the file size can become rather large, it is recommended to not load the
whole file into memory but rather cache it somewhere on the disk using a
fast-access data format. This can be specified using the -readOptions cache
parameter (this is recommended!).

The assay library (transition list) is provided through the @p -tr parameter and can be in one of the following formats:

  <ul>
    <li> @ref OpenMS::TraMLFile "TraML" </li>
    <li> @ref OpenMS::TransitionTSVFile "OpenSWATH TSV transition lists" </li>
    <li> @ref OpenMS::TransitionPQPFile "OpenSWATH PQP SQLite files" </li>
    <li> SpectraST MRM transition lists </li>
    <li> Skyline transition lists </li>
    <li> Spectronaut transition lists </li>
  </ul>

<h3>Parameters</h3>
The current parameters are optimized for 2 hour gradients on SCIEX 5600 /
6600 TripleTOF instruments with a peak width of around 30 seconds using iRT
peptides.  If your chromatography differs, please consider adjusting
@p -Scoring:TransitionGroupPicker:min_peak_width  to allow for smaller or larger
peaks and adjust the @p -rt_extraction_window to use a different extraction
window for the retention time. In m/z domain, consider adjusting
@p -mz_extraction_window to your instrument resolution, which can be in Th or
ppm.

Furthermore, if you wish to use MS1 information, use the @p -use_ms1_traces flag
and provide an MS1 map in addition to the SWATH data.

If you encounter issues with peak picking, try to disable peak filtering by
setting @p -Scoring:TransitionGroupPicker:compute_peak_quality false which will
disable the filtering of peaks by chromatographic quality. Furthermore, you
can adjust the smoothing parameters for the peak picking, by adjusting
@p -Scoring:TransitionGroupPicker:PeakPickerChromatogram:sgolay_frame_length or using a
Gaussian smoothing based on your estimated peak width. Adjusting the signal
to noise threshold will make the peaks wider or smaller.

<h3>Output: Feature list and chromatograms </h3>
The output of the OpenSwathWorkflow is a feature list, either as FeatureXML
or as tsv (use @p -out_features or @p -out_tsv) while the latter is more memory
friendly and can be directly used as input to other tools such as mProphet or
pyProphet. If you analyze large datasets, it is recommended to only use @p
-out_tsv and not @p -out_features. For downstream analysis (e.g. using mProphet
 or pyProphet) also the @p -out_tsv format is recommended.

The feature list generated by @p -out_tsv is a tab-separated file. It can be
used directly as input to the mProphet or pyProphet (a Python
re-implementation of mProphet) software tool, see Reiter et al (2011, Nature
Methods).

In addition, the extracted chromatograms can be written out using the
@p -out_chrom parameter.

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
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> MC </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Int </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Missed cleavages of the sequence (assuming Trypsin as protease) </td>
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
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> masserror_ppm </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float List </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Pairs of fragment masses (m/z) and their associated error in ppm for all transitions</td>
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
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Annotation of fragment ion traces separated by semicolon </td>
    </tr>
  </table>
</CENTER>

<h3>Execution flow:</h3>

The overall execution flow for this tool is implemented in @ref OpenMS::OpenSwathWorkflow "OpenSwathWorkflow".

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_OpenSwathWorkflow.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_OpenSwathWorkflow.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
class TOPPOpenSwathWorkflow
  : public TOPPOpenSwathBase
{
public:

  TOPPOpenSwathWorkflow()
    : TOPPOpenSwathBase("OpenSwathWorkflow", "Complete workflow to run OpenSWATH", true)
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
    setValidFormats_("tr_irt", ListUtils::create<String>("traML,tsv,pqp"));

    // one of the following two needs to be set
    registerInputFile_("tr_irt_nonlinear", "<file>", "", "additional nonlinear transition file ('TraML')", false);
    setValidFormats_("tr_irt_nonlinear", ListUtils::create<String>("traML,tsv,pqp"));

    registerInputFile_("rt_norm", "<file>", "", "RT normalization file (how to map the RTs of this run to the ones stored in the library). If set, tr_irt may be omitted.", false, true);
    setValidFormats_("rt_norm", ListUtils::create<String>("trafoXML"));

    registerInputFile_("swath_windows_file", "<file>", "", "Optional, tab-separated file containing the SWATH windows for extraction: lower_offset upper_offset. Note that the first line is a header and will be skipped.", false);
    registerFlag_("sort_swath_maps", "Sort input SWATH files when matching to SWATH windows from swath_windows_file", true);

    registerStringOption_("enable_ms1", "<name>", "true", "Extract the precursor ion trace(s) and use for scoring if present", false, true);
    setValidStrings_("enable_ms1", ListUtils::create<String>("true,false"));

    registerStringOption_("enable_ipf", "<name>", "true", "Enable additional scoring of identification assays using IPF (see online documentation)", false, true);
    setValidStrings_("enable_ipf", ListUtils::create<String>("true,false"));

    // one of the following two needs to be set
    registerOutputFile_("out_features", "<file>", "", "output file", false);
    setValidFormats_("out_features", ListUtils::create<String>("featureXML"));

    registerOutputFile_("out_tsv", "<file>", "", "TSV output file (mProphet-compatible TSV file)", false);
    setValidFormats_("out_tsv", ListUtils::create<String>("tsv"));

    registerOutputFile_("out_osw", "<file>", "", "OSW output file (PyProphet-compatible SQLite file)", false);
    setValidFormats_("out_osw", ListUtils::create<String>("osw"));

    registerOutputFile_("out_chrom", "<file>", "", "Also output all computed chromatograms output in mzML (chrom.mzML) or sqMass (SQLite format)", false, true);
    setValidFormats_("out_chrom", ListUtils::create<String>("mzML,sqMass"));

    // additional QC data
    registerOutputFile_("out_qc", "<file>", "", "Optional QC meta data (charge distribution in MS1). Only works with mzML input files.", false, true);
    setValidFormats_("out_qc", ListUtils::create<String>("json"));


    // misc options
    registerDoubleOption_("min_upper_edge_dist", "<double>", 0.0, "Minimal distance to the upper edge of a Swath window to still consider a precursor, in Thomson", false, true);
    registerFlag_("pasef", "data is PASEF data");

    // RT, mz and IM windows
    registerDoubleOption_("rt_extraction_window", "<double>", 600.0, "Only extract RT around this value (-1 means extract over the whole range, a value of 600 means to extract around +/- 300 s of the expected elution).", false);
    registerDoubleOption_("extra_rt_extraction_window", "<double>", 0.0, "Output an XIC with a RT-window by this much larger (e.g. to visually inspect a larger area of the chromatogram)", false, true);
    setMinFloat_("extra_rt_extraction_window", 0.0);
    registerDoubleOption_("ion_mobility_window", "<double>", -1, "Extraction window in ion mobility dimension (in 1/k0 or milliseconds depending on library). This is the full window size, e.g. a value of 10 milliseconds would extract 5 milliseconds on either side. -1 means extract over the whole range or ion mobility is not present. (Default for diaPASEF data: 0.06 1/k0)", false);
    registerDoubleOption_("mz_extraction_window", "<double>", 50, "Extraction window in Thomson or ppm (see mz_extraction_window_unit)", false);
    setMinFloat_("mz_extraction_window", 0.0);
    registerStringOption_("mz_extraction_window_unit", "<name>", "ppm", "Unit for mz extraction", false, true);
    setValidStrings_("mz_extraction_window_unit", ListUtils::create<String>("Th,ppm"));

    // MS1 mz windows and ion mobility
    registerDoubleOption_("mz_extraction_window_ms1", "<double>", 50, "Extraction window used in MS1 in Thomson or ppm (see mz_extraction_window_ms1_unit)", false);
    setMinFloat_("mz_extraction_window_ms1", 0.0);
    registerStringOption_("mz_extraction_window_ms1_unit", "<name>", "ppm", "Unit of the MS1 m/z extraction window", false, true);
    setValidStrings_("mz_extraction_window_ms1_unit", ListUtils::create<String>("ppm,Th"));
    registerDoubleOption_("im_extraction_window_ms1", "<double>", -1, "Extraction window in ion mobility dimension for MS1 (in 1/k0 or milliseconds depending on library). -1 means this is not ion mobility data.", false);

    registerStringOption_("use_ms1_ion_mobility", "<name>", "true", "Also perform precursor extraction using the same ion mobility window as for fragment ion extraction", false, true);
    setValidStrings_("use_ms1_ion_mobility", ListUtils::create<String>("true,false"));

    registerStringOption_("matching_window_only", "<name>", "false", "Assume the input data is targeted / PRM-like data with potentially overlapping DIA windows. Will only attempt to extract each assay from the *best* matching DIA window (instead of all matching windows).", false, true);
    setValidStrings_("matching_window_only", ListUtils::create<String>("true,false"));

    // iRT mz and IM windows
    registerDoubleOption_("irt_mz_extraction_window", "<double>", 50, "Extraction window used for iRT and m/z correction in Thomson or ppm (see irt_mz_extraction_window_unit)", false, true);
    setMinFloat_("irt_mz_extraction_window", 0.0);
    registerStringOption_("irt_mz_extraction_window_unit", "<name>", "ppm", "Unit for mz extraction", false, true);
    setValidStrings_("irt_mz_extraction_window_unit", ListUtils::create<String>("Th,ppm"));
    registerDoubleOption_("irt_im_extraction_window", "<double>", -1, "Ion mobility extraction window used for iRT (in 1/K0 or milliseconds depending on library). -1 means do not perform ion mobility calibration", false, true);

    registerDoubleOption_("min_rsq", "<double>", 0.95, "Minimum r-squared of RT peptides regression", false, true);
    registerDoubleOption_("min_coverage", "<double>", 0.6, "Minimum relative amount of RT peptides to keep", false, true);

    registerFlag_("split_file_input", "The input files each contain one single SWATH (alternatively: all SWATH are in separate files)", true);
    registerFlag_("use_elution_model_score", "Turn on elution model score (EMG fit to peak)", true);

    registerStringOption_("readOptions", "<name>", "normal", "Whether to run OpenSWATH directly on the input data, cache data to disk first or to perform a datareduction step first. If you choose cache, make sure to also set tempDirectory", false, true);
    setValidStrings_("readOptions", ListUtils::create<String>("normal,cache,cacheWorkingInMemory,workingInMemory"));

    registerStringOption_("mz_correction_function", "<name>", "none", "Use the retention time normalization peptide MS2 masses to perform a mass correction (linear, weighted by intensity linear or quadratic) of all spectra.", false, true);
    setValidStrings_("mz_correction_function", ListUtils::create<String>("none,regression_delta_ppm,unweighted_regression,weighted_regression,quadratic_regression,weighted_quadratic_regression,weighted_quadratic_regression_delta_ppm,quadratic_regression_delta_ppm"));

    registerStringOption_("tempDirectory", "<tmp>", File::getTempDirectory(), "Temporary directory to store cached files for example", false, true);

    registerStringOption_("extraction_function", "<name>", "tophat", "Function used to extract the signal", false, true);
    setValidStrings_("extraction_function", ListUtils::create<String>("tophat,bartlett"));

    registerIntOption_("batchSize", "<number>", 1000, "The batch size of chromatograms to process (0 means to only have one batch, sensible values are around 250-1000)", false, true);
    setMinInt_("batchSize", 0);
    registerIntOption_("outer_loop_threads", "<number>", -1, "How many threads should be used for the outer loop (-1 use all threads, use 4 to analyze 4 SWATH windows in memory at once).", false, true);

    registerIntOption_("ms1_isotopes", "<number>", 3, "The number of MS1 isotopes used for extraction", false, true);
    setMinInt_("ms1_isotopes", 0);

    registerSubsection_("Scoring", "Scoring parameters section");
    registerSubsection_("Library", "Library parameters section");

    registerSubsection_("RTNormalization", "Parameters for the RTNormalization for iRT petides. This specifies how the RT alignment is performed and how outlier detection is applied. Outlier detection can be done iteratively (by default) which removes one outlier per iteration or using the RANSAC algorithm.");
    registerSubsection_("Calibration", "Parameters for the m/z and ion mobility calibration.");
    registerTOPPSubsection_("Debugging", "Debugging");
    registerOutputFile_("Debugging:irt_mzml", "<file>", "", "Chromatogram mzML containing the iRT peptides", false);
    setValidFormats_("Debugging:irt_mzml", ListUtils::create<String>("mzML"));
    registerOutputFile_("Debugging:irt_trafo", "<file>", "", "Transformation file for RT transform", false);
    setValidFormats_("Debugging:irt_trafo", ListUtils::create<String>("trafoXML"));
  }

  Param getSubsectionDefaults_(const String& name) const override
  {
    if (name == "Scoring")
    {
      // set sensible default parameters
      Param feature_finder_param = MRMFeatureFinderScoring().getDefaults();
      feature_finder_param.remove("rt_extraction_window");
      feature_finder_param.setValue("stop_report_after_feature", 5);
      feature_finder_param.setValue("rt_normalization_factor", 100.0); // for iRT peptides between 0 and 100 (more or less)
      feature_finder_param.setValue("Scores:use_ms1_mi", "true");
      feature_finder_param.setValue("Scores:use_mi_score", "true");

      feature_finder_param.setValue("TransitionGroupPicker:min_peak_width", -1.0);
      feature_finder_param.setValue("TransitionGroupPicker:recalculate_peaks", "true");
      feature_finder_param.setValue("TransitionGroupPicker:compute_peak_quality", "false");
      feature_finder_param.setValue("TransitionGroupPicker:minimal_quality", -1.5);
      feature_finder_param.setValue("TransitionGroupPicker:background_subtraction", "none");
      feature_finder_param.setValue("TransitionGroupPicker:compute_peak_shape_metrics", "false");
      feature_finder_param.remove("TransitionGroupPicker:stop_after_intensity_ratio");

      // Peak Picker
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:use_gauss", "false");
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:sgolay_polynomial_order", 3);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:sgolay_frame_length", 11);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:peak_width", -1.0);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:remove_overlapping_peaks", "true");
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:write_sn_log_messages", "false"); // no log messages
      // TODO it seems that the legacy method produces slightly larger peaks, e.g. it will not cut off peaks too early
      // however the same can be achieved by using a relatively low SN cutoff in the -Scoring:TransitionGroupPicker:PeakPickerChromatogram:signal_to_noise 0.5
      feature_finder_param.setValue("TransitionGroupPicker:recalculate_peaks_max_z", 0.75);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:method", "corrected");
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:signal_to_noise", 0.1);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:gauss_width", 30.0);
      feature_finder_param.setValue("uis_threshold_sn", -1);
      feature_finder_param.setValue("uis_threshold_peak_area", 0);
      feature_finder_param.remove("TransitionGroupPicker:PeakPickerChromatogram:sn_win_len");
      feature_finder_param.remove("TransitionGroupPicker:PeakPickerChromatogram:sn_bin_count");
      feature_finder_param.remove("TransitionGroupPicker:PeakPickerChromatogram:stop_after_feature");

      // EMG Scoring - turn off by default since it is very CPU-intensive
      feature_finder_param.remove("Scores:use_elution_model_score");
      feature_finder_param.setValue("EMGScoring:max_iteration", 10);
      feature_finder_param.remove("EMGScoring:interpolation_step");
      feature_finder_param.remove("EMGScoring:tolerance_stdev_bounding_box");
      feature_finder_param.remove("EMGScoring:deltaAbsError");

      // remove these parameters
      feature_finder_param.remove("EMGScoring:statistics:mean");
      feature_finder_param.remove("EMGScoring:statistics:variance");
      return feature_finder_param;
    }
    else if (name == "RTNormalization")
    {
      Param p;

      p.setValue("alignmentMethod", "linear", "How to perform the alignment to the normalized RT space using anchor points. 'linear': perform linear regression (for few anchor points). 'interpolated': Interpolate between anchor points (for few, noise-free anchor points). 'lowess' Use local regression (for many, noisy anchor points). 'b_spline' use b splines for smoothing.");
      p.setValidStrings("alignmentMethod", {"linear","interpolated","lowess","b_spline"});
      p.setValue("lowess:span", 0.05, "Span parameter for lowess");
      p.setMinFloat("lowess:span", 0.0);
      p.setMaxFloat("lowess:span", 1.0);
      p.setValue("b_spline:num_nodes", 5, "Number of nodes for b spline");
      p.setMinInt("b_spline:num_nodes", 0);

      p.setValue("outlierMethod", "iter_residual", "Which outlier detection method to use (valid: 'iter_residual', 'iter_jackknife', 'ransac', 'none'). Iterative methods remove one outlier at a time. Jackknife approach optimizes for maximum r-squared improvement while 'iter_residual' removes the datapoint with the largest residual error (removal by residual is computationally cheaper, use this with lots of peptides).");
      p.setValidStrings("outlierMethod", {"iter_residual","iter_jackknife","ransac","none"});

      p.setValue("useIterativeChauvenet", "false", "Whether to use Chauvenet's criterion when using iterative methods. This should be used if the algorithm removes too many datapoints but it may lead to true outliers being retained.");
      p.setValidStrings("useIterativeChauvenet", {"true","false"});

      p.setValue("RANSACMaxIterations", 1000, "Maximum iterations for the RANSAC outlier detection algorithm.");
      p.setValue("RANSACMaxPercentRTThreshold", 3, "Maximum threshold in RT dimension for the RANSAC outlier detection algorithm (in percent of the total gradient). Default is set to 3% which is around +/- 4 minutes on a 120 gradient.");
      p.setValue("RANSACSamplingSize", 10, "Sampling size of data points per iteration for the RANSAC outlier detection algorithm.");

      p.setValue("estimateBestPeptides", "false", "Whether the algorithms should try to choose the best peptides based on their peak shape for normalization. Use this option you do not expect all your peptides to be detected in a sample and too many 'bad' peptides enter the outlier removal step (e.g. due to them being endogenous peptides or using a less curated list of peptides).");
      p.setValidStrings("estimateBestPeptides", {"true","false"});

      p.setValue("InitialQualityCutoff", 0.5, "The initial overall quality cutoff for a peak to be scored (range ca. -2 to 2)");
      p.setValue("OverallQualityCutoff", 5.5, "The overall quality cutoff for a peak to go into the retention time estimation (range ca. 0 to 10)");
      p.setValue("NrRTBins", 10, "Number of RT bins to use to compute coverage. This option should be used to ensure that there is a complete coverage of the RT space (this should detect cases where only a part of the RT gradient is actually covered by normalization peptides)");
      p.setValue("MinPeptidesPerBin", 1, "Minimal number of peptides that are required for a bin to counted as 'covered'");
      p.setValue("MinBinsFilled", 8, "Minimal number of bins required to be covered");
      return p;
    }
    else if (name == "Library")
    {
      return TransitionTSVFile().getDefaults();
    }
    else if (name == "Calibration")
    {
      Param p = SwathMapMassCorrection().getDefaults();
      p.remove("mz_extraction_window");
      p.remove("mz_extraction_window_ppm");
      p.remove("im_extraction_window");
      p.remove("mz_correction_function");
      return p;
    }
    else
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown subsection", name);
    }
  }

  ExitCodes main_(int, const char **) override
  {
    ///////////////////////////////////
    // Prepare Parameters
    ///////////////////////////////////
    StringList file_list = getStringList_("in");
    String tr_file = getStringOption_("tr");

    //tr_file input file type
    FileTypes::Type tr_type = FileTypes::nameToType(getStringOption_("tr_type"));
    if (tr_type == FileTypes::UNKNOWN)
    {
      tr_type = FileHandler::getType(tr_file);
      writeDebug_(String("Input file type (-tr): ") + FileTypes::typeToName(tr_type), 2);
    }

    if (tr_type == FileTypes::UNKNOWN)
    {
      writeLogError_("Error: Could not determine input file type for '-tr' !");
      return PARSE_ERROR;
    }

    String out = getStringOption_("out_features");
    String out_tsv = getStringOption_("out_tsv");
    String out_osw = getStringOption_("out_osw");

    String out_qc = getStringOption_("out_qc");


    String irt_tr_file = getStringOption_("tr_irt");
    String nonlinear_irt_tr_file = getStringOption_("tr_irt_nonlinear");
    String trafo_in = getStringOption_("rt_norm");
    String swath_windows_file = getStringOption_("swath_windows_file");

    String out_chrom = getStringOption_("out_chrom");
    bool split_file = getFlag_("split_file_input");
    bool use_emg_score = getFlag_("use_elution_model_score");
    bool force = getFlag_("force");
    bool pasef = getFlag_("pasef");
    bool sort_swath_maps = getFlag_("sort_swath_maps");
    bool use_ms1_traces = getStringOption_("enable_ms1") == "true";
    bool enable_uis_scoring = getStringOption_("enable_ipf") == "true";
    int batchSize = (int)getIntOption_("batchSize");
    int outer_loop_threads = (int)getIntOption_("outer_loop_threads");
    int ms1_isotopes = (int)getIntOption_("ms1_isotopes");
    Size debug_level = (Size)getIntOption_("debug");

    double min_rsq = getDoubleOption_("min_rsq");
    double min_coverage = getDoubleOption_("min_coverage");

    Param debug_params = getParam_().copy("Debugging:", true);

    String readoptions = getStringOption_("readOptions");
    String mz_correction_function = getStringOption_("mz_correction_function");

    // make sure tmp is a directory with proper separator at the end (downstream methods simply do path + filename)
    // (do not use QDir::separator(), since its platform specific (/ or \) while absolutePath() will always use '/')
    String tmp_dir = String(QDir(getStringOption_("tempDirectory").c_str()).absolutePath()).ensureLastChar('/');

    ///////////////////////////////////
    // Parameter validation
    ///////////////////////////////////

    bool load_into_memory = false;
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

    bool is_sqmass_input  = (FileHandler::getTypeByFileName(file_list[0]) == FileTypes::SQMASS);
    if (is_sqmass_input && !load_into_memory)
    {
      std::cout << "When using sqMass input files, it is highly recommended to use the workingInMemory option as otherwise data access will be very slow." << std::endl;
    }

    if (trafo_in.empty() && irt_tr_file.empty())
    {
      std::cout << "Since neither rt_norm nor tr_irt is set, OpenSWATH will " <<
        "not use RT-transformation (rather a null transformation will be applied)" << std::endl;
    }
    if ( int(!out.empty()) + int(!out_tsv.empty()) + int(!out_osw.empty()) != 1 )
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
      OPENMS_LOG_INFO << "Validate provided Swath windows file:" << std::endl;
      std::vector<double> swath_prec_lower;
      std::vector<double> swath_prec_upper;
      SwathWindowLoader::readSwathWindows(swath_windows_file, swath_prec_lower, swath_prec_upper);

      for (Size i = 0; i < swath_prec_lower.size(); i++)
      {
        OPENMS_LOG_DEBUG << "Read lower swath window " << swath_prec_lower[i] << " and upper window " << swath_prec_upper[i] << std::endl;
      }
    }

    double min_upper_edge_dist = getDoubleOption_("min_upper_edge_dist");
    bool use_ms1_im = getStringOption_("use_ms1_ion_mobility") == "true";
    bool prm = getStringOption_("matching_window_only") == "true";

    ChromExtractParams cp;
    cp.min_upper_edge_dist   = min_upper_edge_dist;
    cp.mz_extraction_window  = getDoubleOption_("mz_extraction_window");
    cp.ppm                   = getStringOption_("mz_extraction_window_unit") == "ppm";
    cp.rt_extraction_window  = getDoubleOption_("rt_extraction_window");
    cp.im_extraction_window  = getDoubleOption_("ion_mobility_window");
    cp.extraction_function   = getStringOption_("extraction_function");
    cp.extra_rt_extract      = getDoubleOption_("extra_rt_extraction_window");

    ChromExtractParams cp_irt = cp;
    cp_irt.rt_extraction_window = -1; // extract the whole RT range for iRT measurements
    cp_irt.mz_extraction_window = getDoubleOption_("irt_mz_extraction_window");
    cp_irt.im_extraction_window = getDoubleOption_("irt_im_extraction_window");
    cp_irt.ppm                  = getStringOption_("irt_mz_extraction_window_unit") == "ppm";

    ChromExtractParams cp_ms1 = cp;
    cp_ms1.mz_extraction_window  = getDoubleOption_("mz_extraction_window_ms1");
    cp_ms1.ppm                   = getStringOption_("mz_extraction_window_ms1_unit") == "ppm";
    cp_ms1.im_extraction_window  = (use_ms1_im) ? getDoubleOption_("im_extraction_window_ms1") : -1;

    Param feature_finder_param = getParam_().copy("Scoring:", true);
    feature_finder_param.setValue("use_ms1_ion_mobility", getStringOption_("use_ms1_ion_mobility"));

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

    bool compute_peak_shape_metrics = feature_finder_param.getValue("TransitionGroupPicker:compute_peak_shape_metrics").toBool();
    if (compute_peak_shape_metrics)
    {
      feature_finder_param.setValue("Scores:use_peak_shape_metrics", "true");
    }

    ///////////////////////////////////
    // Load the transitions
    ///////////////////////////////////
    OpenSwath::LightTargetedExperiment transition_exp = loadTransitionList(tr_type, tr_file, tsv_reader_param);
    OPENMS_LOG_INFO << "Loaded " << transition_exp.getProteins().size() << " proteins, " <<
      transition_exp.getCompounds().size() << " compounds with " << transition_exp.getTransitions().size() << " transitions." << std::endl;

    if (tr_type == FileTypes::PQP)
    {
      if (!out_osw.empty())
      { // copy the PQP file and name it OSW file
        std::ifstream  src(tr_file.c_str(), std::ios::binary);
        std::ofstream  dst(out_osw.c_str(), std::ios::binary | std::ios::trunc);
        dst << src.rdbuf();
      }
    }

    // If pasef flag is set, validate that IM is present
    if (pasef)
    {
      auto transitions = transition_exp.getTransitions();

      for ( Size k=0; k < (Size)transitions.size(); k++ )
      {
        if (transitions[k].precursor_im == -1)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Error: Transition " + transitions[k].getNativeID() +  " does not have a valid IM value, this must be set to use the -pasef flag");
        }
      }
    }

    ///////////////////////////////////
    // Load the SWATH files
    ///////////////////////////////////
    boost::shared_ptr<ExperimentalSettings> exp_meta(new ExperimentalSettings);
    std::vector< OpenSwath::SwathMap > swath_maps;

    // collect some QC data
    if (!out_qc.empty())
    {
      OpenSwath::SwathQC qc(30, 0.04);
      MSDataTransformingConsumer qc_consumer; // apply some transformation
      qc_consumer.setSpectraProcessingFunc(qc.getSpectraProcessingFunc());
      qc_consumer.setExperimentalSettingsFunc(qc.getExpSettingsFunc());
      if (!loadSwathFiles(file_list, exp_meta, swath_maps, split_file, tmp_dir, readoptions,
                          swath_windows_file, min_upper_edge_dist, force,
                          sort_swath_maps, prm, pasef, &qc_consumer))
      {
        return PARSE_ERROR;
      }
      qc.storeJSON(out_qc);
    }
    else
    {
      if (!loadSwathFiles(file_list, exp_meta, swath_maps, split_file, tmp_dir, readoptions,
                          swath_windows_file, min_upper_edge_dist, force,
                          sort_swath_maps, prm, pasef))
      {
        return PARSE_ERROR;
      }
    }


    ///////////////////////////////////
    // Get the transformation information (using iRT peptides)
    ///////////////////////////////////
    String irt_trafo_out = debug_params.getValue("irt_trafo").toString();
    String irt_mzml_out = debug_params.getValue("irt_mzml").toString();
    Param irt_detection_param = getParam_().copy("RTNormalization:", true);
    Param calibration_param = getParam_().copy("Calibration:", true);
    calibration_param.setValue("mz_extraction_window", cp_irt.mz_extraction_window);
    calibration_param.setValue("mz_extraction_window_ppm", cp_irt.ppm ? "true" : "false");
    calibration_param.setValue("im_extraction_window", cp_irt.im_extraction_window);
    calibration_param.setValue("mz_correction_function", mz_correction_function);
    TransformationDescription trafo_rtnorm;
    if (nonlinear_irt_tr_file.empty())
    {
      trafo_rtnorm = performCalibration(trafo_in, irt_tr_file, swath_maps,
                                        min_rsq, min_coverage, feature_finder_param,
                                        cp_irt, irt_detection_param, calibration_param,
                                        debug_level, pasef, load_into_memory,
                                        irt_trafo_out, irt_mzml_out);
    }
    else
    {
      ///////////////////////////////////
      // First perform a simple linear transform, then do a second, nonlinear one
      ///////////////////////////////////

      Param linear_irt = irt_detection_param;
      linear_irt.setValue("alignmentMethod", "linear");
      Param no_calibration = calibration_param;
      no_calibration.setValue("mz_correction_function", "none");
      trafo_rtnorm = performCalibration(trafo_in, irt_tr_file, swath_maps,
                                        min_rsq, min_coverage, feature_finder_param,
                                        cp_irt, linear_irt, no_calibration,
                                        debug_level, pasef, load_into_memory,
                                        irt_trafo_out, irt_mzml_out);

      cp_irt.rt_extraction_window = 900; // extract some substantial part of the RT range (should be covered by linear correction)
      cp_irt.rt_extraction_window = 600; // extract some substantial part of the RT range (should be covered by linear correction)

      ///////////////////////////////////
      // Get the secondary transformation (nonlinear)
      ///////////////////////////////////
      OpenSwath::LightTargetedExperiment transition_exp_nl;
      transition_exp_nl = loadTransitionList(FileHandler::getType(nonlinear_irt_tr_file), nonlinear_irt_tr_file, tsv_reader_param);

      std::vector< OpenMS::MSChromatogram > chromatograms;
      OpenSwathCalibrationWorkflow wf;
      wf.setLogType(log_type_);
      wf.simpleExtractChromatograms_(swath_maps, transition_exp_nl, chromatograms,
                                    trafo_rtnorm, cp_irt, pasef, load_into_memory);

      // always use estimateBestPeptides for the nonlinear approach
      Param nonlinear_irt = irt_detection_param;
      nonlinear_irt.setValue("estimateBestPeptides", "true");

      TransformationDescription im_trafo; // exp -> theoretical
      trafo_rtnorm = wf.doDataNormalization_(transition_exp_nl, chromatograms, im_trafo, swath_maps,
                                             min_rsq, min_coverage,
                                             feature_finder_param, nonlinear_irt, calibration_param, pasef);

      TransformationDescription im_trafo_inv = im_trafo;
      im_trafo_inv.invert(); // theoretical -> experimental

      // We now modify the library as this is the easiest thing to do
      for (auto & p : transition_exp.getCompounds())
      {
        p.drift_time = im_trafo_inv.apply(p.drift_time);
      }

    }

    ///////////////////////////////////
    // Set up chromatogram output
    // Either use chrom.mzML or sqliteDB (sqMass)
    ///////////////////////////////////
    Interfaces::IMSDataConsumer* chromatogramConsumer;
    UInt64 run_id = OpenMS::UniqueIdGenerator::getUniqueId();
    prepareChromOutput(&chromatogramConsumer, exp_meta, transition_exp, out_chrom, run_id);

    ///////////////////////////////////
    // Set up peakgroup file output (.tsv or .osw file)
    ///////////////////////////////////
    FeatureMap out_featureFile;
    OpenSwathTSVWriter tsvwriter(out_tsv, file_list[0], use_ms1_traces); // only active if filename not empty
    OpenSwathOSWWriter oswwriter(out_osw, run_id, file_list[0], enable_uis_scoring); // only active if filename not empty

    OpenSwathWorkflow wf(use_ms1_traces, use_ms1_im, prm, pasef, outer_loop_threads);
    wf.setLogType(log_type_);
    wf.performExtraction(swath_maps, trafo_rtnorm, cp, cp_ms1, feature_finder_param, transition_exp,
        out_featureFile, !out.empty(), tsvwriter, oswwriter, chromatogramConsumer, batchSize, ms1_isotopes, load_into_memory);

    if (!out.empty())
    {
      addDataProcessing_(out_featureFile, getProcessingInfo_(DataProcessing::QUANTITATION));
      out_featureFile.ensureUniqueId();
      FileHandler().storeFeatures(out, out_featureFile, {FileTypes::FEATUREXML});
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
