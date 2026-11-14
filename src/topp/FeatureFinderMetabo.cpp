// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Erhan Kenar, Holger Franken $
// --------------------------------------------------------------------------
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/IONMOBILITY/FAIMSHelper.h>
#include <OpenMS/IONMOBILITY/IMDataConverter.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/FEATUREFINDER/MassTraceDetection.h>
#include <OpenMS/FEATUREFINDER/ElutionPeakDetection.h>
#include <OpenMS/FEATUREFINDER/FeatureFindingMetabo.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/PROCESSING/FEATURE/FeatureOverlapFilter.h>

#include <limits>
#include <set>
#include <unordered_set>


using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_FeatureFinderMetabo FeatureFinderMetabo

@brief FeatureFinderMetabo assembles metabolite features from singleton mass traces.

<CENTER>
<table>
<tr>
<th ALIGN = "center"> pot. predecessor tools </td>
<td VALIGN="middle" ROWSPAN=3> &rarr; FeatureFinderMetabo &rarr;</td>
<th ALIGN = "center"> pot. successor tools </td>
</tr>
<tr>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerHiRes </td>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_TextExporter</td>
</tr>
<tr>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerHiRes </td>
</tr>
</table>
</CENTER>

Mass traces alone would allow for further analysis such as metabolite ID or
statistical evaluation. However, in general, monoisotopic mass traces are
accompanied by satellite C13 peaks and thus may render the analysis more
difficult. FeatureFinderMetabo fulfills a further data reduction step by
assembling compatible mass traces to metabolite features (that is, all mass
traces originating from one metabolite). To this end, multiple metabolite
hypotheses are formulated and scored according to how well differences in RT (optional),
m/z or intensity ratios match to those of theoretical isotope patterns.

If the raw data scans contain the scan polarity information, it is stored as
meta value "scan_polarity" in the output file.

Mass trace clustering can be done using either 13C distances or a linear model (Kenar et al) -- see parameter 'ffm:mz_scoring_13C'.
Generally, for lipidomics, use 13C, since lipids contain a lot of 13C.
For general metabolites, the linear model is usually more appropriate.
To decide what is better, the total number of features can be used as indirect measure
- the lower(!) the better (since more mass traces are assembled into single features).
Detailed information is stored in the featureXML output: it contains meta-values for each feature about the 
mass trace differences (inspectable via TOPPView). If you want this in a tabular format, use TextExporter, i.e.,
@code
   TextExporter.exe -feature:add_metavalues 1 -in <ff_metabo.featureXML> -out <ff_metabo.csv>
@endcode
By default, the linear model is used.

@section ffmetabo_bounded_traces Bounded mass traces and parallel extraction

The adaptive m/z window is a statistical estimate, not a physical mass-error
limit: @p algorithm:mtd:mass_error_ppm initializes the standard deviation, and
extension searches within three current standard deviations. Interference can
inflate that estimate, admitting still more interference and moving the centroid.
The new @p algorithm:mtd:max_apex_deviation_ppm parameter therefore additionally
bounds every collected peak relative to the immutable seed/apex m/z when enabled.
The cap is disabled by default (zero), preserving existing extraction behavior
in all callers. Candidate searches intersect both windows before choosing the
nearest peak, so an excluded peak cannot hide an admissible alternative.
Enabling the cap can change trace and feature counts, particularly on weak
chromatographic flanks. Validate an instrument-appropriate bound using tail
recovery, false joins and quantification bias before changing a workflow's
defaults; the 20 ppm benchmark setting below is not a validated recommendation.

For example, unbounded extraction of the public, centroided repository excerpt
@c doc/pyopenms/src/data/Metabolomics_2.mzML produced trace T1930 with 119 points
over 211.0623 seconds, spanning m/z 694.2455444 to 695.4852295 around centroid
694.8168204: a 1784.19 ppm span. Dynamic variance re-estimation had allowed
unrelated signals to join; a tolerance initialized to 20 ppm did not prevent it.
This is an example of pathological assembly, not evidence that these observations
belong to one analyte. The fixed-apex cap prevents such unbounded drift without
changing the variance estimator.

With a positive cap and @c -threads @c 2 (or more), outlier-mode extraction first builds
conservative independent components. All above-noise peaks participate, including
weak peaks that cannot seed a trace. If a trace is confined to seed m/z +/- H,
two admissible observations differ by at most 2H in the seed's mass coordinate.
A symmetric, outward-rounded twice-cap gate conservatively connects such pairs.
Edges cover the same scan and up to @p trace_termination_outliers + 1 nonempty
working scans: empty scans do not increment the existing miss counter. There
is no additional RT cutoff or change to sample-rate accounting.

Each component retains the original greedy seed priority, including intensity
ties, and writes distinct byte-addressed peak-ownership entries. No peak locks
or overlapping-region reconciliation are needed. Results are merged in original
seed order before assigning labels, giving exactly the capped serial result.
The largest component remains a limit on parallel speedup. Additional component
storage is linear in the number of above-noise peaks; spectra are shared rather
than copied per worker.

@c -threads @c 1 remains the default. Extraction stays serial with a disabled
cap, sample-rate termination, a @c max_traces limit in the library API, or inside
an existing OpenMP parallel region. The later elution-peak detection and isotope
assembly stages retain their own threading behavior.

The standalone @c tools/benchmark_mtd utility compares complete ordered traces,
peak coordinates/intensities, centroids, dispersion, areas and available metadata
on every repetition. Timings include filtering, seed sorting, component discovery,
extraction and merging, but exclude file I/O and equality checks. It reports both
first-call and interleaved repeated-run timings. Against an installed OpenMS:
@code
cmake -S tools/benchmark_mtd -B benchmark-mtd-build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=<installation>
cmake --build benchmark-mtd-build -j2
benchmark-mtd-build/benchmark_mtd doc/pyopenms/src/data/Metabolomics_1.mzML 2 9 20
benchmark-mtd-build/benchmark_mtd doc/pyopenms/src/data/Metabolomics_2.mzML 2 9 20
# Reproduce the unbounded trace example:
benchmark-mtd-build/benchmark_mtd doc/pyopenms/src/data/Metabolomics_2.mzML 2 9 0
@endcode

A native Release (-O2, GCC 13.3) run on a shared Linux environment with two
available CPUs gave the following medians over nine interleaved repetitions,
using the 20 ppm cap. These are mass-trace extraction timings, not full-pipeline
speedups:

<table>
<tr><th>Public excerpt</th><th>Traces</th><th>1 thread</th><th>2 threads</th><th>Speedup</th></tr>
<tr><td>Metabolomics_1.mzML</td><td>6806</td><td>0.3805 s</td><td>0.2066 s</td><td>1.84x</td></tr>
<tr><td>Metabolomics_2.mzML</td><td>5423</td><td>0.4289 s</td><td>0.2391 s</td><td>1.79x</td></tr>
</table>

All ordered trace values matched exactly on every repetition. First-call
serial/parallel times were 0.3992/0.2293 s and 0.4449/0.2383 s, respectively;
timings are illustrative, not a scaling guarantee on other machines.
Maximum observed trace spans were 39.90 and 39.91 ppm, rather than the
unbounded example's 1784.19 ppm. Whole benchmark-process peak RSS was about
107 and 114 MiB, including input, reference/comparison outputs and libraries;
this is not incremental partitioning memory or a per-worker allocation.
The complete FeatureFinderMetabo output on the second excerpt also matched
between one and two threads with zero numeric tolerance (ignoring feature IDs).

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_FeatureFinderMetabo.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_FeatureFinderMetabo.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFeatureFinderMetabo :
  public TOPPBase
{
public:
  TOPPFeatureFinderMetabo() :
    TOPPBase("FeatureFinderMetabo", "Assembles metabolite features from centroided (LC-)MS data using the mass trace approach.")
  {
  }

protected:

  /**
   * @brief Process a single PeakMap (or FAIMS CV group) through the full FFMetabo pipeline
   * @param[in,out] ms_peakmap Input peak map (will be modified - sorted)
   * @param[in] common_param Common parameters
   * @param[in] mtd_param MassTraceDetection parameters
   * @param[in] epd_param ElutionPeakDetection parameters
   * @param[in] ffm_param FeatureFindingMetabo parameters
   * @param[out] feat_map Output feature map
   * @param[out] feat_chromatograms Output chromatograms (if enabled)
   * @return True on success, false on error
   */
  bool processOneGroup_(PeakMap& ms_peakmap,
                        const Param& common_param,
                        Param mtd_param,
                        Param epd_param,
                        Param ffm_param,
                        FeatureMap& feat_map,
                        std::vector<std::vector<OpenMS::MSChromatogram>>& feat_chromatograms)
  {
    // make sure the spectra are sorted by m/z
    ms_peakmap.sortSpectra(true);

    vector<MassTrace> m_traces;

    //-------------------------------------------------------------
    // configure and run mass trace detection
    //-------------------------------------------------------------
    MassTraceDetection mtdet;
    mtd_param.insert("", common_param);
    mtd_param.remove("chrom_fwhm");
    mtdet.setParameters(mtd_param);
    mtdet.run(ms_peakmap, m_traces);

    //-------------------------------------------------------------
    // configure and run elution peak detection
    //-------------------------------------------------------------
    std::vector<MassTrace> m_traces_final;
    if (epd_param.getValue("enabled").toBool())
    {
      std::vector<MassTrace> splitted_mtraces;
      epd_param.remove("enabled");
      epd_param.insert("", common_param);
      epd_param.remove("noise_threshold_int");
      ElutionPeakDetection epdet;
      epdet.setParameters(epd_param);
      epdet.detectPeaks(m_traces, splitted_mtraces);
      if (epdet.getParameters().getValue("width_filtering") == "auto")
      {
        m_traces_final.clear();
        epdet.filterByPeakWidth(splitted_mtraces, m_traces_final);
      }
      else
      {
        m_traces_final = splitted_mtraces;
      }
    }
    else
    {
      m_traces_final = m_traces;
      for (Size i = 0; i < m_traces_final.size(); ++i)
      {
        m_traces_final[i].estimateFWHM(false);
      }
      if (ffm_param.getValue("use_smoothed_intensities").toBool())
      {
        OPENMS_LOG_WARN << "Without EPD, smoothing is not supported. Setting 'use_smoothed_intensities' to false!" << std::endl;
        ffm_param.setValue("use_smoothed_intensities", "false");
      }
    }

    //-------------------------------------------------------------
    // configure and run feature finding
    //-------------------------------------------------------------
    ffm_param.insert("", common_param);
    ffm_param.remove("noise_threshold_int");
    ffm_param.remove("chrom_peak_snr");

    FeatureFindingMetabo ffmet;
    ffmet.setParameters(ffm_param);
    ffmet.run(m_traces_final, feat_map, feat_chromatograms);

    Size trace_count(0);
    for (Size i = 0; i < feat_map.size(); ++i)
    {
      OPENMS_PRECONDITION(feat_map[i].metaValueExists(Constants::UserParam::NUM_OF_MASSTRACES),
          "MetaValue 'num_of_masstraces' missing from FFMetabo output!");
      trace_count += (Size) feat_map[i].getMetaValue(Constants::UserParam::NUM_OF_MASSTRACES);
    }

    if (trace_count != m_traces_final.size())
    {
      if (!ffm_param.getValue("remove_single_traces").toBool())
      {
        OPENMS_LOG_ERROR << "FF-Metabo: Internal error. Not all mass traces have been assembled to features!" << std::endl;
        return false;
      }
      else
      {
        OPENMS_LOG_INFO << "FF-Metabo: " << (m_traces_final.size() - trace_count) << " unassembled traces have been removed." << std::endl;
      }
    }

    OPENMS_LOG_INFO << "-- FF-Metabo stats --\n"
             << "Input traces:    " << m_traces_final.size() << "\n"
             << "Output features: " << feat_map.size() << " (total trace count: " << trace_count << ")" << std::endl;

    return true;
  }

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Centroided mzML file");
    setValidFormats_("in", {"mzML",
#ifdef WITH_THERMO_RAW
      "raw",
#endif
    });
    registerOutputFile_("out", "<file>", "", "FeatureXML file with metabolite features");
    setValidFormats_("out", ListUtils::create<std::string>("featureXML"));

    registerOutputFile_("out_chrom", "<file>", "", "Optional mzML file with chromatograms", false);
    setValidFormats_("out_chrom", ListUtils::create<std::string>("mzML"));

    addEmptyLine_();
    registerStringOption_("faims_merge_features", "<true/false>", "true",
      "For FAIMS data with multiple compensation voltages: Merge features representing the same analyte "
      "detected at different CV values into a single feature. Only features with DIFFERENT FAIMS CV values "
      "are merged (same CV = different analytes). Has no effect on non-FAIMS data.", false);
    setValidStrings_("faims_merge_features", {"true", "false"});

    addEmptyLine_();
    registerSubsection_("algorithm", "Algorithm parameters section");
  }

  Param getSubsectionDefaults_(const std::string& /*section*/) const override
  {
    Param combined;
    Param p_com;
    p_com.setValue("noise_threshold_int", 10.0, "Intensity threshold below which peaks are regarded as noise.");
    p_com.setValue("chrom_peak_snr", 3.0, "Minimum signal-to-noise a mass trace should have.");
    p_com.setValue("chrom_fwhm", 5.0, "Expected chromatographic peak width (in seconds).");
    combined.insert("common:", p_com);
    combined.setSectionDescription("common", "Common parameters for all other subsections");

    Param p_mtd = MassTraceDetection().getDefaults();
    p_mtd.remove("noise_threshold_int");
    p_mtd.remove("chrom_peak_snr");
    combined.insert("mtd:", p_mtd);
    combined.setSectionDescription("mtd", "Mass Trace Detection parameters");

    Param p_epd;
    p_epd.setValue("enabled", "true", "Enable splitting of isobaric mass traces by chromatographic peak detection. Disable for direct injection.");
    p_epd.setValidStrings("enabled", {"true","false"});
    p_epd.insert("", ElutionPeakDetection().getDefaults());
    p_epd.remove("chrom_peak_snr");
    p_epd.remove("chrom_fwhm");

    combined.insert("epd:", p_epd);
    combined.insert("epd:", p_epd);
    combined.setSectionDescription("epd", "Elution Profile Detection (to separate isobaric Mass Traces by elution time).");

    Param p_ffm = FeatureFindingMetabo().getDefaults();
    p_ffm.remove("chrom_fwhm");
    p_ffm.remove("report_chromatograms");
    combined.insert("ffm:", p_ffm);
    combined.setSectionDescription("ffm", "FeatureFinder parameters (assembling mass traces to charged features)");

    return combined;
  }

  ExitCodes main_(int, const char**) override
  {

    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    std::string in = getStringOption_("in");
    std::string out = getStringOption_("out");
    std::string out_chrom = getStringOption_("out_chrom");

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------
    FileHandler mz_data_file;
    PeakMap ms_peakmap;
    std::vector<Int> ms_level(1, 1);
    mz_data_file.getOptions().setMSLevels(ms_level);
    mz_data_file.loadExperiment(in, ms_peakmap, {FileTypes::MZML, FileTypes::RAW}, log_type_);

    if (ms_peakmap.empty())
    {
      OPENMS_LOG_WARN << "The given file does not contain any conventional peak data, but might"
                  " contain chromatograms. This tool currently cannot handle them, sorry.";
      return INCOMPATIBLE_INPUT_DATA;
    }

    // determine type of spectral data (profile or centroided)
    SpectrumSettings::SpectrumType spectrum_type = ms_peakmap[0].getType();

    if (spectrum_type == SpectrumSettings::SpectrumType::PROFILE)
    {
      if (!getFlag_("force"))
      {
        throw OpenMS::Exception::FileEmpty(__FILE__, __LINE__, __FUNCTION__,
            "Error: Profile data provided but centroided spectra expected. To enforce processing of the data set the -force flag.");
      }
    }

    //-------------------------------------------------------------
    // set parameters
    //-------------------------------------------------------------

    Param common_param = getParam_().copy("algorithm:common:", true);
    writeDebug_("Common parameters passed to sub-algorithms (mtd and ffm)", common_param, 3);

    Param mtd_param = getParam_().copy("algorithm:mtd:", true);
    writeDebug_("Parameters passed to MassTraceDetection", mtd_param, 3);

    Param epd_param = getParam_().copy("algorithm:epd:", true);
    writeDebug_("Parameters passed to ElutionPeakDetection", epd_param, 3);

    Param ffm_param = getParam_().copy("algorithm:ffm:", true);
    writeDebug_("Parameters passed to FeatureFindingMetabo", ffm_param, 3);

    std::string report_chromatograms = out_chrom.empty() ? "false" : "true";
    ffm_param.setValue("report_chromatograms", report_chromatograms);

    // Store polarities before potentially moving ms_peakmap (needed later for output annotation)
    set<IonSource::Polarity> polarities;
    for (const auto& spec : ms_peakmap)
    {
      polarities.insert(spec.getInstrumentSettings().getPolarity());
    }

    //-------------------------------------------------------------
    // Split by FAIMS CV (returns single NaN-keyed element for non-FAIMS data)
    //-------------------------------------------------------------
    auto faims_groups = IMDataConverter::splitByFAIMSCV(std::move(ms_peakmap));

    const bool has_faims = faims_groups.size() > 1 || !std::isnan(faims_groups[0].first);
    if (has_faims)
    {
      OPENMS_LOG_INFO << "FAIMS data detected with " << faims_groups.size() << " compensation voltage(s)." << endl;
    }

    FeatureMap feat_map;
    std::vector<std::vector<OpenMS::MSChromatogram>> feat_chromatograms;

    // Process each FAIMS CV group (or single group for non-FAIMS data)
    for (auto& [group_cv, faims_group] : faims_groups)
    {
      if (has_faims)
      {
        OPENMS_LOG_INFO << "Processing FAIMS CV group: " << group_cv << " V (" << faims_group.size() << " spectra)" << endl;
      }

      // Process this group
      FeatureMap feat_map_cv;
      std::vector<std::vector<OpenMS::MSChromatogram>> feat_chromatograms_cv;

      if (!processOneGroup_(faims_group, common_param, mtd_param, epd_param, ffm_param,
                            feat_map_cv, feat_chromatograms_cv))
      {
        if (has_faims)
        {
          OPENMS_LOG_ERROR << "Processing failed for FAIMS CV " << group_cv << " V. Aborting." << endl;
        }
        return UNEXPECTED_RESULT;
      }

      // Annotate features with FAIMS CV (if FAIMS data) and add to combined results
      for (auto& feat : feat_map_cv)
      {
        if (has_faims)
        {
          feat.setMetaValue(Constants::UserParam::FAIMS_CV, group_cv);
        }
        feat_map.push_back(feat);
      }

      // Combine chromatograms
      for (auto& chrom_group : feat_chromatograms_cv)
      {
        feat_chromatograms.push_back(std::move(chrom_group));
      }
    }

    if (has_faims)
    {
      OPENMS_LOG_INFO << "Combined " << feat_map.size() << " features from all FAIMS CV groups." << endl;

      // Optionally merge features representing the same analyte at different CV values
      if (getStringOption_("faims_merge_features") == "true")
      {
        Size before_merge = feat_map.size();
        FeatureOverlapFilter::mergeFAIMSFeatures(feat_map, 5.0, 0.05);
        OPENMS_LOG_INFO << "FAIMS feature merge: " << before_merge << " -> " << feat_map.size()
                        << " features (merged " << (before_merge - feat_map.size()) << ")" << endl;
      }
    }

    // filter features with zero intensity (this can happen if the FWHM is zero (bc of overly skewed shape) and no peaks end up being summed up)
    auto intensity_zero = [&](Feature& f) { return f.getIntensity() == 0; };
    feat_map.erase(remove_if(feat_map.begin(),feat_map.end(),intensity_zero),feat_map.end());

    // keep the chromatograms of the remaining features (the FAIMS merge and the filter removed some). Each
    // chromatogram carries the unique id of its feature in its native ID ("<id>_<trace>", see
    // FeatureHypothesis::getChromatograms); the order does not pair them, FeatureFindingMetabo sorts features by m/z
    if (!out_chrom.empty())
    {
      std::unordered_set<UInt64> kept_ids;
      for (const Feature& f : feat_map) kept_ids.insert(f.getUniqueId());
      std::erase_if(feat_chromatograms, [&kept_ids](const std::vector<MSChromatogram>& chroms)
      {
        if (chroms.empty()) return true;
        const std::string& id = chroms.front().getNativeID();
        return !kept_ids.contains(std::stoull(id.substr(0, id.find('_'))));
      });
    }

    // store chromatograms
    if (!out_chrom.empty())
    {
      if (feat_chromatograms.size() == feat_map.size())
        {
          MSExperiment out_exp;
            for (Size i = 0; i < feat_chromatograms.size(); ++i)
            {
                for (Size j = 0; j < feat_chromatograms[i].size(); ++j)
                {
                  out_exp.addChromatogram(feat_chromatograms[i][j]);
                }
            }
          FileHandler().storeExperiment(out_chrom, out_exp, {FileTypes::MZML});
        }
        else
        {
            OPENMS_LOG_ERROR << "FF-Metabo: Internal error. The number of features (" << feat_chromatograms.size() << ") and chromatograms (" << feat_map.size() << ") are different! Aborting." << std::endl;
            return UNEXPECTED_RESULT;
        }
    }

    // store ionization mode of spectra (useful for post-processing by AccurateMassSearch tool)
    // Note: polarities were collected before ms_peakmap was potentially moved
    if (!feat_map.empty() && !polarities.empty())
    {
      StringList sl_pols;
      for (const auto& pol : polarities)
      {
        sl_pols.push_back(IonSource::polarityToString(pol));
      }
      feat_map[0].setMetaValue("scan_polarity", ListUtils::concatenate(sl_pols, ";"));
    }

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    // ensure unique IDs for the combined feature map
    feat_map.ensureUniqueId();

    // annotate output with data processing info
    addDataProcessing_(feat_map, getProcessingInfo_(DataProcessing::QUANTITATION));

    // annotate "spectra_data" metavalue
    // Note: use simple form since ms_peakmap may have been moved for FAIMS processing
    if (getFlag_("test"))
    {
      // if test mode set, add file without path so we can compare it
      feat_map.setPrimaryMSRunPath({"file://" + File::basename(in)});
    }
    else
    {
      feat_map.setPrimaryMSRunPath({in});
    }    

    FileHandler().storeFeatures(out, feat_map, {FileTypes::FEATUREXML});
  
    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPFeatureFinderMetabo tool;
  return tool.main(argc, argv);
}

/// @endcond
