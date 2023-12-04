// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Jihyung Kim $
// $Authors: Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHQuantAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHQuantHelper.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/FILTERING/DATAREDUCTION/ElutionPeakDetection.h>
#include <OpenMS/FILTERING/DATAREDUCTION/FeatureFindingMetabo.h>
#include <OpenMS/FILTERING/DATAREDUCTION/MassTraceDetection.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MassTrace.h>


using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------
/**
    @page TOPP_FLASHQuant TOPP_FLASHQuant

    @brief TOPP_FLASHQuant The intact protein feature detection for quantification (centroided).
 */

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFLASHQuant :
    public TOPPBase,
    public ProgressLogger
{
public:
  typedef FLASHQuantHelper::FeatureGroup FeatureGroup;
  typedef FLASHQuantHelper::FeatureSeed FeatureSeed;

  TOPPFLASHQuant():
    TOPPBase("FLASHQuant", "The intact protein feature detection for quantification", false, {}), ProgressLogger()
  {
    this->setLogType(CMD);
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "MzML input file", true);
    setValidFormats_("in", ListUtils::create<String>("mzML"));

    registerOutputFile_("out", "<file>", "", "Tsv output file with quantified feature groups (putative proteoform)", true);
    setValidFormats_("out", ListUtils::create<String>("tsv"));

    registerOutputFile_("out_feat", "<file>", "", "FeatureXML output file with quantified feature groups (putative proteoform)", false);
    setValidFormats_("out_feat", ListUtils::create<String>("featureXML"));

    registerOutputFile_("out_detail", "<file>", "", "Tsv output file with mass trace information per feature group", false);
    setValidFormats_("out_detail", ListUtils::create<String>("tsv"));

    addEmptyLine_();
    registerSubsection_("algorithm", "Algorithm parameters section");
  }

  Param getSubsectionDefaults_(const String& /*section*/) const override
  {
    Param combined;

    Param p_mtd = MassTraceDetection().getDefaults();
    p_mtd.setValue("noise_threshold_int", 0.0);
    p_mtd.setValue("chrom_peak_snr", 0.0);
    p_mtd.setValue("mass_error_ppm", 5.0);
    combined.insert("mtd:", p_mtd);
    combined.setSectionDescription("mtd", "Mass Trace Detection parameters");

    Param p_epd = ElutionPeakDetection().getDefaults();
    p_epd.setValue("width_filtering", "auto");
    combined.insert("epd:", p_epd);
    combined.setSectionDescription("epd", "Elution Profile Detection (to separate isobaric Mass Traces by elution time).");

    Param p_ffi = FLASHQuantAlgorithm().getDefaults();
    combined.insert("fdq:", p_ffi);
    combined.setSectionDescription("fdq", "FLASHQuant parameters (assembling mass traces to charged features)");

    return combined;
  }

  void storeFeatureGroupInOpenMSFeature(std::vector<FeatureGroup> &feature_groups, FeatureMap &out_featmap) const
  {
    out_featmap.clear();
    for (auto &fgroup: feature_groups)
    {
      // create OpenMS::Feature per charge
      std::vector<Feature> feat_vec;
      for (auto &cs : fgroup.getChargeSet())
      {
        Feature feat;
        feat.setCharge(cs);
        feat.setOverallQuality(fgroup.getIsotopeCosineOfCharge(cs));
        feat.setIntensity(fgroup.getIntensityOfCharge(cs));
        feat.setMetaValue("monoisotopic_mass_of_feature", fgroup.getMonoisotopicMass());
        feat.setMetaValue("feature_group_score", fgroup.getFeatureGroupScore());

        std::vector<ConvexHull2D> tmp_hulls;
        std::vector<std::vector<double>> intensity_of_hulls;
        FeatureSeed* apex_ptr;
        double fwhm_start = LONG_MAX;
        double fwhm_end = .0;
        double max_intensity = .0;
        for (auto &seed: fgroup)
        {
          if (seed.getCharge() != cs)
          {
            continue;
          }

          // get apex information
          if (max_intensity < seed.getIntensity())
          {
            max_intensity = seed.getIntensity();
            apex_ptr = &seed;
          }

          // get fwhm information
          if (seed.getFwhmStart() < fwhm_start)
          {
            fwhm_start = seed.getFwhmStart();
          }
          if (seed.getFwhmEnd() > fwhm_end)
          {
            fwhm_end = seed.getFwhmEnd();
          }

          // generate ConvexHull2D from FeatureSeed
          const MassTrace& mt_ptr = seed.getMassTrace();
          ConvexHull2D::PointArrayType hull_points(mt_ptr.getSize());
          std::vector<double> intensities;

          Size i = 0;
          for (MassTrace::const_iterator l_it = mt_ptr.begin(); l_it != mt_ptr.end(); ++l_it)
          {
            hull_points[i][0] = (*l_it).getRT();
            hull_points[i][1] = (*l_it).getMZ();
            intensities.push_back((*l_it).getIntensity());
            ++i;
          }

          ConvexHull2D hull;
          hull.addPoints(hull_points);
          tmp_hulls.push_back(hull);
          intensity_of_hulls.push_back(intensities);
        }
        if (tmp_hulls.empty()) // if this feature is empty
        {
          continue;
        }

        // store calculated information
        feat.setConvexHulls(tmp_hulls);
        feat.setMZ(apex_ptr->getCentroidMz());
        feat.setRT(apex_ptr->getMassTrace().getCentroidRT());
        feat.setWidth(fwhm_end-fwhm_start);
        feat.setMetaValue("num_of_masstraces", intensity_of_hulls.size());

        int i = 1;
        for (auto& inty_vec: intensity_of_hulls)
        {
          String meta_label = "masstrace_intensity_" + std::to_string(i);
          feat.setMetaValue(meta_label, inty_vec);
          ++i;
        }
        feat.applyMemberFunction(&UniqueIdInterface::setUniqueId);

        // add features to output FeatureMap
        out_featmap.push_back(feat);
      }
    }
    out_featmap.setUniqueId(UniqueIdGenerator::getUniqueId());
    out_featmap.sortByRT();
  }

  void writeFeatureGroupsInTsvFile(std::vector<FeatureGroup> &fgroups, String infile_path, String outfile_path) const
  {
    std::fstream out_stream;
    out_stream.open(outfile_path, std::fstream::out);

    // header
    out_stream << "FeatureGroupIndex\tFileName\tMonoisotopicMass\tAverageMass\t"
                  "StartRetentionTime(FWHM)\tEndRetentionTime(FWHM)\tHighestApexRetentionTime\tMedianApexRetentionTime\t" // centroid_rt_of_apices
                  "FeatureGroupQuantity\tAllAreaUnderTheCurve\tSumIntensity\tMinCharge\tMaxCharge\tChargeCount\tMostAbundantFeatureCharge\t"
                  "IsotopeCosineScore\tFeatureScore\n"; // mass_trace_ids\n";

//    bool use_smoothed_intensities = FLASHQuantuantAlgorithm().getDefaults().getValue("use_smoothed_intensities").toBool();
    int fg_index = 0;
    for (auto &fg : fgroups)
    {
      // intensities
      double feature_quant = .0; // "bulk" (until 10% of maximum) area under the curve
      double all_area = .0; // all area under the curve

      // centroid rt of apices from all MassTraces
      std::vector<double> apex_rts;
      apex_rts.reserve(fg.size());

      // mass trace labels (ids)
      std::vector<String> mass_trace_labels;
      mass_trace_labels.reserve(fg.size());

      // getting information while looping through mass traces in FeatureGroup
      for (auto &lmt: fg)
      {
        if (lmt.getIsotopeIndex() < 0)
        {
          continue;
        }
        auto &lmt_ptr = lmt.getMassTrace();
        mass_trace_labels.push_back(lmt_ptr.getLabel());

        // find apex
        Size max_idx = lmt_ptr.findMaxByIntPeak(false);
        apex_rts.push_back(lmt_ptr[max_idx].getRT());

//        if (use_smoothed_intensities)
//        {
//          feature_quant += lmt_ptr.computeFwhmAreaSmooth();
//        }
//        else
//        {
//          feature_quant += lmt.computeBulkPeakArea();
//        }
        // calculate bulk area
        feature_quant += lmt.computeBulkPeakArea();

        // to calculate area
        double previous_peak_inty = lmt_ptr[0].getIntensity();
        double previous_peak_rt = lmt_ptr[0].getRT();
        for (auto &peaks: lmt_ptr)
        {
          all_area += (previous_peak_inty + peaks.getIntensity()) / 2 * (peaks.getRT() - previous_peak_rt);
          previous_peak_inty = peaks.getIntensity();
          previous_peak_rt = peaks.getRT();
        }
      }

      // get most abundant charge
      std::vector<float> per_charge_inty = fg.getChargeIntensities();
      int most_abundant_cs = std::distance(per_charge_inty.begin(), std::max_element(per_charge_inty.begin(), per_charge_inty.end()));

      // calculate centroid value
      double centroid_rt_of_apices;
      std::sort(apex_rts.begin(), apex_rts.end());
      Size mts_count = apex_rts.size();
      if (mts_count % 2 == 0) {
        // Find the average of value at index N/2 and (N-1)/2
        centroid_rt_of_apices = (double)(apex_rts[(mts_count-1) / 2] + apex_rts[mts_count / 2]) / 2.0;
      }
      else
      {
        centroid_rt_of_apices = (double) apex_rts[mts_count / 2];
      }

      out_stream << fg_index++ << "\t" << infile_path << "\t"
                 << std::to_string(fg.getMonoisotopicMass()) << "\t" << std::to_string(fg.getAverageMass()) << "\t"
                 << std::to_string(fg.getFwhmRange().first) << "\t" << std::to_string(fg.getFwhmRange().second) << "\t"
                 << std::to_string(fg.getRtOfMostAbundantMT()) << "\t" << std::to_string(centroid_rt_of_apices) << "\t"
                 << std::to_string(feature_quant) << "\t" << std::to_string(all_area) << "\t" << std::to_string(fg.getIntensity()) << "\t"
                 << fg.getMinCharge() << "\t" << fg.getMaxCharge() << "\t" << fg.getChargeSet().size() << "\t" << most_abundant_cs << "\t"
                 << std::to_string(fg.getIsotopeCosine()) << "\t" << std::to_string(fg.getFeatureGroupScore())
                 << std::endl;
      out_stream.flush();
    }
    out_stream.close();
  }

  void writeFeaturSeedsOfFeatureGroupInTsvFile(std::vector<FeatureGroup> &fgroups, String outfile_path) const
  {
    std::fstream out_stream;
    out_stream.open(outfile_path, std::fstream::out);
    out_stream << "FeatureGroupID\tMass\tCharge\tIsotopeIndex\tQuantValue\tCentroidMz\tRTs\tMZs\tIntensities\n"; // header

    Size fg_index = 0;
    for (const auto &fgroup : fgroups)
    {
      for (const auto &trace : fgroup)
      {
        stringstream rts;
        stringstream mzs;
        stringstream intys;
        for (const auto &peak: trace.getMassTrace())
        {
          mzs << std::to_string(peak.getMZ()) << ",";
          rts << std::to_string(peak.getRT()) << ",";
          intys << std::to_string(peak.getIntensity()) << ",";
        }
        std::string peaks = rts.str();
        peaks.pop_back();
        peaks = peaks + "\t" + mzs.str();
        peaks.pop_back();
        peaks = peaks + "\t" + intys.str();
        peaks.pop_back();

        out_stream << fg_index << "\t" << std::to_string(fgroup.getMonoisotopicMass()) << "\t"
                   << trace.getCharge() << "\t" << trace.getIsotopeIndex() << "\t"
                   << std::to_string(trace.getIntensity()) << "\t" << std::to_string(trace.getCentroidMz()) << "\t"
                   << peaks + "\n";
      }
      ++fg_index;
      out_stream.flush();
    }
    out_stream.close();
  }

public:
  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------
    String in = getStringOption_("in");
    String out = getStringOption_("out");
    String out_feat = getStringOption_("out_feat");
    String out_detail = getStringOption_("out_detail");

    MzMLFile mz_data_file;
    mz_data_file.setLogType(log_type_);
    PeakMap ms_peakmap;
    std::vector<Int> ms_level(1, 1);
    mz_data_file.getOptions().setMSLevels(ms_level);
    mz_data_file.load(in, ms_peakmap);

    if (ms_peakmap.empty())
    {
      OPENMS_LOG_WARN << "The given file does not contain any conventional peak data, but might"
                         " contain chromatograms. This tool currently cannot handle them, sorry.";
      return INCOMPATIBLE_INPUT_DATA;
    }
    OPENMS_LOG_INFO << "using " << ms_peakmap.getNrSpectra() << " MS1 spectra" << endl;

    // determine type of spectral data (profile or centroided)
    SpectrumSettings::SpectrumType spectrum_type = ms_peakmap[0].getType();

    if (spectrum_type == SpectrumSettings::PROFILE)
    {
      if (!getFlag_("force"))
      {
        throw OpenMS::Exception::FileEmpty(__FILE__, __LINE__, __FUNCTION__,
                                           "Error: Profile data provided but centroided spectra expected. To enforce processing of the data set the -force flag.");
      }
    }

    // make sure the spectra are sorted by m/z
    ms_peakmap.sortSpectra(true);

    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------
    Param mtd_param = getParam_().copy("algorithm:mtd:", true);
    writeDebug_("Parameters passed to MassTraceDetection", mtd_param, 3);

    Param epd_param = getParam_().copy("algorithm:epd:", true);
    writeDebug_("Parameters passed to ElutionPeakDetection", epd_param, 3);

    Param fdq_param = getParam_().copy("algorithm:fdq:", true);
    writeDebug_("Parameters passed to FLASHQuant", fdq_param, 3);

    //-------------------------------------------------------------
    // Mass traces detection
    //-------------------------------------------------------------
//    Param p_mtd = MassTraceDetection().getDefaults();
//    p_mtd.setValue("noise_threshold_int" , 0.0);
//    p_mtd.setValue("chrom_peak_snr" , 0.0);
//    p_mtd.setValue("mass_error_ppm", 5.0);
//    p_mtd.setValue("trace_termination_criterion", "sample_rate");
//    p_mtd.setValue("min_sample_rate", 0.2);

    vector<MassTrace> m_traces;
    MassTraceDetection mtdet;
    mtdet.setParameters(mtd_param);
    mtdet.run(ms_peakmap, m_traces);
    OPENMS_LOG_INFO << "# initial input mass traces : " << m_traces.size() << endl;

    //-------------------------------------------------------------
    // Elution peak detection
    //-------------------------------------------------------------
//    Param p_epd = ElutionPeakDetection().getDefaults();
//    p_epd.setValue("width_filtering", "off");

    std::vector<MassTrace> m_traces_final;
    ElutionPeakDetection epdet;
    epdet.setParameters(epd_param);
    // fill mass traces with smoothed data as well .. bad design..
    epdet.detectPeaks(m_traces, m_traces_final);

    OPENMS_LOG_INFO << "# final input mass traces : " << m_traces_final.size() << endl;

    //-------------------------------------------------------------
    // Feature finding
    //-------------------------------------------------------------
    FLASHQuantAlgorithm fdq;
    fdq.setParameters(fdq_param);
    std::vector<FeatureGroup> out_fgroups;

    fdq.output_file_path_ = out;
    fdq.run(m_traces_final, out_fgroups);

    //-------------------------------------------------------------
    // writing featureXML output
    //-------------------------------------------------------------

    OPENMS_LOG_INFO << "writing output..." << out << endl;
    writeFeatureGroupsInTsvFile(out_fgroups, in, out);
    if (!out_feat.empty())
    {
      OPENMS_LOG_INFO << "writing output..." << out_feat << endl;

      FeatureMap out_map;
      storeFeatureGroupInOpenMSFeature(out_fgroups, out_map);

      out_map.setPrimaryMSRunPath({in});
      addDataProcessing_(out_map, getProcessingInfo_(DataProcessing::QUANTITATION));
      FeatureXMLFile().store(out_feat, out_map);
    }
    if (!out_detail.empty())
    {
      OPENMS_LOG_INFO << "writing output..." << out_detail << endl;
      writeFeaturSeedsOfFeatureGroupInTsvFile(out_fgroups, out_detail);
    }
    OPENMS_LOG_INFO << "----- output writing done -----" << endl;

    return EXECUTION_OK;
  }
};

int main(int argc, const char** argv)
{
  TOPPFLASHQuant tool;
  return tool.main(argc, argv);
}

/// @endcond
