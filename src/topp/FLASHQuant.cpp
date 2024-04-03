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
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>


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

    registerOutputFile_("out_cons", "<file>", "", "ConsensusXML output file with mass trace information per feature group", false);
    setValidFormats_("out_cons", ListUtils::create<String>("ConsensusXML"));

    registerOutputFile_("out_detail", "<file>", "", "Tsv output file with mass trace information per feature group", false);
    setValidFormats_("out_detail", ListUtils::create<String>("tsv"));

    addEmptyLine_();
    registerSubsection_("algorithm", "Algorithm parameters section");
  }

  Param getSubsectionDefaults_(const String& /*section*/) const override
  {
    Param combined;

    Param p_mtd = MassTraceDetection().getDefaults();
//    p_mtd.setValue("noise_threshold_int", 0.0); // if zero, all peaks are considered for apex
    p_mtd.setValue("chrom_peak_snr", 3.0);
    p_mtd.setValue("mass_error_ppm", 10.0);
    p_mtd.setValue("trace_termination_outliers", 2);
    combined.insert("mtd:", p_mtd);
    combined.setSectionDescription("mtd", "Mass Trace Detection parameters");

    Param p_epd = ElutionPeakDetection().getDefaults();
    p_epd.setValue("width_filtering", "auto");
    p_epd.setValue("smoothing_polynomial", 3);
    p_epd.setValue("min_num_of_peaks", 4);
    combined.insert("epd:", p_epd);
    combined.setSectionDescription("epd", "Elution Profile Detection (to separate isobaric Mass Traces by elution time).");

    Param p_ffi = FLASHQuantAlgorithm().getDefaults();
    combined.insert("fq:", p_ffi);
    combined.setSectionDescription("fq", "FLASHQuant parameters (assembling mass traces to charged features)");

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

  void writeFeatureGroupsInTsvFile(std::vector<FeatureGroup> &fgroups, String infile_path, String outfile_path, bool use_smoothed_intensity) const
  {
    std::fstream out_stream;
    out_stream.open(outfile_path, std::fstream::out);

    // header
    out_stream << "FeatureGroupIndex\tFileName\tMonoisotopicMass\tAverageMass\t"
                  "StartRetentionTime(FWHM)\tEndRetentionTime(FWHM)\tHighestApexRetentionTime\tMedianApexRetentionTime\t" // centroid_rt_of_apices
                  "FeatureGroupQuantity\tAllAreaUnderTheCurve\tSumIntensity\tMinCharge\tMaxCharge\tChargeCount\tMostAbundantFeatureCharge\t"
                  "IsotopeCosineScore\n";

    int fg_index = 0;
    for (auto &fg : fgroups)
    {
      // intensities
      double feature_quant = .0; // "bulk" (until 10% of maximum) area under the curve
      double all_area = .0; // all area under the curve

      // centroid rt of apices from all MassTraces
      std::vector<double> apex_rts;
      apex_rts.reserve(fg.size());
//      std::pair<double, double> medians_of_fwhms = fg.getMedianValuesOfFWHMs();

      // getting information while looping through mass traces in FeatureGroup
      for (auto &lmt: fg)
      {
        if (lmt.getIsotopeIndex() < 0)
        {
          continue;
        }
        auto &lmt_ptr = lmt.getMassTrace();

        // find apex
        Size max_idx = lmt_ptr.findMaxByIntPeak(use_smoothed_intensity);
        apex_rts.push_back(lmt_ptr[max_idx].getRT());

        // calculate bulk area
        feature_quant += lmt.computeBulkPeakArea(use_smoothed_intensity);

        // to calculate area
        if (use_smoothed_intensity)
        {
          auto smoothed_inty = lmt_ptr.getSmoothedIntensities();
          double previous_peak_inty = smoothed_inty[0];
          double previous_peak_rt = lmt_ptr[0].getRT();
          for (Size i = 0; i < lmt_ptr.getSize(); ++i)
          {
            all_area += (previous_peak_inty + smoothed_inty[i]) / 2 * (lmt_ptr[i].getRT() - previous_peak_rt);
            previous_peak_inty = smoothed_inty[i];
            previous_peak_rt = lmt_ptr[i].getRT();
          }
        }
        else
        {
          double previous_peak_inty = lmt_ptr[0].getIntensity();
          double previous_peak_rt = lmt_ptr[0].getRT();
          for (auto &peaks: lmt_ptr)
          {
//            if ((previous_peak_rt >= medians_of_fwhms.first) && (peaks.getRT() <= medians_of_fwhms.second))
//            {
//              feature_quant += (previous_peak_inty + peaks.getIntensity()) / 2 * (peaks.getRT() - previous_peak_rt);
//            }
            all_area += (previous_peak_inty + peaks.getIntensity()) / 2 * (peaks.getRT() - previous_peak_rt);
            previous_peak_inty = peaks.getIntensity();
            previous_peak_rt = peaks.getRT();
          }
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
                 << std::to_string(fg.getIsotopeCosine())
                 << std::endl;
      out_stream.flush();
    }
    out_stream.close();
  }

  void writeFeaturSeedsOfFeatureGroupInTsvFile(std::vector<FeatureGroup> &fgroups, String outfile_path) const
  {
    std::fstream out_stream;
    out_stream.open(outfile_path, std::fstream::out);
    out_stream << "FeatureGroupID\tMass\tCharge\tIsotopeIndex\tQuantValue\tCentroidMz\tisTheoretical\tTraceLabel\tRTs\tMZs\tIntensities\n"; // header

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
                   << 0 << "\t" << trace.getMassTrace().getLabel() << "\t" << peaks + "\n";
      }

      // theoretical
      for (const auto &shape : fgroup.getTheoreticalShapes())
      {
        stringstream rts;
        stringstream mzs;
        stringstream intys;
        for (const auto &peak: shape.getMassTrace())
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
                   << shape.getCharge() << "\t" << shape.getIsotopeIndex() << "\t"
                   << std::to_string(shape.getIntensity()) << "\t" << std::to_string(shape.getCentroidMz()) << "\t"
                   << 1 << "\t" << shape.getMassTrace().getLabel() << "\t" << peaks + "\n";
      }

      ++fg_index;
      out_stream.flush();
    }
    out_stream.close();
  }

  void storeFeatureGroupInConsensusMap(std::vector<FeatureGroup> &fgroups, ConsensusMap& consensus_map) const
  {
    for (Size fg_index = 0; fg_index < fgroups.size(); ++fg_index)
    {
      // traces
      insertTracesInConsensusFeature(fgroups[fg_index].getSeeds(), fg_index, consensus_map);

      // theoretical shape
      insertTracesInConsensusFeature(fgroups[fg_index].getTheoreticalShapes(), fg_index, consensus_map, true);
    }
  }

  void insertTracesInConsensusFeature(vector<FeatureSeed> const& seed_group, Size fg_index, ConsensusMap& out_map, bool isTheoreticalShape = false) const
  {
    for (const auto &trace: seed_group)
    {
      ConsensusFeature fcons;
      int k = 0;
      for (const Peak2D& peak : trace.getMassTrace())
      {
        FeatureHandle fhandle;
        fhandle.setRT(peak.getRT());
        fhandle.setMZ(peak.getMZ());
        fhandle.setIntensity(peak.getIntensity());
        fhandle.setUniqueId(++k);
        fcons.insert(fhandle);
      }

      fcons.setMetaValue("FeatureGroupIndex", fg_index);
      fcons.setMetaValue("Mass", trace.getMass());
      fcons.setMetaValue("IsotopeIndex", trace.getIsotopeIndex());
      if (isTheoreticalShape)
      {
        fcons.setMetaValue("theoretical_shape", 1);
      }

      fcons.setCharge(trace.getCharge());
      fcons.setWidth(trace.getFwhmEnd() - trace.getFwhmStart());
      fcons.setQuality(1 - (1.0 / trace.getMassTrace().getSize()));

      fcons.setRT(trace.getCentroidRT());
      fcons.setMZ(trace.getCentroidMz());
      fcons.setIntensity(trace.getIntensity());
      out_map.push_back(fcons);
    }
  }

  double getNoiseIntensityFromTheInitialScans(PeakMap &map, Size num_of_scans = 3)
  {
    std::vector<double> median_intensities;
    for (PeakMap::const_iterator iter=map.begin(); iter != map.begin()+num_of_scans; ++iter)
    {
      std::vector<double> intensities;
      intensities.reserve(iter->size());
      for (auto &peaks : *iter)
      {
        intensities.push_back(peaks.getIntensity());
      }
      median_intensities.push_back(OpenMS::Math::median(intensities.begin(), intensities.end()));
    }
    return OpenMS::Math::median(median_intensities.begin(), median_intensities.end());
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
    String out_cons = getStringOption_("out_cons");
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
    // get rt information from the original spectra
    std::vector<double> rts_from_original_spec;
    rts_from_original_spec.reserve(ms_peakmap.getNrSpectra());
    for (const auto& scan : ms_peakmap)
    {
      rts_from_original_spec.push_back(scan.getRT());
    }

    // noise estimation
//    double noise = getNoiseIntensityFromTheInitialScans(ms_peakmap, 3);
//    cout << "calculated noise: " + to_string(noise) << endl;

    // make sure the spectra are sorted by m/z
    ms_peakmap.sortSpectra(true);

    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------
    Param mtd_param = getParam_().copy("algorithm:mtd:", true);
    writeDebug_("Parameters passed to MassTraceDetection", mtd_param, 3);

    Param epd_param = getParam_().copy("algorithm:epd:", true);
    writeDebug_("Parameters passed to ElutionPeakDetection", epd_param, 3);

    Param fq_param = getParam_().copy("algorithm:fq:", true);
    writeDebug_("Parameters passed to FLASHQuant", fq_param, 3);

    //-------------------------------------------------------------
    // Mass traces detection
    //-------------------------------------------------------------
    std::vector<MassTrace> m_traces;
    MassTraceDetection mtdet;
//    mtd_param.setValue("noise_threshold_int", noise);
    mtdet.setParameters(mtd_param);
    mtdet.run(ms_peakmap, m_traces);
    OPENMS_LOG_INFO << "# initial input mass traces : " << m_traces.size() << endl;

    //-------------------------------------------------------------
    // Elution peak detection
    //-------------------------------------------------------------
    std::vector<MassTrace> m_traces_final;
    ElutionPeakDetection epdet;
    epdet.setParameters(epd_param);
    // fill mass traces with smoothed data as well .. bad design..
    epdet.detectPeaks(m_traces, m_traces_final);
//    std::vector<MassTrace> tmp_traces (m_traces.begin()+500000, m_traces.begin()+500100);
//    epdet.detectPeaks(tmp_traces, m_traces_final);

    OPENMS_LOG_INFO << "# mass traces after elution peak detection : " << m_traces_final.size() << endl;

    //-------------------------------------------------------------
    // Feature finding
    //-------------------------------------------------------------
    FLASHQuantAlgorithm fq_algo;
    fq_algo.setParameters(fq_param);
    std::vector<FeatureGroup> out_fgroups;

    fq_algo.output_file_path_ = out;
    fq_algo.rts_from_org_scans_ = rts_from_original_spec;
    fq_algo.run(m_traces_final, out_fgroups);

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    bool use_smoothed_intensity = getParam_().getValue("algorithm:fq:use_smoothed_intensities").toBool();
    OPENMS_LOG_INFO << "writing output..." << out << endl;
    writeFeatureGroupsInTsvFile(out_fgroups, in, out, use_smoothed_intensity);
    if (!out_feat.empty())
    {
      OPENMS_LOG_INFO << "writing output..." << out_feat << endl;

      FeatureMap out_map;
      storeFeatureGroupInOpenMSFeature(out_fgroups, out_map);

      out_map.setPrimaryMSRunPath({in});
      addDataProcessing_(out_map, getProcessingInfo_(DataProcessing::QUANTITATION));
      FeatureXMLFile().store(out_feat, out_map);
    }
    if (!out_cons.empty())
    {
      OPENMS_LOG_INFO << "writing output..." << out_cons << endl;
      ConsensusMap out_consensus_map;
      out_consensus_map.reserve(m_traces_final.size());
      out_consensus_map.setPrimaryMSRunPath({in}, ms_peakmap);
      storeFeatureGroupInConsensusMap(out_fgroups, out_consensus_map);

      out_consensus_map.applyMemberFunction(&UniqueIdInterface::setUniqueId);
      addDataProcessing_(out_consensus_map, getProcessingInfo_(DataProcessing::QUANTITATION));
      ConsensusXMLFile().store(out_cons, out_consensus_map);
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
