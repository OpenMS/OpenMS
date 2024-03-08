// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricChannelExtractor.h>
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantifier.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqEightPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqFourPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTEighteenPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTElevenPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTSixPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTSixteenPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTTenPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/TOPDOWN/TopDownIsobaricQuantifier.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/METADATA/SpectrumLookup.h>

namespace OpenMS
{
  TopDownIsobaricQuantifier::TopDownIsobaricQuantifier() : DefaultParamHandler("TopDownIsobaricQuantifier")
  {
    setDefaultParams_();
  }

  TopDownIsobaricQuantifier::TopDownIsobaricQuantifier(const TopDownIsobaricQuantifier& other) : DefaultParamHandler(other)
  {
  }

  TopDownIsobaricQuantifier& TopDownIsobaricQuantifier::operator=(const TopDownIsobaricQuantifier& rhs)
  {
    if (this == &rhs)
      return *this;

    DefaultParamHandler::operator=(rhs);
    return *this;
  }

  void TopDownIsobaricQuantifier::setDefaultParams_()
  {
    defaults_.setValue("type", "none", "Isobaric Quantitation method used in the experiment.");
    defaults_.setValidStrings("type", {"none", "itraq4plex", "itraq8plex", "tmt10plex", "tmt11plex", "tmt16plex", "tmt18plex", "tmt6plex"});
    defaults_.setValue("isotope_correction", "true",
                       "Enable isotope correction (highly recommended). "
                       "Note that you need to provide a correct isotope correction matrix "
                       "otherwise the tool will fail or produce invalid results.");
    defaults_.setValidStrings("isotope_correction", {"true", "false"});
    defaults_.setValue("reporter_mz_tol", 2e-3, "m/z tolerance in Th from the expected position of reporter ion m/zs.");
    defaults_.setValue("only_fully_quantified", "false", "Use only the fully quantified spectra in which non-zero intensity report ions are found for all channels.");
    defaults_.setValidStrings("only_fully_quantified", {"true", "false"});
    defaultsToParam_();
  }

  void TopDownIsobaricQuantifier::updateMembers_()
  {
    addMethod_(std::make_unique<ItraqFourPlexQuantitationMethod>());
    addMethod_(std::make_unique<ItraqEightPlexQuantitationMethod>());
    addMethod_(std::make_unique<TMTSixPlexQuantitationMethod>());
    addMethod_(std::make_unique<TMTTenPlexQuantitationMethod>());
    addMethod_(std::make_unique<TMTElevenPlexQuantitationMethod>());
    addMethod_(std::make_unique<TMTSixteenPlexQuantitationMethod>());
    addMethod_(std::make_unique<TMTEighteenPlexQuantitationMethod>());
    only_fully_quantified_ = param_.getValue("only_fully_quantified").toString() == "true";
  }

  void TopDownIsobaricQuantifier::quantify(const MSExperiment& exp, std::vector<DeconvolvedSpectrum>& deconvolved_spectra, const std::vector<FLASHDeconvHelperStructs::MassFeature>& mass_features)
  {
    // set the parameters for this method
    // Param extract_param(getParam_().copy("extraction:", true));
    String type = getParameters().getValue("type").toString();

    if (type == "none")
    {
      return;
    }

    const auto& quant_method = quant_methods_[type];

    IsobaricChannelExtractor channel_extractor(quant_method.get());
    Param extract_param = channel_extractor.getDefaults();
    extract_param.setValue("reporter_mass_shift", getParameters().getValue("reporter_mz_tol"));
    channel_extractor.setParameters(extract_param);

    IsobaricQuantifier quantifier(quant_method.get());
    Param quant_param = quantifier.getDefaults();
    quant_param.setValue("isotope_correction", getParameters().getValue("isotope_correction"));
    quant_param.setValue("normalization", "false"); // here use its own normalization for the same precursor masses.
    quantifier.setParameters(quant_param);

    ConsensusMap consensus_map_raw, consensus_map_quant;

    // extract channel information
    channel_extractor.extractChannels(exp, consensus_map_raw);
    quantifier.quantify(consensus_map_raw, consensus_map_quant);

    std::map<double, int> rt_scan_map;
    std::map<int, std::set<PeakGroup>> scan_precursors_map; // MS1 scan to precursor peak groups
    std::vector<std::vector<PeakGroup>> precursor_clusters; // clusters of the precursor peak groups
    std::vector<std::vector<std::vector<double>>> intensity_clusters;
    std::vector<std::vector<double>> merged_intensity_clusters;
    std::map<PeakGroup, int> precursor_cluster_index; // precursor to cluster index

    std::map<int, std::vector<int>> precursor_scan_ms2_scans; // from precursor scan to ms2 scans
    std::map<int, int> ms2_scan_precursor_scan;               // from ms2 scan to precursor scan
    std::map<int, double> ms2_scan_precursor_mz;              // from ms2 scan to precursor mz

    int pre_scan = 0;
    for (auto it = exp.begin(); it != exp.end(); ++it)
    {
      int scan_number = exp.getSourceFiles().empty() ? -1 : SpectrumLookup::extractScanNumber(it->getNativeID(), exp.getSourceFiles()[0].getNativeIDTypeAccession());

      if (scan_number < 0)
      {
        scan_number = (int)std::distance(exp.begin(), it) + 1;
      }
      rt_scan_map[it->getRT()] = scan_number;

      if (it->getMSLevel() == 1)
      {
        pre_scan = scan_number;
        precursor_scan_ms2_scans[pre_scan] = std::vector<int>();
      }
      else
      {
        precursor_scan_ms2_scans[pre_scan].push_back(scan_number);
        ms2_scan_precursor_scan[scan_number] = pre_scan;
        ms2_scan_precursor_mz[scan_number] = it->getPrecursors()[0].getMZ();
      }
    }

    for (auto& dspec : deconvolved_spectra)
    {
      if (dspec.getOriginalSpectrum().getMSLevel() == 1)
      {
        continue;
      }
      auto& precursor = dspec.getPrecursorPeakGroup();
      if (precursor.empty())
        continue;
      scan_precursors_map[precursor.getScanNumber()].insert(precursor);
    }

    for (auto& mf : mass_features)
    {
      auto& mass_trace = mf.mt;

      std::vector<PeakGroup> cluster;
      cluster.reserve(mass_trace.getSize());

      // each peak = a precursor peak group
      for (auto& p : mass_trace)
      {
        auto trt = *rt_scan_map.lower_bound(p.getRT());
        if (abs(trt.first - p.getRT()) > .01)
          continue;
        int scan = trt.second;
        if (scan_precursors_map.find(scan) == scan_precursors_map.end())
          continue;
        for (auto& pg : scan_precursors_map[scan])
        {
          if (abs(pg.getMonoMass() - p.getMZ()) > .01)
            continue;
          cluster.push_back(pg);
        }
      }

      if (!cluster.empty())
      {
        precursor_clusters.push_back(cluster);
      }

      for (auto& pg : cluster)
      {
        precursor_cluster_index[pg] = precursor_clusters.size() - 1;
      }
    }

    for (auto& dspec : deconvolved_spectra)
    {
      if (dspec.getOriginalSpectrum().getMSLevel() == 1)
      {
        continue;
      }
      auto& precursor = dspec.getPrecursorPeakGroup();
      if (precursor.empty() || precursor_cluster_index.find(precursor) != precursor_cluster_index.end())
        continue;
      precursor_clusters.push_back(std::vector<PeakGroup> {precursor});
      precursor_cluster_index[precursor] = precursor_clusters.size() - 1;
    }

    std::map<int, std::vector<double>> ms2_ints; // from ms2 scan to intensities

    for (auto& feature : consensus_map_quant)
    {
      std::vector<double> intensities;
      float max_int = 0;
      for (auto& i : feature)
      {
        max_int = std::max(max_int, i.getIntensity());
        intensities.push_back(i.getIntensity());
      }

      if (max_int <= 0)
        continue;

      auto trt = *rt_scan_map.lower_bound(feature.getRT());
      if (abs(trt.first - feature.getRT()) > .01)
        continue;
      int scan = trt.second;
      ms2_ints[scan] = intensities;
    }

    intensity_clusters.resize(precursor_clusters.size());
    merged_intensity_clusters.resize(precursor_clusters.size());

    for (auto& dspec : deconvolved_spectra)
    {
      if (dspec.getOriginalSpectrum().getMSLevel() == 1 || dspec.getPrecursorPeakGroup().empty())
      {
        continue;
      }
      int cluster_index = precursor_cluster_index[dspec.getPrecursorPeakGroup()];
      int scan = dspec.getScanNumber();
      double pre_mz = ms2_scan_precursor_mz[scan];
      std::vector<double> intensities (0);
      for (int ms2_scan : precursor_scan_ms2_scans[ms2_scan_precursor_scan[scan]])
      {
        if (ms2_ints.find(ms2_scan) == ms2_ints.end() || ms2_ints[ms2_scan].empty())
          continue;
        if (ms2_scan_precursor_mz.find(ms2_scan) == ms2_scan_precursor_mz.end() || abs(ms2_scan_precursor_mz[ms2_scan] - pre_mz) > .01)
          continue;
        if (intensities.empty())
        {
          intensities = ms2_ints[ms2_scan];
        }
        else
        {
          for (Size j = 0; j < intensities.size(); j++)
            intensities[j] += ms2_ints[ms2_scan][j];
        }
      }

      if (intensities.empty()) continue;
      double min_intensity = only_fully_quantified_? *std::min_element(intensities.begin(), intensities.end()) : *std::max_element(intensities.begin(), intensities.end());
      if (min_intensity > 0) // at least one channel quantified
      {
        intensity_clusters[cluster_index].push_back(intensities);
        FLASHDeconvHelperStructs::IsobaricQuantities iq;
        iq.scan = scan;
        iq.quantities = intensities;
        iq.merged_quantities = intensities;
        dspec.setQuantities(iq);
      }
    }

    for (Size i = 0; i < intensity_clusters.size(); i++)
    {
      auto intensities = intensity_clusters[i];
      if (intensities.empty())
        continue;
      merged_intensity_clusters[i] = intensities[0];
      for (Size j = 1; j < intensities.size(); j++)
      {
        for (Size k = 0; k < merged_intensity_clusters[i].size(); k++)
        {
          merged_intensity_clusters[i][k] += intensities[j][k];
        }
      }
    }
    for (auto& dspec : deconvolved_spectra)
    {
      if (dspec.getOriginalSpectrum().getMSLevel() == 1 || dspec.getPrecursorPeakGroup().empty())
        continue;
      int cluster_index = precursor_cluster_index[dspec.getPrecursorPeakGroup()];
      if (merged_intensity_clusters[cluster_index].empty())
        continue;

      if (dspec.getQuantities().empty())
        continue;
      auto intensities = merged_intensity_clusters[cluster_index];
      if (intensities.empty()) continue;
      double min_intensity = only_fully_quantified_? *std::min_element(intensities.begin(), intensities.end()) : *std::max_element(intensities.begin(), intensities.end());
      if (min_intensity > 0) // all channel quantified
      {
        FLASHDeconvHelperStructs::IsobaricQuantities iq = dspec.getQuantities();
        iq.merged_quantities = intensities;
        dspec.setQuantities(iq);
      }
      //else quantities_.erase(scan);
    }
  }
} // namespace OpenMS