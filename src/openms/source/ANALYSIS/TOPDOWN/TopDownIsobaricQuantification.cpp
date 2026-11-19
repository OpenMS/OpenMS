// Copyright (c) 2002-2025, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
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
#include <OpenMS/ANALYSIS/TOPDOWN/TopDownIsobaricQuantification.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/METADATA/SpectrumLookup.h>
#include <cmath>

namespace OpenMS
{
  TopDownIsobaricQuantification::TopDownIsobaricQuantification() : DefaultParamHandler("TopDownIsobaricQuantification")
  {
    setDefaultParams_();
  }

  TopDownIsobaricQuantification::TopDownIsobaricQuantification(const TopDownIsobaricQuantification& other) : DefaultParamHandler(other)
  {
  }

  TopDownIsobaricQuantification& TopDownIsobaricQuantification::operator=(const TopDownIsobaricQuantification& rhs)
  {
    if (this == &rhs)
      return *this;

    DefaultParamHandler::operator=(rhs);
    return *this;
  }

  void TopDownIsobaricQuantification::setDefaultParams_()
  {
    defaults_.setValue("type", "none", "Specifies the isobaric quantification method used in the experiment.");
    defaults_.setValidStrings("type", {"none", "itraq4plex", "itraq8plex", "tmt10plex", "tmt11plex", "tmt16plex", "tmt18plex", "tmt6plex"});
    defaults_.setValue("isotope_correction", "true",
                       "Enable isotope correction (highly recommended).");
    defaults_.setValidStrings("isotope_correction", {"true", "false"});
    defaults_.setValue("correction_matrix", std::vector<std::string>(),
                       "The reporter ion isotope impurities of the reagent lot, used by isotope_correction: one entry per channel, in the order of "
                       "the channels, each entry <-2Da>/<-1Da>/<+1Da>/<+2Da> in percent as on the product sheet of the lot (e.g. 0.0/0.0/8.6/0.3). "
                       "Empty (default): the matrix built into the method of 'type' is used (tmt6plex: TMTsixplex lot ZE386964). "
                       "Given: it replaces the built-in matrix. It needs as many entries as the method has channels and the number of values per "
                       "entry of the built-in matrix; otherwise the tool stops with an error before any quantification. "
                       "A reagent that shares the reporter masses of a method but not its impurities (e.g. iodoTMTsixplex with tmt6plex) needs "
                       "the values of its own lot here, or isotope_correction false.");
    defaults_.setValue("reporter_mz_tol", 2e-3, "Specifies the m/z tolerance (in Th) for reporter ion detection.");
    defaults_.setValue("only_fully_quantified", "false", "Restricts analysis to spectra that are fully quantified, meaning all channels have non-zero intensity reporter ions.");
    defaults_.setValidStrings("only_fully_quantified", {"true", "false"});
    defaultsToParam_();
  }

  void TopDownIsobaricQuantification::updateMembers_()
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

  void TopDownIsobaricQuantification::quantify(const MSExperiment& exp, std::vector<DeconvolvedSpectrum>& deconvolved_spectra, const std::vector<FLASHHelperClasses::MassFeature>& mass_features)
  {
    // set the parameters for this method
    std::string type = getParameters().getValue("type").toString();

    if (type == "none")
    {
      return;
    }

    const auto& quant_method = quant_methods_[type];

    // the impurities of the user's reagent lot instead of the built-in ones
    const StringList correction_matrix = ListUtils::toStringList<std::string>(getParameters().getValue("correction_matrix"));
    if (! correction_matrix.empty())
    {
      if (correction_matrix.size() != quant_method->getNumberOfChannels())
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "correction_matrix has " + std::to_string(correction_matrix.size()) + " entries, but the method '" + type + "' has "
                                            + std::to_string(quant_method->getNumberOfChannels()) + " channels (one entry per channel is needed).");
      }
      Param method_param = quant_method->getParameters();
      method_param.setValue("correction_matrix", correction_matrix);
      quant_method->setParameters(method_param);
      quant_method->getIsotopeCorrectionMatrix(); // stops here on a malformed entry
    }

    // The channel extractor needs the precursor of every spectrum it quantifies. It stops with an exception on an MSn (n > 2) spectrum
    // whose MS(n-1) spectrum is not in the file (e.g. a file of MS3 spectra only), and it takes the first precursor of a spectrum
    // without a check that there is one. Such spectra are handled on a copy of the experiment: an MSn spectrum without its MS(n-1)
    // spectrum is passed as an MS2 spectrum, so that its reporter ions are quantified in the spectrum itself and assigned to its own
    // scan (in a file that also holds complete MSn spectra the extractor quantifies only those, as before); a spectrum above MS1
    // without a precursor, and an MSn spectrum whose MS(n-1) spectrum has none, are left out. A complete file is passed as it is.
    const MSExperiment* quant_exp = &exp;
    MSExperiment adapted_exp;
    {
      std::vector<bool> left_out(exp.size(), false), orphan(exp.size(), false);
      Size left_out_count = 0, orphan_count = 0;
      for (auto it = exp.begin(); it != exp.end(); ++it)
      {
        if (it->getMSLevel() < 2) continue;
        const Size i = (Size)std::distance(exp.begin(), it);
        if (it->getPrecursors().empty()) { left_out[i] = true; }
        else if (it->getMSLevel() > 2)
        {
          const auto parent = exp.getPrecursorSpectrum(it);
          if (parent == exp.end()) { orphan[i] = true; }
          else if (parent->getPrecursors().empty()) { left_out[i] = true; }
        }
        left_out_count += left_out[i];
        orphan_count += orphan[i];
      }
      if (left_out_count + orphan_count > 0)
      {
        static_cast<ExperimentalSettings&>(adapted_exp) = exp;
        adapted_exp.reserveSpaceSpectra(exp.size() - left_out_count);
        for (Size i = 0; i < exp.size(); i++)
        {
          if (left_out[i]) continue;
          adapted_exp.addSpectrum(exp[i]);
          if (orphan[i]) adapted_exp[adapted_exp.size() - 1].setMSLevel(2);
        }
        quant_exp = &adapted_exp;
        if (orphan_count > 0)
        {
          OPENMS_LOG_WARN << "Isobaric quantification: " << orphan_count
                          << " MSn (n > 2) spectra without their MS(n-1) spectrum in the file are quantified by their own reporter ions.\n";
        }
        if (left_out_count > 0)
        {
          OPENMS_LOG_WARN << "Isobaric quantification: " << left_out_count << " spectra without precursor information are not quantified.\n";
        }
      }
    }

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
    channel_extractor.extractChannels(*quant_exp, consensus_map_raw);
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
        if (! it->getPrecursors().empty()) ms2_scan_precursor_mz[scan_number] = it->getPrecursors()[0].getMZ(); // an MSn without precursor: none
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

    // the spectrum of the retention time nearest to rt (end() if there is none): lower_bound alone was dereferenced past the last spectrum
    const auto nearest_rt = [&rt_scan_map](double rt) {
      auto it = rt_scan_map.lower_bound(rt);
      if (it == rt_scan_map.begin()) return it;
      auto prev = std::prev(it);
      if (it == rt_scan_map.end() || std::abs(prev->first - rt) < std::abs(it->first - rt)) return prev;
      return it;
    };

    for (auto& mf : mass_features)
    {
      auto& mass_trace = mf.mt;

      std::vector<PeakGroup> cluster;
      cluster.reserve(mass_trace.getSize());

      // each peak = a precursor peak group
      for (auto& p : mass_trace)
      {
        const auto trt_it = nearest_rt(p.getRT());
        if (trt_it == rt_scan_map.end() || std::abs(trt_it->first - p.getRT()) > .01)
          continue;
        int scan = trt_it->second;
        if (scan_precursors_map.find(scan) == scan_precursors_map.end())
          continue;
        for (auto& pg : scan_precursors_map[scan])
        {
          if (std::abs(pg.getMonoMass() - p.getMZ()) > .01)
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
        precursor_cluster_index[pg] = (int)precursor_clusters.size() - 1;
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
      precursor_cluster_index[precursor] = (int)precursor_clusters.size() - 1;
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

      const auto trt_it = nearest_rt(feature.getRT());
      if (trt_it == rt_scan_map.end() || std::abs(trt_it->first - feature.getRT()) > .01)
        continue;
      int scan = trt_it->second;
      ms2_ints[scan] = intensities;
    }

    intensity_clusters.resize(precursor_clusters.size());
    merged_intensity_clusters.resize(precursor_clusters.size());

    for (auto& dspec : deconvolved_spectra)
    {
      if (dspec.getOriginalSpectrum().getMSLevel() == 1)
      {
        continue;
      }

      int scan = dspec.getScanNumber();
      double pre_mz = ms2_scan_precursor_mz[scan];
      std::vector<double> intensities (0);
      for (int ms2_scan : precursor_scan_ms2_scans[ms2_scan_precursor_scan[scan]])
      {
        if (ms2_ints.find(ms2_scan) == ms2_ints.end() || ms2_ints[ms2_scan].empty())
          continue;
        if (ms2_scan_precursor_mz.find(ms2_scan) == ms2_scan_precursor_mz.end() || std::abs(ms2_scan_precursor_mz[ms2_scan] - pre_mz) > .01)
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

      int cluster_index = dspec.getPrecursorPeakGroup().empty()? -1 : precursor_cluster_index[dspec.getPrecursorPeakGroup()];
      double min_intensity = only_fully_quantified_? *std::min_element(intensities.begin(), intensities.end()) : *std::max_element(intensities.begin(), intensities.end());
      if (min_intensity > 0) // at least one channel quantified
      {
        if (cluster_index >= 0) intensity_clusters[cluster_index].push_back(intensities);
        FLASHHelperClasses::IsobaricQuantities iq;
        iq.scan = scan;
        iq.quantities = intensities;
        iq.merged_quantities = intensities;
        dspec.setQuantities(iq);
      }
    }

    for (Size i = 0; i < intensity_clusters.size(); i++)
    {
      const auto& intensities = intensity_clusters[i];
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
        FLASHHelperClasses::IsobaricQuantities iq = dspec.getQuantities();
        iq.merged_quantities = intensities;
        dspec.setQuantities(iq);
      }
    }
  }
} // namespace OpenMS