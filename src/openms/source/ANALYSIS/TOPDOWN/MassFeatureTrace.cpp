//--------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/MassFeatureTrace.h>

#include <utility>

namespace OpenMS
{
  MassFeatureTrace::MassFeatureTrace() :
      DefaultParamHandler("MassFeatureTrace")
  {
    Param mtd_defaults = MassTraceDetection().getDefaults();
    mtd_defaults.setValue("min_sample_rate",
                          .1,
                          "Minimum fraction of scans along the feature trace that must contain a peak. To raise feature detection sensitivity, lower this value close to 0.");
    mtd_defaults.setValue("min_trace_length", 10.0);

    mtd_defaults.setValue("chrom_peak_snr", .0);
    mtd_defaults.addTag("chrom_peak_snr", "advanced");
    mtd_defaults.setValue("reestimate_mt_sd", "false");
    mtd_defaults.addTag("reestimate_mt_sd", "advanced");
    mtd_defaults.setValue("noise_threshold_int", .0);
    mtd_defaults.addTag("noise_threshold_int", "advanced");

    mtd_defaults.setValue("quant_method", "area");
    mtd_defaults.addTag("quant_method", "advanced"); // hide entry

    defaults_.insert("", mtd_defaults);
    defaults_.setValue("min_isotope_cosine", .75, "cosine threshold between avg. and observed isotope pattern for MS1");
    defaultsToParam_();
  }

  std::vector<FLASHDeconvHelperStructs::MassFeature> MassFeatureTrace::findFeatures(const PrecalculatedAveragine& averagine)
  {
    MSExperiment map;
    std::map<int, MSSpectrum> index_spec_map;
    int min_abs_charge = INT_MAX;
    int max_abs_charge = INT_MIN;
    bool is_positive = true;
    std::vector<FLASHDeconvHelperStructs::MassFeature> mass_features;
    for (auto& item: peak_group_map_)
    {
      double rt = item.first;
      MSSpectrum deconv_spec;
      deconv_spec.setRT(rt);
      for (auto& pg: item.second)
      {
        is_positive = pg.second.isPositive();
        auto crange = pg.second.getAbsChargeRange();
        max_abs_charge = max_abs_charge > std::get<1>(crange) ? max_abs_charge : std::get<1>(crange);
        min_abs_charge = min_abs_charge < std::get<0>(crange) ? min_abs_charge : std::get<0>(crange);

        Peak1D tp(pg.first, (float) pg.second.getIntensity());
        deconv_spec.push_back(tp);
      }
      map.addSpectrum(deconv_spec);
    }

    // when map size is less than 3, MassTraceDetection aborts - too few spectra for mass tracing.
    if (map.size() < 3)
    {
      return mass_features;
    }

    MassTraceDetection mtdet;
    Param mtd_param = getParameters().copy("");
    mtd_param.remove("min_isotope_cosine");

    mtdet.setParameters(mtd_param);
    std::vector<MassTrace> m_traces;

    mtdet.run(map, m_traces);  // m_traces : output of this function
    int charge_range = max_abs_charge - min_abs_charge + 1;

    for (auto& mt: m_traces)
    {
      double max_qscore = .0;
      int min_feature_abs_charge = INT_MAX; // min feature charge
      int max_feature_abs_charge = INT_MIN; // max feature charge

      auto per_isotope_intensity = std::vector<float>(averagine.getMaxIsotopeIndex(), .0f);
      auto per_charge_intensity = std::vector<float>(charge_range + min_abs_charge + 1, .0f);

      double max_iso = 0;
      boost::dynamic_bitset<> charges(charge_range + 1);
      std::vector<PeakGroup> pgs;
      pgs.reserve(mt.getSize());

      for (auto& p2: mt)
      {
        auto& pg_map = peak_group_map_[p2.getRT()];
        auto& pg = pg_map[p2.getMZ()];
        auto crange = pg.getAbsChargeRange();

        min_feature_abs_charge =
            min_feature_abs_charge < std::get<0>(crange) ? min_feature_abs_charge : std::get<0>(crange);
        max_feature_abs_charge =
            max_feature_abs_charge > std::get<1>(crange) ? max_feature_abs_charge : std::get<1>(crange);

        if (pg.getIsotopeCosine() > max_iso)
        {
          max_iso = pg.getIsotopeCosine();
        }
        pgs.push_back(pg);
      }

      for(auto& pg: pgs)
      {
        for (auto& p: pg)
        {
          if (p.isotopeIndex < 0 || p.isotopeIndex >= averagine.getMaxIsotopeIndex() || p.abs_charge < min_abs_charge ||
              p.abs_charge >= charge_range + min_abs_charge + 1)
          {
            continue;
          }

          charges[p.abs_charge - min_abs_charge] = true;
          per_charge_intensity[p.abs_charge] += (p.intensity);
          per_isotope_intensity[p.isotopeIndex] += (p.intensity);
        }

        max_qscore = max_qscore < pg.getQScore() ? pg.getQScore() : max_qscore;
      }

      int offset = 0;
      double mass = mt.getCentroidMZ();
      double isotope_score = FLASHDeconvAlgorithm::getIsotopeCosineAndDetermineIsotopeIndex(mass,
                                                                                            per_isotope_intensity,
                                                                                            offset, averagine, 0);

      if (isotope_score < min_isotope_cosine_)
      {
        continue;
      }

      FLASHDeconvHelperStructs::MassFeature mass_feature;
      mass_feature.iso_offset = offset;
      mass += offset * Constants::ISOTOPE_MASSDIFF_55K_U;

      mass_feature.avg_mass = averagine.getAverageMassDelta(mass) + mass;
      mass_feature.mt = mt;
      mass_feature.charge_count = charges.count();
      mass_feature.isotope_score = isotope_score;
      mass_feature.min_charge = (is_positive ? min_feature_abs_charge : -max_feature_abs_charge) ;
      mass_feature.max_charge = (is_positive ? max_feature_abs_charge : -min_feature_abs_charge) ;
      mass_feature.qscore = max_qscore;

      if(offset != 0)
      {
        per_isotope_intensity = std::vector<float>(averagine.getMaxIsotopeIndex(), .0f); // recalculate with updated offset
        per_charge_intensity = std::vector<float>(charge_range + min_abs_charge + 1, .0f);

        for (auto& pg : pgs)
        {
          for (auto& p : pg)
          {
            if (p.isotopeIndex < offset || p.isotopeIndex >= averagine.getMaxIsotopeIndex() + offset || p.abs_charge < min_abs_charge || p.abs_charge >= charge_range + min_abs_charge + 1)
            {
              continue;
            }

            charges[p.abs_charge - min_abs_charge] = true;
            per_charge_intensity[p.abs_charge] += (p.intensity);
            per_isotope_intensity[p.isotopeIndex - offset] += (p.intensity);
          }
        }
      }

      mass_feature.per_charge_intensity = per_charge_intensity;
      mass_feature.per_isotope_intensity = per_isotope_intensity;

      auto apex = mt[mt.findMaxByIntPeak()];
      auto& sub_pg_map = peak_group_map_[apex.getRT()];
      auto& rep_pg = sub_pg_map[apex.getMZ()];
      mass_feature.rep_mz = apex.getMZ();
      mass_feature.scan_number = rep_pg.getScanNumber();
      mass_feature.rep_charge = rep_pg.getRepAbsCharge();
      mass_features.push_back(mass_feature);
    }
    return mass_features;
 }

  void MassFeatureTrace::storeInformationFromDeconvolvedSpectrum(DeconvolvedSpectrum& deconvolved_spectrum)
  {
    double rt = deconvolved_spectrum.getOriginalSpectrum().getRT();
    if (deconvolved_spectrum.getOriginalSpectrum().getMSLevel() != 1)
    {
      return;
    }
    else
    {
      peak_group_map_[rt] = std::map<double, PeakGroup>();
      auto& sub_pg_map = peak_group_map_[rt];
      for (auto& pg: deconvolved_spectrum)
      {
        sub_pg_map[pg.getMonoMass()] = pg;
      }
    }
  }

  void MassFeatureTrace::updateMembers_()
  {
    min_isotope_cosine_ = param_.getValue("min_isotope_cosine");
  }
}
