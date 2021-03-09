//--------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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

#include "OpenMS/ANALYSIS/TOPDOWN/MassFeatureTrace.h"
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <utility>

namespace OpenMS
{
  MassFeatureTrace::MassFeatureTrace() :
      DefaultParamHandler("MassFeatureTrace")
  {
    Param mtd_defaults = MassTraceDetection().getDefaults();
    mtd_defaults.setValue("min_trace_length", 10.0);
    mtd_defaults.setValue("min_sample_rate",
                          .2,
                          "Minimum fraction of scans along the feature trace that must contain a peak. To raise feature detection sensitivity, lower this value close to 0.");

    mtd_defaults.setValue("mass_error_da",
                          1.5,
                          "da tolerance for feature tracing. Due to frequent isotope errer, 1.5 Da is recommended.");
    mtd_defaults.setValue("min_trace_length", 10.0);//

    mtd_defaults.setValue("chrom_peak_snr", .0);
    mtd_defaults.addTag("chrom_peak_snr", "advanced");
    mtd_defaults.setValue("reestimate_mt_sd",
                          "false");
    mtd_defaults.addTag("reestimate_mt_sd", "advanced");
    mtd_defaults.setValue("noise_threshold_int",
                          .0);
    mtd_defaults.addTag("noise_threshold_int", "advanced");
    //mtd_defaults.setValue("min_isotope_cosine", -1.0, "if not set, controlled by -Algorithm:min_isotope_cosine_ option");
    //mtd_defaults.addTag("min_isotope_cosine", "advanced");

    mtd_defaults.setValue("quant_method", "area");
    mtd_defaults.addTag("quant_method", "advanced"); // hide entry

    defaults_.insert("", mtd_defaults);
    defaults_.setValue("min_isotope_cosine", .75, "cosine threshold between avg. and observed isotope pattern for MS1");
    defaultsToParam_();
  }

  MassFeatureTrace::~MassFeatureTrace()
  {
    for (auto& item : peak_group_map_)
    {
      std::unordered_map<double, PeakGroup>().swap(item.second);
    }
    std::unordered_map<double, std::unordered_map<double, PeakGroup>>().swap(peak_group_map_);
  }


  void MassFeatureTrace::findFeatures(const String& file_name, const bool promex_out,
                                      int& feature_cntr,
                                      int& feature_index,
                                      std::fstream& fsf,
                                      std::fstream& fsp,
                                      const PrecalculatedAveragine& averagine)
  {
    MSExperiment map;
    std::map<int, MSSpectrum> index_spec_map;
    int min_abs_charge = INT_MAX;
    int max_abs_charge = INT_MIN;
    bool is_positive = true;
    for (auto &item : peak_group_map_)
    {
      double rt = item.first;
      MSSpectrum deconv_spec;
      deconv_spec.setRT(rt);
      for (auto &pg : item.second)
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

    if (map.size() < 3)
    {
      return;
    }

    map.sortSpectra();
    MassTraceDetection mtdet;
    Param mtd_param = getParameters().copy("");
    mtd_param.remove("min_isotope_cosine");

    mtdet.setParameters(mtd_param);
    std::vector<MassTrace> m_traces;

    mtdet.run(map, m_traces);  // m_traces : output of this function

    int charge_range = max_abs_charge - min_abs_charge + 1;

    for (auto& mt : m_traces)
    {
      double max_qscore = .0;
      int min_feature_abs_charge = INT_MAX; // min feature charge
      int max_feature_abs_charge = INT_MIN; // max feature charge

      //int minFIso = INT_MAX; // min feature isotope index
      //int maxFIso = INT_MIN; // max feature isotope index

      auto per_charge_intensity = std::vector<double>(charge_range + 1, 0);
      auto per_charge_max_intensity = std::vector<double>(charge_range + 1, 0);
      auto per_charge_mz = std::vector<double>(charge_range + 1, 0);
      auto per_isotope_intensity = std::vector<double>(averagine.getMaxIsotopeIndex(), 0);

      int min_scan_num = (int) map.size() + 1000;
      int max_scan_num = 0;

      int rep_scan = 0, rep_charge = 0;
      double max_intensity = 0;
      double max_mass = .0;
      double max_iso = 0;
      boost::dynamic_bitset<> charges(charge_range + 1);

      for (auto& p2 : mt)
      {
        auto &pg_map = peak_group_map_[p2.getRT()];
        auto &pg = pg_map[p2.getMZ()];
        int scan_number = pg.getScanNumber();
        auto crange = pg.getAbsChargeRange();

        min_feature_abs_charge =
            min_feature_abs_charge < std::get<0>(crange) ? min_feature_abs_charge : std::get<0>(crange);
        max_feature_abs_charge =
            max_feature_abs_charge > std::get<1>(crange) ? max_feature_abs_charge : std::get<1>(crange);

        min_scan_num = min_scan_num < scan_number ? min_scan_num : scan_number;
        max_scan_num = max_scan_num > scan_number ? max_scan_num : scan_number;

        if (pg.getIntensity() > max_intensity)
        {
          max_intensity = pg.getIntensity();
          rep_scan = scan_number;

        }

        if (pg.getIsotopeCosine() > max_iso)
        {
          max_iso = pg.getIsotopeCosine();
          max_mass = pg.getMonoMass();
        }

        for (auto& p : pg)
        {
          if (p.isotopeIndex < 0 || p.isotopeIndex >= averagine.getMaxIsotopeIndex() || p.abs_charge < min_abs_charge ||
              p.abs_charge >= charge_range + min_abs_charge + 1)
          {
            continue;
          }

          charges[p.abs_charge - min_abs_charge] = true;
          per_charge_intensity[p.abs_charge - min_abs_charge] += p.intensity;
          per_isotope_intensity[p.isotopeIndex] += p.intensity;
          if (per_charge_max_intensity[p.abs_charge - min_abs_charge] > p.intensity)
          {
            continue;
          }
          per_charge_max_intensity[p.abs_charge - min_abs_charge] = p.intensity;
          per_charge_mz[p.abs_charge - min_abs_charge] = p.mz;
        }

        max_qscore = max_qscore < pg.getQScore() ? pg.getQScore() : max_qscore;
      }

      //double charge_score_ = FLASHDeconvAlgorithm::getChargeFitScore(per_charge_intensity,
      //                                                            charge_range);
      //if (charge_score_ < minChargeCosine) //
      //{
      //  continue;
      //}

      int offset = 0;

      double mass = mt.getCentroidMZ();
      double isotope_score = FLASHDeconvAlgorithm::getIsotopeCosineAndDetermineIsotopeIndex(mass,
                                                                                            per_isotope_intensity,
                                                                                            offset, averagine);
      if (isotope_score < min_isotope_cosine_)
      {
        continue;
      }
      if (offset != 0)
      {
        mass += offset * Constants::ISOTOPE_MASSDIFF_55K_U;
      }

      double sum_intensity = .0;

      for (auto& p : mt)
      {
        sum_intensity += p.getIntensity();
      }
      double avg_mass = averagine.getAverageMassDelta(mass) + mass;
      ++feature_cntr;
      fsf << feature_index++ << "\t" << file_name << "\t" << std::to_string(mass) << "\t"
          << std::to_string(avg_mass) << "\t" // massdiff
          << mt.getSize() << "\t"
          << mt.begin()->getRT() << "\t"
          << mt.rbegin()->getRT() << "\t"
          << mt.getTraceLength() << "\t"
          << mt[mt.findMaxByIntPeak()].getRT() << "\t"
          << sum_intensity << "\t"
          << mt.getMaxIntensity(false) << "\t"
          << mt.computePeakArea() << "\t"
          << (is_positive ? min_feature_abs_charge : -max_feature_abs_charge) << "\t"
          << (is_positive ? max_feature_abs_charge : -min_feature_abs_charge) << "\t"
          << charges.count() << "\t"
          << isotope_score << "\t"
          << max_qscore << "\t";

      if (is_positive)
      {
        for (int i = min_feature_abs_charge; i <= max_feature_abs_charge; i++)
        {
          fsf << per_charge_intensity[i - min_abs_charge];
          if (i < max_feature_abs_charge)
          {
            fsf << ";";
          }
        }
      }
      else
      {
        for (int i = max_feature_abs_charge; i >= min_feature_abs_charge; i--)
        {
          fsf << per_charge_intensity[i - min_abs_charge];
          if (i > min_feature_abs_charge)
          {
            fsf << ";";
          }
        }
      }
      fsf << "\t";
      int iso_end_index = 0;

      for (int i = 0; i < averagine.getMaxIsotopeIndex(); i++)
      {
        if (per_isotope_intensity[i] == 0)
        {
          continue;
        }
        iso_end_index = i;
      }
      for (int i = 0; i <= iso_end_index; i++)
      {
        fsf << per_isotope_intensity[i];
        if (i < iso_end_index)
        {
          fsf << ";";
        }
      }

      fsf << "\n";

      if (promex_out)
      {
        double max_charge_intensity = 0;
        for (int c = 0; c < charge_range; c++) // c is charge range!!
        {
          if (per_charge_intensity[c] > max_charge_intensity)
          {
            max_charge_intensity = per_charge_intensity[c];
            rep_charge = c + min_abs_charge;
          }
        }
        auto apex = mt[mt.findMaxByIntPeak()];

        //int si = rtSpecMap[(float) apex.getRT()];
        auto &sub_pg_map = peak_group_map_[apex.getRT()];
        auto &spg = sub_pg_map[apex.getMZ()];

        fsp << feature_index << "\t" << min_scan_num << "\t" << max_scan_num << "\t"
            << (is_positive ? min_feature_abs_charge : -max_feature_abs_charge) << "\t"
            << (is_positive ? max_feature_abs_charge : -min_feature_abs_charge) << "\t" << std::to_string(mass) << "\t"
            << std::fixed << std::setprecision(2)
            << rep_scan << "\t" << (is_positive ? rep_charge : -rep_charge) << "\t" << per_charge_mz[rep_charge] << "\t"
            << sum_intensity << "\t"
            << spg.getScanNumber() << "\t" << spg.getIntensity() << "\t"
            << mt.begin()->getRT() / 60.0 << "\t"
            << mt.rbegin()->getRT() / 60.0 << "\t"
            << mt.getTraceLength() / 60.0 << "\t";


        for (int j = 0; j < averagine.getMaxIsotopeIndex(); ++j)
        {
          if (per_isotope_intensity[j] <= 0)
          {
            continue;
          }
          fsp << j << "," << per_isotope_intensity[j] << ";";
        }
        fsp << "\t" << isotope_score << "\n";
        fsp << std::setprecision(0);
      }
    }
  }

  void MassFeatureTrace::storeInformationFromDeconvolutedSpectrum(DeconvolutedSpectrum& deconvoluted_spectrum)
  {
    if (deconvoluted_spectrum.getOriginalSpectrum().getMSLevel() != 1)
    {
      return;
    }
    double rt = deconvoluted_spectrum.getOriginalSpectrum().getRT();
    peak_group_map_[rt] = std::unordered_map<double, PeakGroup>();
    auto& sub_pg_map = peak_group_map_[rt];
    for (auto& pg : deconvoluted_spectrum)
    {
      sub_pg_map[pg.getMonoMass()] = pg;
    }
  }

  void MassFeatureTrace::writeHeader(std::fstream& fs)
  {
    fs << "FeatureIndex\tFileName\tMonoisotopicMass\tAverageMass\tMassCount\tStartRetentionTime"
          "\tEndRetentionTime\tRetentionTimeDuration\tApexRetentionTime"
          "\tSumIntensity\tMaxIntensity\tFeatureQuantity\tMinCharge\tMaxCharge\tChargeCount\tIsotopeCosineScore\tMaxQScore\tPerChargeIntensity\tPerIsotopeIntensity"
          "\n";
  }


  void MassFeatureTrace::writePromexHeader(std::fstream& fs)
  {
    fs << "FeatureID\tMinScan\tMaxScan\tMinCharge\tMaxCharge\t"
          "MonoMass\tRepScan\tRepCharge\tRepMz\tAbundance\tApexScanNum\tApexIntensity\tMinElutionTime\tMaxElutionTime\t"
          "ElutionLength\tEnvelope\tLikelihoodRatio"
          "\n";
  }

  void MassFeatureTrace::updateMembers_()
  {
    //tol_ = param_.getValue("mass_error_ppm");
    //minChargeCosine = param_.getValue("min_charge_cosine");
    min_isotope_cosine_ = param_.getValue("min_isotope_cosine");
  }
}
