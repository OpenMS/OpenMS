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
    //mtd_defaults.setValue("min_trace_length", 10.0);
    mtd_defaults.setValue("min_sample_rate",
                          .1,
                          "Minimum fraction of scans along the feature trace that must contain a peak. To raise feature detection sensitivity, lower this value close to 0.");
    //mtd_defaults.setValue("mass_error_da",
    //                      1.5,
    //                      "da tolerance for feature tracing. Due to frequent isotope error, 1.5 Da is recommended.");
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

  void MassFeatureTrace::findFeatures(const String& file_name, const bool promex_out, const bool topfd_feature_out,
                                      const std::unordered_map<int, PeakGroup>& precursor_peak_groups,
                                      const PrecalculatedAveragine& averagine,
                                      int& feature_cntr,
                                      int& feature_index,
                                      std::fstream& fsf,
                                      std::fstream& fsp,
                                      std::vector<std::fstream>& fst)
  {
    MSExperiment map;
    std::map<int, MSSpectrum> index_spec_map;
    int min_abs_charge = INT_MAX;
    int max_abs_charge = INT_MIN;
    bool is_positive = true;
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

    for (auto& mt: m_traces)
    {
      double max_qscore = .0;
      int min_feature_abs_charge = INT_MAX; // min feature charge
      int max_feature_abs_charge = INT_MIN; // max feature charge

      //int minFIso = INT_MAX; // min feature isotope index
      //int maxFIso = INT_MIN; // max feature isotope index

      auto per_charge_intensity = std::vector<double>(charge_range + 1, .0);
      auto per_charge_max_intensity = std::vector<double>(charge_range + 1, .0);
      auto per_isotope_intensity = std::vector<double>(averagine.getMaxIsotopeIndex(), .0);

      int min_scan_num = -1;
      int max_scan_num = 0;
      double max_intensity = 0;
      double max_iso = 0;
      boost::dynamic_bitset<> charges(charge_range + 1);

      for (auto& p2: mt)
      {
        auto& pg_map = peak_group_map_[p2.getRT()];
        auto& pg = pg_map[p2.getMZ()];
        int scan_number = pg.getScanNumber();
        auto crange = pg.getAbsChargeRange();

        min_feature_abs_charge =
            min_feature_abs_charge < std::get<0>(crange) ? min_feature_abs_charge : std::get<0>(crange);
        max_feature_abs_charge =
            max_feature_abs_charge > std::get<1>(crange) ? max_feature_abs_charge : std::get<1>(crange);

        if (pg.getIsotopeCosine() > max_iso)
        {
          max_iso = pg.getIsotopeCosine();
        }

        for (auto& p: pg)
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
        }

        max_qscore = max_qscore < pg.getQScore() ? pg.getQScore() : max_qscore;
      }

      int offset = 0;

      double mass = mt.getCentroidMZ();
      double isotope_score = FLASHDeconvAlgorithm::getIsotopeCosineAndDetermineIsotopeIndex(mass,
                                                                                            per_isotope_intensity,
                                                                                            offset, averagine);
      if (isotope_score < min_isotope_cosine_)
      {
        continue;
      }

      double sum_intensity = .0;

      for (auto& p: mt)
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

    }


    if (topfd_feature_out || promex_out)
    {
      std::set<int> selected_traces_indices;

      int topid = 1;
      int promexid = 1;

      auto per_isotope_intensity = std::vector<double>(averagine.getMaxIsotopeIndex(), .0);

      for (auto& precursor: precursor_peak_groups)
      {
        int ms1_scan_number = precursor.second.getScanNumber();
        int ms2_scan_number = precursor.first;
        auto rt = scan_rt_map[ms2_scan_number];
        bool selected = false;
        MassTrace smt;
        double max_intensity = .0;
        int selected_index = -1;
        for (int l = 0; l < m_traces.size(); l++)
        {
          auto mt = m_traces[l];
          if (abs(precursor.second.getMonoMass() - mt.getCentroidMZ()) > 1.5)
          {
            continue;
          }
          if (rt < mt.begin()->getRT() || rt > mt.rbegin()->getRT())
          {
            continue;
          }

          double sum_intensity = .0;

          for (auto& p: mt)
          {
            sum_intensity += p.getIntensity();
          }

          if (max_intensity < sum_intensity)
          {
            max_intensity = sum_intensity;
            smt = mt;
            selected_index = l;
          }
          selected = true;
        }

        if (selected&&  selected_traces_indices.find(selected_index) == selected_traces_indices.end())
        {
          selected_traces_indices.insert(selected_index);
          double sum_intensity = .0;
          int min_feature_abs_charge = INT_MAX;
          int max_feature_abs_charge = INT_MIN;
          int min_scan_num = -1;
          int max_scan_num = 0;
          int rep_scan = 0, rep_charge = 0;
          max_intensity = .0;
          for (auto& p: smt)
          {
            sum_intensity += p.getIntensity();
            auto& pg_map = peak_group_map_[p.getRT()];
            auto& pg = pg_map[p.getMZ()];
            auto crange = pg.getAbsChargeRange();
            int scan_number = pg.getScanNumber();

            min_feature_abs_charge =
                min_feature_abs_charge < std::get<0>(crange) ? min_feature_abs_charge : std::get<0>(
                    crange);
            max_feature_abs_charge =
                max_feature_abs_charge > std::get<1>(crange) ? max_feature_abs_charge : std::get<1>(
                    crange);

            for (auto& pp: pg)
            {
              if (pp.isotopeIndex < 0 || pp.isotopeIndex >= averagine.getMaxIsotopeIndex() ||
                  pp.abs_charge < min_abs_charge ||
                  pp.abs_charge >= charge_range + min_abs_charge + 1)
              {
                continue;
              }

              per_isotope_intensity[pp.isotopeIndex] += pp.intensity;
            }

            if (min_scan_num < 0)
            {
              min_scan_num = scan_number;
            }
            else
            {
              min_scan_num = min_scan_num < scan_number ? min_scan_num : scan_number;
            }
            max_scan_num = max_scan_num > scan_number ? max_scan_num : scan_number;

            if (pg.getIntensity() > max_intensity)
            {
              max_intensity = pg.getIntensity();
              rep_scan = scan_number;
              rep_charge = pg.getRepAbsCharge();
            }
          }

          for (int i = 0; i < fst.size(); i++)
          {
            if (topfd_feature_out)
            {
              if (i == 0)
              {
                fst[i] << "0\t" << topid << "\t" << smt.getCentroidMZ() << "\t" << sum_intensity << "\t"
                       << smt.begin()->getRT() << "\t" << smt.rbegin()->getRT()
                       << "\t" << (is_positive ? min_feature_abs_charge : -max_feature_abs_charge) << "\t"
                       << (is_positive ? max_feature_abs_charge : -min_feature_abs_charge) << "\t1\t1\n";
              }
              else
              {
                fst[i] << ms2_scan_number << "\t0\t" << file_name << "\t" << ms2_scan_number << "\t" << ms1_scan_number
                       << "\t" << ms1_scan_number << "\t" << precursor.second.getMonoMass() << "\t"
                       << precursor.second.getIntensity() << "\t" << topid << "\t" << sum_intensity << "\t-1000\t"
                       << topid << "\t" << sum_intensity << "\n";
              }
              topid++;
            }
            if (promex_out&&  i == 0)
            {
              auto apex = smt[smt.findMaxByIntPeak()];
              auto& sub_pg_map = peak_group_map_[apex.getRT()];
              auto& spg = sub_pg_map[apex.getMZ()];
              auto mz_range = spg.getMzRange(rep_charge);
              double rmz = (std::get<0>(mz_range) + std::get<1>(mz_range)) / 2.0;

              fsp << promexid << "\t" << min_scan_num << "\t" << max_scan_num << "\t"
                  << (is_positive ? min_feature_abs_charge : -max_feature_abs_charge) << "\t"
                  << (is_positive ? max_feature_abs_charge : -min_feature_abs_charge) << "\t"
                  << std::to_string(smt.getCentroidMZ())
                  << "\t"
                  << std::fixed << std::setprecision(2)
                  << rep_scan << "\t" << (is_positive ? rep_charge : -rep_charge) << "\t" << rmz
                  << "\t"
                  << sum_intensity << "\t"
                  << spg.getScanNumber() << "\t" << spg.getIntensity() << "\t"
                  << smt.begin()->getRT() / 60.0 << "\t"
                  << smt.rbegin()->getRT() / 60.0 << "\t"
                  << smt.getTraceLength() / 60.0 << "\t";

              int offset = 0;
              double isotope_score = FLASHDeconvAlgorithm::getIsotopeCosineAndDetermineIsotopeIndex(smt.getCentroidMZ(),
                                                                                                    per_isotope_intensity,
                                                                                                    offset,
                                                                                                    averagine);

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

              promexid++;
            }
          }
          continue;
        }
        auto crange = precursor.second.getAbsChargeRange();
        for (int i = 0; i < fst.size(); i++)
        {
          if (topfd_feature_out)
          {
            if (i == 0)
            {
              fst[i] << "0\t" << topid << "\t" << precursor.second.getMonoMass() << "\t"
                     << precursor.second.getIntensity() << "\t" << rt - 1 << "\t" << rt + 1
                     << "\t" << (is_positive ? std::get<0>(crange) : -std::get<1>(crange)) << "\t"
                     << (is_positive ? std::get<1>(crange) : -std::get<0>(crange)) << "\t1\t1\n";
            }
            else
            {
              fst[i] << ms2_scan_number << "\t0\t" << file_name << "\t" << ms2_scan_number << "\t" << ms1_scan_number
                     << "\t" << ms1_scan_number << "\t" << precursor.second.getMonoMass() << "\t"
                     << precursor.second.getIntensity() << "\t" << topid << "\t" << precursor.second.getIntensity()
                     << "\t-1000\t"
                     << topid << "\t" << precursor.second.getIntensity() << "\n";


            }
            topid++;
          }
          if (promex_out&&  i == 0)
          {
            auto mz_range = precursor.second.getMaxQScoreMzRange();
            double rmz = (std::get<0>(mz_range) + std::get<1>(mz_range)) / 2.0;

            for (auto& pp: precursor.second)
            {
              if (pp.isotopeIndex < 0 || pp.isotopeIndex >= averagine.getMaxIsotopeIndex() ||
                  pp.abs_charge < min_abs_charge ||
                  pp.abs_charge >= charge_range + min_abs_charge + 1)
              {
                continue;
              }

              per_isotope_intensity[pp.isotopeIndex] += pp.intensity;
            }

            fsp << promexid << "\t" << precursor.second.getScanNumber() << "\t" << precursor.second.getScanNumber()
                << "\t"
                << (is_positive ? std::get<0>(crange) : -std::get<1>(crange)) << "\t"
                << (is_positive ? std::get<1>(crange) : -std::get<0>(crange)) << "\t"
                << std::to_string(precursor.second.getMonoMass())
                << "\t"
                << std::fixed << std::setprecision(2)
                << precursor.second.getScanNumber() << "\t"
                << (is_positive ? precursor.second.getRepAbsCharge() : -precursor.second.getRepAbsCharge())
                << "\t" << rmz << "\t"
                << precursor.second.getIntensity() << "\t"
                << precursor.second.getScanNumber() << "\t" << precursor.second.getIntensity() << "\t"
                << (rt - 1) / 60.0 << "\t"
                << (rt + 1) / 60.0 << "\t"
                << 2.0 / 60.0 << "\t";

            int offset = 0;
            double isotope_score = precursor.second.getIsotopeCosine();

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
            promexid++;
          }
        }
      }


    }

  }

  void MassFeatureTrace::storeInformationFromDeconvolvedSpectrum(DeconvolvedSpectrum& deconvolved_spectrum)
  {
    double rt = deconvolved_spectrum.getOriginalSpectrum().getRT();
    scan_rt_map[deconvolved_spectrum.getScanNumber()] = rt;
    if (deconvolved_spectrum.getOriginalSpectrum().getMSLevel() != 1)
    {
      return;
    }
    else
    {

      peak_group_map_[rt] = std::unordered_map<double, PeakGroup>();
      auto& sub_pg_map = peak_group_map_[rt];
      for (auto& pg: deconvolved_spectrum)
      {
        sub_pg_map[pg.getMonoMass()] = pg;
      }
      //scan_rt_map[deconvolved_spectrum.getScanNumber()] = rt;
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

  void MassFeatureTrace::writeTopFDFeatureHeader(std::vector<std::fstream>& fs)
  {
    for (int i = 0; i < fs.size(); ++i)
    {
      if (i == 0)
      {
        fs[i]
            << "Sample_ID\tID\tMass\tIntensity\tTime_begin\tTime_end\tMinimum_charge_state\tMaximum_charge_state\tMinimum_fraction_id\tMaximum_fraction_id\n";
      }
      else
      {
        fs[i]
            << "Spec_ID\tFraction_ID\tFile_name\tScans\tMS_one_ID\tMS_one_scans\tPrecursor_mass\tPrecursor_intensity\tFraction_feature_ID\tFraction_feature_intensity\tFraction_feature_score\tSample_feature_ID\tSample_feature_intensity\n";
      }
    }


  }

  void MassFeatureTrace::updateMembers_()
  {
    //tol_ = param_.getValue("mass_error_ppm");
    //minChargeCosine = param_.getValue("min_charge_cosine");
    min_isotope_cosine_ = param_.getValue("min_isotope_cosine");
  }
}
