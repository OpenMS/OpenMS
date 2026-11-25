// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong$
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FLASHDeconvSpectrumFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <algorithm>
#include <atomic>
#include <exception>
#include <memory>
#include <mutex>
#include <sstream>
#include <thread>
#ifdef _OPENMP
  #include <omp.h>
#endif

namespace OpenMS
{
  /**
   * @brief FLASHDeconv Spectrum level output *.tsv, *.msalign (for TopPIC) file formats
     @ingroup FileIO

   */


  namespace
  {
    // The number of target masses of a deconvolved spectrum. With report_FDR, a spectrum also holds the decoy masses; the counts and
    // limits of the target output are taken over the target masses only, so that it does not change with report_FDR.
    Size countTargetMasses(const DeconvolvedSpectrum& dspec)
    {
      return (Size)std::count_if(dspec.begin(), dspec.end(),
                                 [](const PeakGroup& pg) { return pg.getTargetDecoyType() == PeakGroup::TargetDecoyType::target; });
    }

    // The masses of the isolation window of an MSn spectrum (DeconvolvedSpectrum::getPrecursorMassIntensityMap) as two space separated
    // lists, the masses and their shares of the total window intensity, the precursor first; pg_intensity is the window intensity of the
    // precursor. The window masses are those of the one survey spectrum the precursor was found in. The precursor mass is the mass of
    // its MS1 feature at the isotope label chosen with the fragments, so it differs from its window mass by a few ppm (1-99 % of 222,663
    // precursors of 15 iodoTMT runs: -5.8..+5.5 ppm) and at most one isotope. The precursor is therefore the window mass nearest to the
    // precursor mass within 10 ppm, the same isotope before a neighbouring one. A precursor that is not among the window masses (taken
    // from another survey spectrum or given from outside) leaves pg_intensity 0 and the lists without a leading precursor entry.
    void windowMasses(const DeconvolvedSpectrum& dspec, double total_int, double& pg_intensity, std::string& mass_str, std::string& int_str)
    {
      const auto& window = dspec.getPrecursorMassIntensityMap();
      const double pm = dspec.getPrecursorPeakGroup().getMonoMass();
      const double iso = dspec.getPrecursorPeakGroup().getIsotopeDaDistance();
      int best = -1;
      double best_cost = 2.0; // the isotope distance (0 or 1) + the residual relative to 10 ppm (below 1)
      for (Size i = 0; i < window.size(); i++)
      {
        const double m = std::get<0>(window[i]);
        if (m <= 0) continue; // the total intensity entry
        const double k = std::round((m - pm) / iso);
        if (std::abs(k) > 1) continue;
        const double residual = std::abs(m - pm - k * iso) / (pm * 1e-5);
        if (residual >= 1) continue;
        const double cost = std::abs(k) + residual;
        if (cost < best_cost)
        {
          best_cost = cost;
          best = (int)i;
        }
      }
      std::stringstream mass_ss, int_ss;
      if (best >= 0)
      {
        pg_intensity = std::get<1>(window[best]);
        mass_ss << std::to_string(std::get<0>(window[best])) << " ";
        int_ss << pg_intensity / total_int << " ";
      }
      for (Size i = 0; i < window.size(); i++)
      {
        if ((int)i == best || std::get<0>(window[i]) <= 0) continue;
        mass_ss << std::to_string(std::get<0>(window[i])) << " ";
        int_ss << std::get<1>(window[i]) / total_int << " ";
      }
      mass_str = mass_ss.str();
      int_str = int_ss.str();
      if (! mass_str.empty()) mass_str.pop_back();
      if (! int_str.empty()) int_str.pop_back();
    }

    /// TopPIC's name of an activation method: TopPIC knows CID, HCD, ETD, ETCID, ETHCD, MPD and UVPD (res/base_data/activation_base.xml
    /// of TopPIC 2.1) and stops at any other name; ECD gives the c/z-dot ions of ETD, photodissociation those of UVPD, infrared
    /// multiphoton dissociation those of MPD, the other collisional methods those of CID. Empty if unknown.
    std::string topPICActivation_(const Precursor::ActivationMethod method)
    {
      switch (method)
      {
        case Precursor::ActivationMethod::HCD: return "HCD";
        case Precursor::ActivationMethod::ETD:
        case Precursor::ActivationMethod::ECD: return "ETD";
        case Precursor::ActivationMethod::ETciD: return "ETCID";
        case Precursor::ActivationMethod::EThcD: return "ETHCD";
        case Precursor::ActivationMethod::PHD: return "UVPD";
        case Precursor::ActivationMethod::IMD: return "MPD";
        case Precursor::ActivationMethod::CID:
        case Precursor::ActivationMethod::HCID:
        case Precursor::ActivationMethod::LCID:
        case Precursor::ActivationMethod::TRAP:
        case Precursor::ActivationMethod::PQD:
        case Precursor::ActivationMethod::SID:
        case Precursor::ActivationMethod::SORI:
        case Precursor::ActivationMethod::BIRD:
        case Precursor::ActivationMethod::PSD:
        case Precursor::ActivationMethod::INSOURCE: return "CID";
        default: return "";
      }
    }

    /// the isolation window target m/z of @p precursor (mzML), its m/z if not given
    double isolationTargetMz_(const Precursor& precursor)
    {
      return precursor.metaValueExists("isolation window target m/z") ? (double)precursor.getMetaValue("isolation window target m/z") : precursor.getMZ();
    }

    /// the row index of each MS level in the spectrum tsv files (1-based), continued over all calls of the writers
    std::vector<uint>& writerIndices_()
    {
      static std::vector<uint> indices {};
      return indices;
    }

    /// With report_decoy, the noise decoy masses are written noise_decoy_weight times each on average, without randomness: each adds the weight
    /// to @p carry, and is written as many times as whole units are in it (taken out). Returns how often peak group @p pg is written.
    Size noiseDecoyCopies_(const PeakGroup& pg, const double noise_decoy_weight, double& carry)
    {
      if (pg.getTargetDecoyType() != PeakGroup::TargetDecoyType::noise_decoy) return 1;
      carry += noise_decoy_weight;
      const Size copies = (Size)std::floor(carry);
      carry -= (double)copies;
      return copies;
    }

    /// the number of rows writeDeconvolvedMasses_ writes for @p dspec with report_decoy, advancing @p carry exactly as it does
    Size countDecoyReportRows_(const DeconvolvedSpectrum& dspec, const double noise_decoy_weight, double& carry)
    {
      Size rows = 0;
      for (const auto& pg : dspec) rows += noiseDecoyCopies_(pg, noise_decoy_weight, carry);
      return rows;
    }

    /// the noise decoy carry of each MS level (by MS level - 1) of the single-spectrum writer, continued over its calls (as writerIndices_)
    std::vector<double>& writerCarries_()
    {
      static std::vector<double> carries {};
      return carries;
    }

    /// the weight of MS level @p ms_level in @p weights (by MS level - 1), 1 if absent
    double levelWeight_(const std::vector<double>& weights, const uint ms_level)
    {
      return ms_level >= 1 && ms_level <= weights.size() ? weights[ms_level - 1] : 1.0;
    }

    /// the rows of the non-empty spectrum @p dspec into @p fs, numbered from @p index on (advanced past them); noise decoys are written as
    /// noiseDecoyCopies_ says, from @p carry on (advanced)
    void writeDeconvolvedMasses_(DeconvolvedSpectrum& dspec, std::ostream& fs, const std::string& file_name, const FLASHHelperClasses::PrecalculatedAveragine& avg,
                                 const FLASHHelperClasses::PrecalculatedAveragine& decoy_avg, double tol, const bool write_detail, const bool report_decoy,
                                 const double noise_decoy_weight, uint& index, double& carry,
                                 const std::vector<std::pair<double, Size>>* abundances = nullptr)
    {
    std::stringstream ss;

    std::stringstream precursor_ss;
    if (dspec.getOriginalSpectrum().getMSLevel() > 1)
    {
      // the decoy precursor mass of the isolation window (FLASHTnT -precursor_FDR) and its weight in the FDR
      const std::string decoy_precursor = dspec.getDecoyPrecursorMass() > 0
                                       ? std::to_string(dspec.getDecoyPrecursorMass()) + "\t" + StringUtils::toStr(dspec.getDecoyPrecursorWeight()) + "\t"
                                       : std::string("NaN\tNaN\t");
      precursor_ss << dspec.getPrecursorScanNumber() << "\t"
                   << (dspec.getPrecursorPeakGroup().getFeatureIndex() == 0 ? "NaN" : std::to_string(dspec.getPrecursorPeakGroup().getFeatureIndex()))
                   << "\t" << std::to_string(dspec.getPrecursor().getMZ()) << "\t"  << dspec.getPrecursor().getIntensity() << "\t";

      double total_int = 0;
      for (const auto& [m, it] : dspec.getPrecursorMassIntensityMap())
      {
        if (m <= 0)
        {
          total_int = it;
          break;
        }
      }

      if (dspec.getPrecursorPeakGroup().empty())
      {
        precursor_ss << "NaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\t0\tNaN\tNaN\t" << std::to_string(total_int) << "\t";
        if (report_decoy) precursor_ss << "NaN\t" << decoy_precursor;
      }
      else
      {
        double pg_intensity = 0;
        std::string mass_str, int_str;
        windowMasses(dspec, total_int, pg_intensity, mass_str, int_str);

        precursor_ss << dspec.getPrecursor().getCharge() << "\t" << dspec.getPrecursorPeakGroup().getChargeSNR(dspec.getPrecursor().getCharge()) << "\t"
                     << std::to_string(dspec.getPrecursorPeakGroup().getMonoMass()) << "\t"
                     << std::to_string(dspec.getPrecursorPeakGroup().getQscore()) << "\t" << dspec.getPrecursorPeakGroup().getQscore2D() << "\t";

        precursor_ss <<std::to_string(dspec.getDeconvolvedPrecursor().getMZ()) << "\t" << dspec.getDeconvolvedPrecursor().getIntensity()  << "\t" << pg_intensity << "\t"
                     << dspec.getPrecursorMassIntensityMap().size() - 1 << "\t" << mass_str << "\t" << int_str << "\t" << total_int << "\t";

        if (report_decoy) { precursor_ss << dspec.getPrecursorPeakGroup().getQvalue() << "\t" << decoy_precursor; }
      }
    }

    const Size target_mass_count = countTargetMasses(dspec);
    // log m/z of the original peaks, computed once for the noisy peaks of all masses of the spectrum
    std::vector<double> log_mzs;
    bool log_mzs_positive = true;
    std::vector<FLASHHelperClasses::LogMzPeak> noisy_peaks;
    Size remaining = 0; // the rows of the current mass still to be written
    for (int i = 0; i < dspec.size(); i++)
    {
      auto& pg = dspec[i];
      if (!report_decoy && pg.getTargetDecoyType() != PeakGroup::TargetDecoyType::target) continue;

      if (remaining == 0) remaining = report_decoy ? noiseDecoyCopies_(pg, noise_decoy_weight, carry) : 1;
      if (remaining == 0) continue;
      if (--remaining > 0) i--; // written again

      const auto& avg_ = pg.getTargetDecoyType() == PeakGroup::TargetDecoyType::noise_decoy? decoy_avg : avg;
      const double mono_mass = pg.getMonoMass();
      const double avg_mass = pg.getMonoMass() + avg_.getAverageMassDelta(mono_mass);
      const double intensity = pg.getIntensity();

      const auto& charge_range = pg.getAbsChargeRange();
      int min_charge = pg.isPositive() ? std::get<0>(charge_range) : -std::get<1>(charge_range);
      int max_charge = pg.isPositive() ? std::get<1>(charge_range) : -std::get<0>(charge_range);

      pg.setIndex(index);
      ss << index++ << "\t" << file_name << "\t" << pg.getScanNumber() << "\t" << (pg.getFeatureIndex() == 0 ? "nan" : std::to_string(pg.getFeatureIndex())) << "\t";

      if (report_decoy)
      {
        ss << pg.getTargetDecoyType() << "\t";
      }
      ss << std::to_string(dspec.getOriginalSpectrum().getRT()/60.0) << "\t" << target_mass_count << "\t" << std::to_string(avg_mass) << "\t" << std::to_string(mono_mass) << "\t" << intensity << "\t"
         << min_charge << "\t" << max_charge << "\t" << pg.size() << "\t";

      if (write_detail)
      {
        const auto& spec = dspec.getOriginalSpectrum();
        if (log_mzs.size() != spec.size() || log_mzs_positive != pg.isPositive())
        {
          log_mzs_positive = pg.isPositive();
          log_mzs.resize(spec.size());
          for (Size k = 0; k < spec.size(); k++) log_mzs[k] = FLASHHelperClasses::getLogMz(spec[k].getMZ(), log_mzs_positive);
        }
        pg.recruitAllPeaksInSpectrum(spec, tol * 1e-6, avg_, pg.getMonoMass(), noisy_peaks, false, &log_mzs);

        std::sort(noisy_peaks.begin(), noisy_peaks.end());
        ss << std::fixed << std::setprecision(2);
        for (auto& p : pg)
        {
          ss << std::to_string(p.mz) << " ";
        }

        ss << "\t";
        ss << std::fixed << std::setprecision(1);
        for (auto& p : pg)
        {
          ss << p.intensity << " ";
        }

        ss << "\t";
        ss << std::setprecision(-1);
        for (auto& p : pg)
        {
          ss << (p.is_positive ? p.abs_charge : -p.abs_charge) << " ";
        }

        ss << "\t";
        for (auto& p : pg)
        {
          ss << p.getUnchargedMass() << " ";
        }

        ss << "\t";
        for (auto& p : pg)
        {
          ss << p.isotopeIndex << " ";
        }

        ss << "\t";
        ss << std::setprecision(2);
        for (auto& p : pg)
        {
          double average_mass = pg.getMonoMass() + p.isotopeIndex * pg.getIsotopeDaDistance();
          double mass_error = (average_mass / p.abs_charge + FLASHHelperClasses::getChargeMass(p.is_positive) - p.mz) / p.mz;
          ss << 1e6 * mass_error << " ";
        }
        ss << std::setprecision(-1);
        ss << "\t";
        ss << std::fixed << std::setprecision(2);
        for (auto& np : noisy_peaks)
        {
          ss << std::to_string(np.mz) << " ";
        }

        ss << "\t";
        ss << std::fixed << std::setprecision(1);
        for (auto& np : noisy_peaks)
        {
          ss << np.intensity << " ";
        }

        ss << "\t";
        ss << std::setprecision(-1);
        for (auto& np : noisy_peaks)
        {
          ss << (np.is_positive ? np.abs_charge : -np.abs_charge) << " ";
        }

        ss << "\t";
        for (auto& np : noisy_peaks)
        {
          ss << np.getUnchargedMass() << " ";
        }

        ss << "\t";
        for (auto& np : noisy_peaks)
        {
          ss << np.isotopeIndex << " ";
        }

        ss << "\t";
        ss << std::setprecision(2);
        for (auto& np : noisy_peaks)
        {
          double average_mass = pg.getMonoMass() + np.isotopeIndex * pg.getIsotopeDaDistance();
          double mass_error = (average_mass / np.abs_charge + FLASHHelperClasses::getChargeMass(np.is_positive) - np.mz) / np.mz;
          ss << 1e6 * mass_error << " ";
        }
        ss << std::setprecision(-1);
        ss << "\t";
      }

      if (dspec.getOriginalSpectrum().getMSLevel() > 1)
      {
        ss << precursor_ss.str();
      }
      ss << pg.getIsotopeCosine() << "\t" << pg.getChargeIsotopeCosine(pg.getRepAbsCharge()) << "\t" << pg.getChargeScore() << "\t";

      const auto& max_qscore_mz_range = pg.getRepMzRange();
      ss << pg.getSNR() << "\t" << pg.getChargeSNR(pg.getRepAbsCharge()) << "\t" << pg.getAvgPPMError() << "\t" << (pg.isPositive() ? pg.getRepAbsCharge() : -pg.getRepAbsCharge()) << "\t"
         << std::to_string(std::get<0>(max_qscore_mz_range)) << "\t" << std::to_string(std::get<1>(max_qscore_mz_range)) << "\t" << std::to_string(pg.getQscore()) << "\t"
         << std::to_string(pg.getQscore2D());

      if (report_decoy)
      {
        ss << "\t" << pg.getQvalue();
      }

      if (write_detail)
      {
        ss << "\t" << std::setprecision(-1);

        for (int j = min_charge; j <= max_charge; j++)
        {
          ss << pg.getChargeIntensity(j);

          if (j < max_charge)
          {
            ss << ";";
          }
        }
        ss << "\t";

        const auto& iso_intensities = pg.getIsotopeIntensities();
        for (size_t j = 0; j < iso_intensities.size(); j++)
        {
          ss << iso_intensities[j];
          if (j < iso_intensities.size() - 1)
          {
            ss << ";";
          }
        }
      }
      if (abundances != nullptr && dspec.getOriginalSpectrum().getMSLevel() == 1)
      {
        if (pg.getTargetDecoyType() == PeakGroup::TargetDecoyType::target && i >= 0 && (Size)i < abundances->size())
        {
          ss << "\t" << std::to_string((*abundances)[i].first) << "\t" << (*abundances)[i].second;
        }
        else { ss << "\tnan\t0"; }
      }
      ss << "\n";
    }
    fs << ss.str();
  }
  } // namespace

  void FLASHDeconvSpectrumFile::writeDeconvolvedMasses(DeconvolvedSpectrum& dspec, std::fstream& fs, const std::string& file_name, const FLASHHelperClasses::PrecalculatedAveragine& avg, const FLASHHelperClasses::PrecalculatedAveragine& decoy_avg, double tol,
                                                       const bool write_detail, const bool report_decoy, const double noise_decoy_weight,
                                                       const std::vector<std::pair<double, Size>>* abundances)
  {
    if (!report_decoy && dspec.isDecoy()) return;
    if (dspec.empty()) return;
    auto& indices = writerIndices_();
    auto& carries = writerCarries_();
    while (indices.size() <= dspec.getOriginalSpectrum().getMSLevel())
    {
      indices.push_back(1);
      carries.push_back(.5);
    }
    writeDeconvolvedMasses_(dspec, fs, file_name, avg, decoy_avg, tol, write_detail, report_decoy, noise_decoy_weight, indices[dspec.getOriginalSpectrum().getMSLevel() - 1],
                            carries[dspec.getOriginalSpectrum().getMSLevel() - 1], abundances);
  }

  void FLASHDeconvSpectrumFile::writeDeconvolvedMasses(std::vector<DeconvolvedSpectrum>& spectra, std::vector<std::fstream>& streams,
                                                       const std::vector<std::string>& out_files, const std::string& file_name,
                                                       const FLASHHelperClasses::PrecalculatedAveragine& avg,
                                                       const FLASHHelperClasses::PrecalculatedAveragine& decoy_avg, const std::vector<double>& tols,
                                                       const bool write_detail, const bool report_decoy, const std::vector<double>& noise_decoy_weights,
                                                       const std::vector<std::vector<std::pair<double, Size>>>* abundances)
  {
    // the spectra of the files to be written, in input order
    std::vector<Size> jobs;
    for (Size s = 0; s < spectra.size(); s++)
    {
      const uint ms_level = spectra[s].getOriginalSpectrum().getMSLevel();
      if (ms_level == 0 || ms_level > out_files.size() || out_files[ms_level - 1].empty()) continue;
      jobs.push_back(s);
    }
    // A spectrum's rows depend only on the spectrum, the index of its first row and, with report_decoy, the noise decoy carry of its MS level
    // (one carry per MS level in spectrum order). Both are advanced up front per spectrum, in input order (by one per target mass, or with
    // report_decoy by replaying the carry). The spectra are then formatted in parallel, and the formatted rows are written in
    // input order as they become ready, by whichever thread holds the write lock: no thread waits at a barrier, and a buffer is released as
    // soon as it is written.
    auto& indices = writerIndices_();
    const Size n = jobs.size();
    std::vector<uint> first_index(n, 0);
    std::vector<double> first_carry(n, .5); // the noise decoy carry each spectrum starts from (with report_decoy), per MS level from 0.5 on
    std::vector<double> carries; // by MS level - 1
    std::vector<char> skip(n, 0);
    for (Size k = 0; k < n; k++)
    {
      const auto& dspec = spectra[jobs[k]];
      if ((! report_decoy && dspec.isDecoy()) || dspec.empty())
      {
        skip[k] = 1;
        continue;
      }
      const uint ms_level = dspec.getOriginalSpectrum().getMSLevel();
      while (indices.size() <= ms_level) indices.push_back(1);
      first_index[k] = indices[ms_level - 1];
      while (carries.size() < ms_level) carries.push_back(.5);
      first_carry[k] = carries[ms_level - 1];
      indices[ms_level - 1] += (uint)(report_decoy ? countDecoyReportRows_(dspec, levelWeight_(noise_decoy_weights, ms_level), carries[ms_level - 1])
                                                   : countTargetMasses(dspec));
    }

    int threads = 1;
#ifdef _OPENMP
    threads = omp_get_max_threads();
#endif
    // a spectrum is formatted at most this many spectra ahead of the next one to be written: bounds the formatted text held at once
    const Size window = 16 * (Size)std::max(1, threads);
    std::vector<std::string> buffers(n);
    std::vector<std::exception_ptr> errors(n);
    std::unique_ptr<std::atomic<bool>[]> ready(new std::atomic<bool>[n]);
    for (Size k = 0; k < n; k++) ready[k].store(false, std::memory_order_relaxed);
    std::atomic<Size> next_write {0};  // the next spectrum to be written
    std::atomic<bool> failed {false};  // a spectrum failed (its exception is rethrown after the loop; nothing after it is written)
    std::mutex write_mutex;
    // writes the formatted spectra in input order as far as they are ready; the caller holds write_mutex
    const auto write_ready = [&]() {
      while (! failed.load() && next_write.load() < n && ready[next_write.load()].load(std::memory_order_acquire))
      {
        const Size k = next_write.load();
        if (errors[k])
        {
          failed = true;
          break;
        }
        try
        {
          streams[spectra[jobs[k]].getOriginalSpectrum().getMSLevel() - 1] << buffers[k];
        }
        catch (...)
        {
          errors[k] = std::current_exception();
          failed = true;
          break;
        }
        std::string().swap(buffers[k]);
        next_write.store(k + 1);
      }
    };
// the spectra are taken in input order (the write window relies on it): a nonmonotonic schedule (the OpenMP >= 5.0 default, e.g. the work
// stealing of clang/libomp) could hand out spectra far ahead of an unassigned one
#if defined(_OPENMP) && _OPENMP >= 201511
  #pragma omp parallel for schedule(monotonic: dynamic, 1)
#else
  #pragma omp parallel for schedule(dynamic, 1)
#endif
    for (SignedSize k = 0; k < (SignedSize)n; k++)
    {
      // wait (writing what is ready) while this spectrum is too far ahead of the next one to be written
      while ((Size)k >= next_write.load() + window && ! failed.load())
      {
        std::unique_lock<std::mutex> lock(write_mutex, std::try_to_lock);
        if (lock.owns_lock()) { write_ready(); }
        else { std::this_thread::yield(); }
      }
      if (! skip[k] && ! failed.load())
      {
        try
        {
          auto& dspec = spectra[jobs[k]];
          std::ostringstream os;
          uint index = first_index[k];
          double carry = first_carry[k];
          writeDeconvolvedMasses_(dspec, os, file_name, avg, decoy_avg, tols[dspec.getOriginalSpectrum().getMSLevel() - 1], write_detail, report_decoy,
                                  levelWeight_(noise_decoy_weights, dspec.getOriginalSpectrum().getMSLevel()), index, carry,
                                  abundances != nullptr && jobs[k] < abundances->size() ? &(*abundances)[jobs[k]] : nullptr);
          buffers[k] = os.str();
        }
        catch (...)
        {
          errors[k] = std::current_exception();
        }
      }
      ready[k].store(true, std::memory_order_release);
      std::unique_lock<std::mutex> lock(write_mutex, std::try_to_lock);
      if (lock.owns_lock()) { write_ready(); }
    }
    {
      std::lock_guard<std::mutex> lock(write_mutex);
      write_ready(); // what is not written yet, after the join
    }
    if (failed) { std::rethrow_exception(errors[next_write.load()]); }
  }

  void FLASHDeconvSpectrumFile::writeDeconvolvedMassesHeader(std::fstream& fs, const uint ms_level, const bool detail, const bool report_decoy,
                                                             const bool write_abundance)
  {
    if (detail)
    {
      if (ms_level == 1)
      {
        fs << "Index\tFileName\tScanNum\tFeatureIndex\t";
        if (report_decoy)
        {
          fs << "TargetDecoyType\t";
        }
        fs << "RetentionTime\tMassCountInSpec\tAverageMass\tMonoisotopicMass\t"
              "SumIntensity\tMinCharge\tMaxCharge\t"
              "PeakCount\tPeakMZs\tPeakIntensities\tPeakCharges\tPeakMasses\tPeakIsotopeIndices\tPeakPPMErrors\t"
              "NoisePeakMZs\tNoisePeakIntensities\tNoisePeakCharges\tNoisePeakMasses\tNoisePeakIsotopeIndices\tNoisePeakPPMErrors\t"
              "IsotopeCosine\tChargeCosine\tChargeScore\tMassSNR\tChargeSNR\tAveragePPMError\tRepresentativeCharge\tRepresentativeMzStart\tRepresentativeMzEnd\tQscore\tQscore2D\t";
        if (report_decoy)
        {
          fs << "Qvalue\t";
        }
        fs << "PerChargeIntensity\tPerIsotopeIntensity";
        if (write_abundance) { fs << "\tAbundance\tAbundanceTraceCount"; }
        fs << "\n";
      }
      else
      {
        fs << "Index\tFileName\tScanNum\tFeatureIndex\t";
        if (report_decoy)
        {
          fs << "TargetDecoyType\t";
        }
        fs << "RetentionTime\tMassCountInSpec\tAverageMass\tMonoisotopicMass\t"
              "SumIntensity\tMinCharge\tMaxCharge\t"
              "PeakCount\tPeakMZs\tPeakIntensities\tPeakCharges\tPeakMasses\tPeakIsotopeIndices\tPeakPPMErrors\t"
              "NoisePeakMZs\tNoisePeakIntensities\tNoisePeakCharges\tNoisePeakMasses\tNoisePeakIsotopeIndices\tNoisePeakPPMErrors\t"
              "PrecursorScanNum\tPrecursorFeatureIndex\tPrecursorMz\tPrecursorPeakIntensity\tPrecursorCharge\tPrecursorSNR\tPrecursorMonoisotopicMass\tPrecursorQscore\tPrecursorQscore2D\t"
              "DeconvolvedPrecursorPeakMzWithinIsolationWindow\tDeconvolvedPrecursorPeakIntensityWithinIsolationWindow\tPrecursorMassIntensityWithinIsolationWindow\tMassCountWithinIsolationWindow\tMassesWithinIsolationWindow\t"
              "MassIntensityRatiosWithinIsolationWindow\tTotalPeakIntensityWithinIsolationWindow\t";
        if (report_decoy)
        {
          fs << "PrecursorQvalue\tDecoyPrecursorMass\tDecoyPrecursorWeight\t";
        }
        fs << "IsotopeCosine\tChargeCosine\tChargeScore\tMassSNR\tChargeSNR\tAveragePPMError\tRepresentativeCharge\tRepresentativeMzStart\tRepresentativeMzEnd\tQscore\tQscore2D\t";
        if (report_decoy)
        {
          fs << "Qvalue\t";
        }
        fs << "PerChargeIntensity\tPerIsotopeIntensity\n";
      }
    }
    else
    {
      if (ms_level == 1)
      {
        fs << "Index\tFileName\tScanNum\tFeatureIndex\t";
        if (report_decoy)
        {
          fs << "TargetDecoyType\t";
        }
        fs << "RetentionTime\tMassCountInSpec\tAverageMass\tMonoisotopicMass\t"
              "SumIntensity\tMinCharge\tMaxCharge\t"
              "PeakCount\t"
              "IsotopeCosine\tChargeCosine\tChargeScore\tMassSNR\tChargeSNR\tAveragePPMError\tRepresentativeCharge\tRepresentativeMzStart\tRepresentativeMzEnd\tQscore\tQscore2D\t";
        if (report_decoy)
        {
          fs << "Qvalue";
        }
        if (write_abundance) { fs << (report_decoy ? "\t" : "") << "Abundance\tAbundanceTraceCount"; }
        fs << "\n";
      }
      else
      {
        fs << "Index\tFileName\tScanNum\tFeatureIndex\t";
        if (report_decoy)
        {
          fs << "TargetDecoyType\t";
        }
        fs << "RetentionTime\tMassCountInSpec\tAverageMass\tMonoisotopicMass\t"
              "SumIntensity\tMinCharge\tMaxCharge\t"
              "PeakCount\t"
              "PrecursorScanNum\tPrecursorFeatureIndex\tPrecursorMz\tPrecursorPeakIntensity\tPrecursorCharge\tPrecursorSNR\tPrecursorMonoisotopicMass\tPrecursorQscore\tPrecursorQscore2D\t"
              "DeconvolvedPrecursorPeakMzWithinIsolationWindow\tDeconvolvedPrecursorPeakIntensityWithinIsolationWindow\tPrecursorMassIntensityWithinIsolationWindow\tMassCountWithinIsolationWindow\tMassesWithinIsolationWindow\t"
              "MassIntensityRatiosWithinIsolationWindow\tTotalPeakIntensityWithinIsolationWindow\t";
        if (report_decoy)
        {
          fs << "PrecursorQvalue\tDecoyPrecursorMass\tDecoyPrecursorWeight\t";
        }
        fs << "IsotopeCosine\tChargeCosine\tChargeScore\tMassSNR\tChargeSNR\tAveragePPMError\tRepresentativeCharge\tRepresentativeMzStart\tRepresentativeMzEnd\tQscore\tQscore2D\t";
        if (report_decoy)
        {
          fs << "Qvalue\t";
        }
        fs << "\n";
      }
    }
  }

  void FLASHDeconvSpectrumFile::writeIsobaricQuantification(std::fstream& fs, std::vector<DeconvolvedSpectrum>& deconvolved_spectra)
  {
    std::stringstream ss;
    ss << "ScanNum\tPrecursorScanNum\tPrecursorFeatureIndex\tPrecursorMZ\tPrecursorPeakIntensity\tPrecursorCharge\tPrecursorSNR\tPrecursorMonoisotopicMass\tPrecursorQscore\t"
          "IsolationLeft\tIsolationRight\tDeconvolvedPrecursorPeakMzWithinIsolationWindow\tDeconvolvedPrecursorPeakIntensityWithinIsolationWindow\tPrecursorMassIntensityWithinIsolationWindow\tMassCountWithinIsolationWindow\tMassesWithinIsolationWindow\t"
          "MassIntensityRatiosWithinIsolationWindow\tTotalPeakIntensityWithinIsolationWindow\t"
          "RetentionTime\tMassCountInSpec\tQuantifiedChannelCount\tActivation\tActivationEnergy";

    Size channel_count = 0;
    for (auto& dspec : deconvolved_spectra)
    {
      if (dspec.isDecoy() || dspec.getOriginalSpectrum().getMSLevel() == 1)
        continue;
      auto quant = dspec.getQuantities();
      if (!quant.empty())
      {
        channel_count = quant.quantities.size();
        for (Size i = 0; i < quant.quantities.size(); i++)
        {
          ss << "\tQuantForCh" << (i + 1);
        }
        for (Size i = 0; i < quant.quantities.size(); i++)
        {
          ss << "\tNormalizedQuantForCh" << (i + 1);
        }
        for (Size i = 0; i < quant.quantities.size(); i++)
        {
          ss << "\tMergedQuantForCh" << (i + 1);
        }
        for (Size i = 0; i < quant.quantities.size(); i++)
        {
          ss << "\tNormalizedMergedQuantForCh" << (i + 1);
        }
        ss << "\n";
        break;
      }
    }
    if (channel_count == 0) return;

    for (auto& dspec : deconvolved_spectra)
    {
      bool precursor_found = ! dspec.getPrecursorPeakGroup().empty();
      if (dspec.isDecoy() || dspec.getOriginalSpectrum().getMSLevel() == 1) continue;
      int scan = dspec.getScanNumber();
      auto quant = dspec.getQuantities();
      // if (quant.empty())
      //   continue;
      int ch_count = 0;
      for (auto q : quant.quantities) if (q != 0) ch_count++;

      ss << scan << "\t" << dspec.getPrecursorScanNumber() << "\t" << (dspec.getPrecursorPeakGroup().getFeatureIndex() == 0 ? "NaN" : std::to_string(dspec.getPrecursorPeakGroup().getFeatureIndex()))
         << "\t" <<  std::to_string(dspec.getPrecursor().getMZ()) << "\t" <<  dspec.getPrecursor().getIntensity() << "\t" << (precursor_found ? std::to_string(dspec.getPrecursor().getCharge()) : "NaN") << "\t"
         << (precursor_found ? std::to_string(dspec.getPrecursorPeakGroup().getChargeSNR(dspec.getPrecursorPeakGroup().getRepAbsCharge())) : "NaN") << "\t"
         << (precursor_found ? std::to_string(dspec.getPrecursorPeakGroup().getMonoMass()) : "NaN") << "\t" << (precursor_found ? std::to_string(dspec.getPrecursorPeakGroup().getQscore2D()) : "NaN")
         << "\t" << dspec.getPrecursor().getMZ() - dspec.getPrecursor().getIsolationWindowLowerOffset() << "\t"
         << dspec.getPrecursor().getMZ() + dspec.getPrecursor().getIsolationWindowUpperOffset() << "\t";

      double total_int = 0;
      for (const auto& [m, it] : dspec.getPrecursorMassIntensityMap())
      {
        if (m <= 0)
        {
          total_int = it;
          break;
        }
      }

      if (precursor_found)
      {
        double pg_intensity = 0;
        std::string mass_str, int_str;
        windowMasses(dspec, total_int, pg_intensity, mass_str, int_str);

        ss << std::to_string(dspec.getDeconvolvedPrecursor().getMZ()) << "\t" << dspec.getDeconvolvedPrecursor().getIntensity() << "\t"
           << pg_intensity << "\t" << dspec.getPrecursorMassIntensityMap().size() - 1 << "\t" << mass_str << "\t" << int_str << "\t";
      }
      else
      {
        ss << "NaN\tNaN\tNaN\t0\tNaN\tNaN\t";
      }

      ss << total_int << "\t" << dspec.getOriginalSpectrum().getRT() << "\t" << countTargetMasses(dspec) << "\t" << ch_count << "\t"
         << (dspec.getActivationMethod() < Precursor::ActivationMethod::SIZE_OF_ACTIVATIONMETHOD
               ? Precursor::NamesOfActivationMethodShort[static_cast<Size>(dspec.getActivationMethod())]
               : "N/A")
         << "\t" << dspec.getPrecursor().getActivationEnergy();

      double sum = 0;

      if (quant.empty())
      {
        for (Size i = 0; i < channel_count * 4; i++)
        {
          ss << "\t" << 0;
        }
        ss << "\n";
      }
      else
      {
        for (auto q : quant.quantities)
        {
          ss << "\t" << std::to_string(q);
          sum += q;
        }
        for (auto q : quant.quantities)
        {
          ss << "\t" << std::to_string(q / sum);
        }
        sum = 0;
        for (auto q : quant.merged_quantities)
        {
          ss << "\t" << std::to_string(q);
          sum += q;
        }
        for (auto q : quant.merged_quantities)
        {
          ss << "\t" << std::to_string(q / sum);
        }
        ss << "\n";
      }
    }
    fs << ss.str();
  }

  void FLASHDeconvSpectrumFile::toDeconvolvedMap(const MSExperiment& map, std::vector<DeconvolvedSpectrum>& deconvolved_spectra, int mzml_charge,
                                                 const DoubleList& tols, MSExperiment& deconvolved_map)
  {
    // the experimental settings of the input map without its spectra (what a copy of it cleared of its spectra and chromatograms keeps but
    // the ranges, which the writer does not read), without copying the spectra
    deconvolved_map = static_cast<const ExperimentalSettings&>(map);
    // the deconvolved spectra in parallel (each from its own deconvolved spectrum), added below in the input order
    std::vector<MSSpectrum> deconvolved_mzMLs(deconvolved_spectra.size());
    std::vector<std::exception_ptr> errors(deconvolved_spectra.size());
#pragma omp parallel for schedule(dynamic, 16)
    for (SignedSize i = 0; i < (SignedSize)deconvolved_spectra.size(); i++)
    {
      auto& deconvolved_spectrum = deconvolved_spectra[i];
      if (deconvolved_spectrum.empty() || deconvolved_spectrum.isDecoy()) continue;
      try
      {
        deconvolved_mzMLs[i] = deconvolved_spectrum.toSpectrum(mzml_charge, tols[deconvolved_spectrum.getOriginalSpectrum().getMSLevel() - 1], false);
      }
      catch (...)
      {
        errors[i] = std::current_exception();
      }
    }
    for (Size i = 0; i < deconvolved_spectra.size(); i++)
    {
      if (errors[i]) std::rethrow_exception(errors[i]);
      if (deconvolved_mzMLs[i].empty()) continue;
      deconvolved_map.addSpectrum(std::move(deconvolved_mzMLs[i]));
    }
  }

  void FLASHDeconvSpectrumFile::writeMzML(const MSExperiment& map, std::vector<DeconvolvedSpectrum>& deconvolved_spectra, const std::string& deconvolved_mzML_file, const std::string& annotated_mzML_file,
                                          int mzml_charge, DoubleList tols)
  {
    if (deconvolved_mzML_file.empty() && annotated_mzML_file.empty())
      return;

    MSExperiment deconvolved_map;
    MSExperiment annotated_map;
    if (!deconvolved_mzML_file.empty())
    {
      toDeconvolvedMap(map, deconvolved_spectra, mzml_charge, tols, deconvolved_map);
    }
    if (!annotated_mzML_file.empty())
    {
      annotated_map = static_cast<const ExperimentalSettings&>(map);
      for (Size i = 0; i < deconvolved_spectra.size(); i++)
      {
        auto& deconvolved_spectrum = deconvolved_spectra[i];
        if (deconvolved_spectrum.empty()) continue;
        if (deconvolved_spectrum.isDecoy()) continue;
        auto anno_spec = MSSpectrum(deconvolved_spectrum.getOriginalSpectrum());
        anno_spec.sortByPosition();
        std::stringstream val {};

        for (auto& pg : deconvolved_spectrum)
        {
          if (pg.empty() || pg.getTargetDecoyType() != PeakGroup::TargetDecoyType::target)
          {
            continue;
          }
          val << std::to_string(pg.getMonoMass()) << ":";
          size_t pindex = 0;
          for (size_t k = 0; k < pg.size(); k++)
          {
            auto& p = pg[k];

            // Increment the pindex while it points to values less than p.mz
            while (pindex < anno_spec.size() - 1 && anno_spec[pindex].getMZ() < p.mz)
            {
              pindex++;
            }
            size_t nearest_index = pindex;
            if (anno_spec[pindex].getMZ() != p.mz)
            {
              if (pindex > 0 && std::abs(anno_spec[pindex - 1].getMZ() - p.mz) < std::abs(anno_spec[pindex].getMZ() - p.mz))
              {
                nearest_index = pindex - 1;
              }
            }

            val << nearest_index;

            if (k < pg.size() - 1)
            {
              val << ",";
            }
          }
          val << ";";
        }
        if (anno_spec.empty()) continue;
        anno_spec.setMetaValue("DeconvMassPeakIndices", val.str());
        annotated_map.addSpectrum(anno_spec);
      }
    }

    if (!deconvolved_mzML_file.empty())
    {
      MzMLFile mzml_file;
      mzml_file.store(deconvolved_mzML_file, deconvolved_map);
    }

    if (!annotated_mzML_file.empty())
    {
      MzMLFile mzml_file;
      mzml_file.store(annotated_mzML_file, annotated_map);
    }
  }

  void FLASHDeconvSpectrumFile::writeTopFDHeader(std::fstream& fs, const Param& param)
  {
    fs << "#FLASHDeconv generated msalign file\n";
    fs << "####################### Parameters ######################\n";
    for (const auto& p : param)
    {
      fs << "#" << p.name << ": " << p.value << "\n";
    }
    fs << "####################### Parameters ######################\n";
  }

  void FLASHDeconvSpectrumFile::writeTopFD(DeconvolvedSpectrum& dspec, std::fstream& fs, const std::string& filename, const double qval_threshold, const uint min_ms_level)
  {
    std::stringstream ss;
    UInt ms_level = dspec.getOriginalSpectrum().getMSLevel();
    if (ms_level > min_ms_level)
    {
      if (dspec.getPrecursorPeakGroup().empty() ||
          dspec.getPrecursorPeakGroup().getQvalue() > qval_threshold)
      {
        return;
      }
    }

    // the peak count limits and the Qscore cut below are taken over the target masses (with report_FDR the spectrum also holds decoys)
    const Size target_mass_count = countTargetMasses(dspec);
    if (target_mass_count < topFD_min_peak_count_)
    {
      return;
    }
    // the retention time in minutes, as TopPIC reads it (since TopPIC 1.8; seconds gave retention times 60 times too large)
    ss << std::fixed << std::setprecision(2);
    ss << "BEGIN IONS\n"
       << "FILE_NAME=" << filename << "\n"
       << "SPECTRUM_ID=" << dspec.getScanNumber() - 1 << "\n"
       << "TITLE=" << (dspec.isDecoy()? "DScan_" : "Scan_") << dspec.getScanNumber() << "\n"
       << "SCANS=" << dspec.getScanNumber() << "\n"
       << "RETENTION_TIME=" << std::setprecision(4) << dspec.getOriginalSpectrum().getRT() / 60.0 << std::setprecision(2) << "\n"
       << "LEVEL=" << dspec.getOriginalSpectrum().getMSLevel() << "\n";


    if (ms_level > 1)
    {
      double precursor_mass = dspec.getPrecursorPeakGroup().getMonoMass();

      ss << "MS_ONE_ID=" << dspec.getPrecursorScanNumber() - 1 << "\n"
         << "MS_ONE_SCAN=" << dspec.getPrecursorScanNumber() << "\n"
        // the isolation window around its target m/z (the offsets are relative to it; MzMLHandler keeps the selected ion m/z as the m/z)
        << "PRECURSOR_WINDOW_BEGIN=" << -dspec.getPrecursor().getIsolationWindowLowerOffset() + isolationTargetMz_(dspec.getPrecursor()) << "\n"
         << "PRECURSOR_WINDOW_END=" << dspec.getPrecursor().getIsolationWindowUpperOffset() + isolationTargetMz_(dspec.getPrecursor()) << "\n";

      // TopPIC's activation name (TopPIC stops at a name it does not know, e.g. ECD or EThcD)
      const std::string activation = topPICActivation_(dspec.getActivationMethod());
      if (!activation.empty())
      {
        ss << "ACTIVATION=" << activation << "\n";
      }
       ss << "PRECURSOR_MZ=" << std::to_string(dspec.getPrecursor().getMZ()) << "\n"
         << "PRECURSOR_CHARGE=" << (int)(dspec.getPrecursorCharge()) << "\n"
         << "PRECURSOR_MASS=" << std::to_string(precursor_mass) << "\n"
         << "PRECURSOR_INTENSITY=" << dspec.getPrecursor().getIntensity() << "\n"
         << "PRECURSOR_FEATURE_ID=" << dspec.getPrecursorPeakGroup().getFeatureIndex() << "\n"
         << "PRECURSOR_QSCORE=" << dspec.getPrecursorPeakGroup().getQscore2D() << "\n";
    }

    ss << std::setprecision(-1);

    double qscore_threshold = 0;
    std::vector<double> qscores;

    if (target_mass_count > topFD_max_peak_count_) // max peak count for TopPic = 500
    {
      qscores.reserve(target_mass_count);
      for (auto& pg : dspec)
      {
        if (pg.getTargetDecoyType() != PeakGroup::TargetDecoyType::target) continue;
        qscores.push_back(pg.getQscore2D());
      }
      std::sort(qscores.begin(), qscores.end());
      qscore_threshold = qscores[qscores.size() - topFD_max_peak_count_];
      std::vector<double>().swap(qscores);
    }

    int size = 0;
    for (auto& pg : dspec)
    {
      if (pg.getQscore2D() < qscore_threshold || pg.getTargetDecoyType() != PeakGroup::TargetDecoyType::target)
      {
        continue;
      }

      ss << std::fixed << std::setprecision(2);
      ss << std::to_string(pg.getMonoMass()) << "\t" << pg.getIntensity() << "\t"
         << (pg.isPositive() ? std::get<1>(pg.getAbsChargeRange()) : -std::get<1>(pg.getAbsChargeRange())) << "\n";
      ss << std::setprecision(-1);
      if (++size >= topFD_max_peak_count_)
      {
        break;
      }
    }
    ss << "END IONS\n\n";
    fs << ss.str();
  }
} // namespace OpenMS
