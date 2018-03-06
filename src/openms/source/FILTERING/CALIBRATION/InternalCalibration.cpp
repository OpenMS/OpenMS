// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------


#include <OpenMS/FILTERING/CALIBRATION/InternalCalibration.h>

#include <OpenMS/FORMAT/SVOutStream.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/SYSTEM/RWrapper.h>

#include <QtCore/QStringList>

#include <cstdio>

namespace OpenMS
{

  InternalCalibration::InternalCalibration()
    : ProgressLogger()
  {
  }

  void InternalCalibration::applyTransformation(std::vector<Precursor>& pcs, const MZTrafoModel& trafo)
  {
    // calibrate the precursor mass
    if (!pcs.empty())
    {
      for (Size i = 0; i < pcs.size(); ++i)
      {
        pcs[i].setMZ(trafo.predict(pcs[i].getMZ()));
      }
    }
  }

  void InternalCalibration::applyTransformation_(PeakMap::SpectrumType& spec, const MZTrafoModel& trafo)
  {
    typedef PeakMap::SpectrumType::Iterator SpecIt;

    // calibrate the spectrum itself
    for (SpecIt it = spec.begin(); it != spec.end(); ++it)
    {
      it->setMZ(trafo.predict(it->getMZ()));
    }
  }

  void InternalCalibration::applyTransformation(PeakMap::SpectrumType& spec, const IntList& target_mslvl, const MZTrafoModel& trafo)
  {
    // calibrate the peaks?
    if (ListUtils::contains(target_mslvl, spec.getMSLevel()))
    {
      applyTransformation_(spec, trafo);
    }
    // apply PC correction (only if target is MS1, and current spec is MS2; or target is MS2 and cs is MS3,...)
    if (ListUtils::contains(target_mslvl, spec.getMSLevel() - 1))
    {
      applyTransformation(spec.getPrecursors(), trafo);
    }     
  }

  void InternalCalibration::applyTransformation(PeakMap& exp, const IntList& target_mslvl, const MZTrafoModel& trafo)
  {
    for (PeakMap::Iterator it = exp.begin(); it != exp.end(); ++it)
    {
      applyTransformation(*it, target_mslvl, trafo);
    }
  }

  Size InternalCalibration::fillCalibrants(const PeakMap exp,
                                           const std::vector<InternalCalibration::LockMass>& ref_masses,
                                           double tol_ppm,
                                           bool lock_require_mono,
                                           bool lock_require_iso,
                                           CalibrationData& failed_lock_masses,
                                           bool verbose /*= true*/)
  {
    cal_data_.clear();

    //
    // find lock masses in data and build calibrant table
    //
    std::map<Size, Size> stats_cal_per_spectrum;
    typedef PeakMap::ConstIterator ExpCIt;
    for (ExpCIt it = exp.begin(); it != exp.end(); ++it)
    {
      // empty spectrum
      if (it->empty()) {
        ++stats_cal_per_spectrum[0];
        continue;
      }

      Size cnt_cd = cal_data_.size();
      // iterate over calibrants
      for (std::vector<InternalCalibration::LockMass>::const_iterator itl = ref_masses.begin(); itl != ref_masses.end(); ++itl)
      {
        // calibrant meant for this MS level?
        if (it->getMSLevel() != itl->ms_level) continue;

        Size s = it->findNearest(itl->mz);
        const double mz_obs = (*it)[s].getMZ();
        if (Math::getPPMAbs(mz_obs, itl->mz) > tol_ppm)
        {
          failed_lock_masses.insertCalibrationPoint(it->getRT(), itl->mz, 0.0, itl->mz, 0.0, std::distance(ref_masses.begin(), itl));
        }
        else
        {
          if (lock_require_mono)
          {
            // check if its the monoisotopic .. discard otherwise
            const double mz_iso_left = mz_obs - (Constants::C13C12_MASSDIFF_U / itl->charge);
            Size s_left = it->findNearest(mz_iso_left);
            if (Math::getPPMAbs(mz_iso_left, (*it)[s_left].getMZ()) < 0.5) // intra-scan ppm should be very good!
            { // peak nearby lock mass was not the monoisotopic
              if (verbose) LOG_INFO << "peak at [RT, m/z] " << it->getRT() << ", " << (*it)[s].getMZ() << " is NOT monoisotopic. Skipping it!\n";
              failed_lock_masses.insertCalibrationPoint(it->getRT(), itl->mz, 1.0, itl->mz, 0.0, std::distance(ref_masses.begin(), itl));
              continue;
            }
          }
          if (lock_require_iso)
          {
            // require it to have a +1 isotope?!
            const double mz_iso_right = mz_obs + Constants::C13C12_MASSDIFF_U / itl->charge;
            Size s_right = it->findNearest(mz_iso_right);
            if (!(Math::getPPMAbs(mz_iso_right, (*it)[s_right].getMZ()) < 0.5)) // intra-scan ppm should be very good!
            { // peak has no +1iso.. weird
              if (verbose) LOG_INFO << "peak at [RT, m/z] " << it->getRT() << ", " << (*it)[s].getMZ() << " has no +1 isotope (ppm to closest: " << Math::getPPM(mz_iso_right, (*it)[s_right].getMZ()) << ")... Skipping it!\n";
              failed_lock_masses.insertCalibrationPoint(it->getRT(), itl->mz, 2.0, itl->mz, 0.0, std::distance(ref_masses.begin(), itl));
              continue;
            }
          }
          cal_data_.insertCalibrationPoint(it->getRT(), mz_obs, (*it)[s].getIntensity(), itl->mz, std::log((*it)[s].getIntensity()), std::distance(ref_masses.begin(), itl));
        }
      }
      // how many locks found in this spectrum?!
      ++stats_cal_per_spectrum[cal_data_.size()-cnt_cd];
    }

    LOG_INFO << "Lock masses found across viable spectra:\n";
    for (std::map<Size, Size>::const_iterator its = stats_cal_per_spectrum.begin(); its != stats_cal_per_spectrum.end(); ++its)
    {
      LOG_INFO << "  " << its->first << " [of " << ref_masses.size() << "] lock masses: " << its->second << "x\n";
    }
    LOG_INFO << std::endl;

    // sort CalData by RT
    cal_data_.sortByRT();

    return cal_data_.size();
  }

  Size InternalCalibration::fillCalibrants( const FeatureMap& fm, double tol_ppm )
  {
    cal_data_.clear();
    for (FeatureMap::ConstIterator it = fm.begin(); it != fm.end(); ++it)
    {
      const std::vector<PeptideIdentification>& ids = it->getPeptideIdentifications();
      if (ids.empty() || ids[0].empty()) continue;

      PeptideIdentification pid = ids[0];
      pid.sort();
      double mz_ref = pid.getHits()[0].getSequence().getMonoWeight(OpenMS::Residue::Full, pid.getHits()[0].getCharge());
      if (tol_ppm < Math::getPPMAbs(it->getMZ(), mz_ref)) continue;
      cal_data_.insertCalibrationPoint(it->getRT(), it->getMZ(), it->getIntensity(), mz_ref, log(it->getIntensity()));
    }

    // unassigned peptide IDs
    fillIDs_(fm.getUnassignedPeptideIdentifications(), tol_ppm);

    LOG_INFO << "Found " << cal_data_.size() << " calibrants (incl. unassigned) in FeatureMap." << std::endl;

    // sort CalData by RT
    cal_data_.sortByRT();

    return cal_data_.size();
  }

  void InternalCalibration::fillIDs_( const std::vector<PeptideIdentification>& pep_ids, double tol_ppm )
  {
    Size cnt_nomz(0);
    Size cnt_nort(0);

    for (std::vector<PeptideIdentification>::const_iterator it = pep_ids.begin(); it != pep_ids.end(); ++it)
    {
      if (it->empty()) continue;
      if (!it->hasMZ())
      {
        ++cnt_nomz;
        continue;
      }
      if (!it->hasRT())
      {
        ++cnt_nort;
        continue;
      }
      PeptideIdentification pid = *it;
      pid.sort();
      int q = pid.getHits()[0].getCharge();
      double mz_ref = pid.getHits()[0].getSequence().getMonoWeight(OpenMS::Residue::Full, q) / q;
      if (tol_ppm < Math::getPPMAbs(it->getMZ(), mz_ref)) continue;

      const double weight = 1.0;
      const double intensity = 1.0;
      cal_data_.insertCalibrationPoint(it->getRT(), it->getMZ(), intensity, mz_ref, weight);
    }
    LOG_INFO << "Found " << cal_data_.size() << " calibrants in peptide IDs." << std::endl;
    if (cnt_nomz > 0) LOG_WARN << "Warning: " << cnt_nomz << "/" << pep_ids.size() << " were skipped, since they have no m/z value set! They cannot be used as calibration point." << std::endl;
    if (cnt_nort > 0) LOG_WARN << "Warning: " << cnt_nort << "/" << pep_ids.size() << " were skipped, since they have no RT value set! They cannot be used as calibration point." << std::endl;
  }

  Size InternalCalibration::fillCalibrants( const std::vector<PeptideIdentification>& pep_ids, double tol_ppm )
  {
    cal_data_.clear();
    fillIDs_(pep_ids, tol_ppm);
    // sort CalData by RT
    cal_data_.sortByRT();

    return cal_data_.size();
  }

  const CalibrationData& InternalCalibration::getCalibrationPoints() const
  {
    return cal_data_;
  }

  bool InternalCalibration::calibrate(PeakMap& exp, 
                                      const IntList& target_mslvl,
                                      MZTrafoModel::MODELTYPE model_type,
                                      double rt_chunk,
                                      bool use_RANSAC,
                                      double post_ppm_median,
                                      double post_ppm_MAD,
                                      const String& file_models,
                                      const String& file_models_plot,
                                      const String& file_residuals,
                                      const String& file_residuals_plot,
                                      const String& rscript_executable_)
  {
    QString rscript_executable = rscript_executable_.toQString();

    // ensure sorting; required for finding RT ranges and lock masses
    if (!exp.isSorted(true))
    {
      exp.sortSpectra(true);
    }

    startProgress(0, exp.size(), "Applying calibration to data");

    std::vector<MZTrafoModel> tms; // each spectrum gets its own model (params are cheap to store)
    std::map<Size, Size> invalid_models; // indices from tms[] -> exp[]; where model creation failed (e..g, not enough calibration points)
    bool hasValidModels(false); // was at least one model valid?
    bool global_model = (rt_chunk < 0);
    if (global_model)
    { // build one global modal
      LOG_INFO << "Building a global model..." << std::endl;
      tms.push_back(MZTrafoModel());
      tms[0].train(cal_data_, model_type, use_RANSAC);
      if (MZTrafoModel::isValidModel(tms[0]))
      {
        applyTransformation(exp, target_mslvl, tms[0]);
        hasValidModels = true;
      }
    } else
    { // one model per spectrum (not all might be needed, if certain MS levels are excluded from calibration)
      tms.reserve(exp.size());
      // go through spectra and calibrate
      Size i(0), i_mslvl(0);
      for (PeakMap::Iterator it = exp.begin(); it != exp.end(); ++it, ++i)
      {
        setProgress(i);

        // skip this MS level?
        if (!(ListUtils::contains(target_mslvl, it->getMSLevel()) ||     // scan m/z needs correction
              ListUtils::contains(target_mslvl, it->getMSLevel() - 1)))  // precursor m/z needs correction
        {
          continue;
        }

        //
        // build model
        //
        tms.push_back(MZTrafoModel());
        tms.back().train(cal_data_, model_type, use_RANSAC, it->getRT() - rt_chunk, it->getRT() + rt_chunk);
        if (!MZTrafoModel::isValidModel(tms.back())) // model not trained or coefficients are too extreme
        {
          invalid_models[i_mslvl] = i;
        }
        else
        {
          applyTransformation(*it, target_mslvl, tms.back());
        }
        ++i_mslvl;
      } // MSExp::iter

      //////////////////////////////////////////////////////////////////////////
      // CHECK Models -- use neighbors if needed
      //////////////////////////////////////////////////////////////////////////

      hasValidModels = (std::find_if(tms.begin(), tms.end(), MZTrafoModel::isValidModel) != tms.end());
      // did we build any model at all?
      if (hasValidModels && !invalid_models.empty())
      {
        // 2nd attempt to calibrate spectra using neighboring models
        // (will not be entered for global model since could_not_cal is empty)
        LOG_INFO << "\nCalibration failed on " << invalid_models.size() << "/" << tms.size() << " [" <<  invalid_models.size() * 100 / tms.size() << " %] spectra. "
          << "Using the closest successful model on these." << std::endl;

        std::vector<MZTrafoModel> tms_new = tms; // will contain corrected models (this wastes a bit of memory)
        for (std::map<Size, Size>::const_iterator it = invalid_models.begin(); it != invalid_models.end(); ++it)
        {
          Size p = it->first;
          // find model closest valid model to p'th model
          std::vector<MZTrafoModel>::iterator it_center_r = tms.begin() + p; // points to 'p'
          std::vector<MZTrafoModel>::iterator it_right = std::find_if(it_center_r, tms.end(), MZTrafoModel::isValidModel);
          std::vector<MZTrafoModel>::reverse_iterator it_center_l = tms.rbegin() + (tms.size() - p - 1); // points to 'p'
          std::vector<MZTrafoModel>::reverse_iterator it_left = std::find_if(it_center_l, tms.rend(), MZTrafoModel::isValidModel);
          Size dist_right(0), dist_left(0);
          if (it_right != tms.end()) dist_right = std::distance(it_center_r, it_right);
          if (it_left != tms.rend()) dist_left  = std::distance(it_center_l, it_left);
          Size model_index;
          if (((dist_left <= dist_right) || dist_right == 0) && dist_left != 0) // left is closer in #spectra, i.e. time; or is the only valid direction
          {
            model_index = p - dist_left;
          } 
          else
          {
            model_index = p + dist_right;
          }
          applyTransformation(exp[it->second], target_mslvl, tms[model_index]);
          tms_new[p].setCoefficients(tms[model_index]); // overwrite invalid model
        }
        tms_new.swap(tms);
        // consistency check: all models must be valid at this point
        for (Size i = 0; i < tms.size(); ++i) if (!MZTrafoModel::isValidModel(tms[i])) throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "InternalCalibration::calibrate(): Internal error. Not all models are valid!", String(i));
      }
    }
    endProgress();

    // check if Rscript is available
    if (!file_models_plot.empty() || !file_residuals_plot.empty())
    {
      if (!RWrapper::findR(rscript_executable, true))
      {
        LOG_ERROR << "The R interpreter is required to create PNG plot files. To avoid the error, either do not request 'quality_control:*_plot' (not recommended) or fix your R installation." << std::endl;
        return false;
      }
    }

    //
    // write the model parameters to file and/or plot them
    //
    if (!file_models.empty() || !file_models_plot.empty())
    {
      String out_table = File::getTemporaryFile(file_models);
      { // we need this scope, to ensure that SVOutStream writes its cache, before we call RWrapper!
        SVOutStream sv(out_table, ", ", ", ", String::NONE);

        sv << "# model parameters (for all successfully trained models)" << nl
          << "RT" << "A (offset)" << "B (slope)" << "C (power)" << "source" << nl;
        for (Size i = 0; i < tms.size(); ++i)
        {
          sv << tms[i].getRT() << tms[i].toString();
          if (!MZTrafoModel::isValidModel(tms[i])) sv << "invalid"; // this only happens if ALL models are invalid (since otherwise they would use 'neighbour')
          else if (invalid_models.count(i) > 0) sv << "neighbor";
          else sv << "local";
          sv << nl;
        }
      }
      
      // plot it
      if (!file_models_plot.empty())
      {
        if (!RWrapper::runScript("InternalCalibration_Models.R", QStringList() << out_table.toQString() << file_models_plot.toQString(), rscript_executable))
        {
          LOG_ERROR << "R script failed. To avoid the error, either disable the creation of 'quality_control:models_plot' (not recommended) or fix your R installation." << std::endl;
          return false;
        }
      }

    }

    //
    // plot the residual error (after calibration)
    // go through Calibration data points
    //
    SVOutStream* sv = nullptr;      
    String out_table_residuals;
    if (!file_residuals.empty() || !file_residuals_plot.empty())
    {
      out_table_residuals = File::getTemporaryFile(file_residuals);
      sv = new SVOutStream(out_table_residuals, ", ", ", ", String::NONE);
    }

    std::vector<double> vec_ppm_before, vec_ppm_after;
    vec_ppm_before.reserve(cal_data_.size());
    vec_ppm_after.reserve(cal_data_.size());
    if (sv != nullptr) *sv << "# residual error after calibration" << nl
                        << "RT" << "intensity" << "mz ref" << "mz before" << "mz after" << "ppm before" << "ppm after" << nl;
    Size ii(0);
    for (CalibrationData::const_iterator itc = cal_data_.begin(); itc != cal_data_.end(); ++itc, ++ii)
    {
      double rt = itc->getRT();
      // find closest model in RT
      Size idx = (global_model ? 0 : MZTrafoModel::findNearest(tms, rt));

      double mz_corrected = std::numeric_limits<double>::quiet_NaN();
      if (MZTrafoModel::isValidModel(tms[idx])) mz_corrected = tms[idx].predict(itc->getMZ());
      double mz_ref = cal_data_.getRefMZ(ii);
      double ppm_before = Math::getPPM(itc->getMZ(), mz_ref);
      double ppm_after = Math::getPPM(mz_corrected, mz_ref);
      vec_ppm_before.push_back(ppm_before);
      vec_ppm_after.push_back(ppm_after);
      if (sv != nullptr)
      {
        *sv << rt 
            << itc->getIntensity()
            << mz_ref
            << itc->getMZ();
        sv->writeValueOrNan(mz_corrected)
            << ppm_before;
        sv->writeValueOrNan(ppm_after)
            << nl;
      }
    }
    delete sv;

    // plot it
    if (!file_residuals_plot.empty())
    {
      if (!RWrapper::runScript("InternalCalibration_Residuals.R", QStringList() << out_table_residuals.toQString() << file_residuals_plot.toQString(), rscript_executable))
      {
        LOG_ERROR << "R script failed. To avoid the error, either disable the creation of 'quality_control:residuals_plot' (not recommended) or fix your R installation." << std::endl;
        return false;
      }
    }


    if (!hasValidModels)
    { // QC tables are done; quit
      LOG_ERROR << "Error: Could not build a single local calibration model! Check your calibrants and/or extend the search window!" << std::endl;
      if (use_RANSAC) LOG_ERROR << "       Since you are using RANSAC, check the parameters as well and test different setups." << std::endl;

      return false;
    }

    // use median and MAD to ignore outliers
    double median_ppm_before = Math::median(vec_ppm_before.begin(), vec_ppm_before.end());
    double MAD_ppm_before =  Math::MAD(vec_ppm_before.begin(), vec_ppm_before.end(), median_ppm_before);
    LOG_INFO << "\n-----\n" <<
      "ppm stats before calibration: median = " << median_ppm_before << "  MAD = " << MAD_ppm_before << "\n";
    double median_ppm_after = Math::median(vec_ppm_after.begin(), vec_ppm_after.end());
    double MAD_ppm_after =  Math::MAD(vec_ppm_after.begin(), vec_ppm_after.end(), median_ppm_after);
    LOG_INFO << "ppm stats after calibration: median = " << median_ppm_after << "  MAD = " << MAD_ppm_after << "\n";

    // check desired limits
    if (post_ppm_median < fabs(median_ppm_after))
    {
      LOG_INFO << "Post calibration median threshold (" << post_ppm_median << " ppm) not reached (median = |" << median_ppm_after << "| ppm). Failed to calibrate!" << std::endl;
      return false;
    }
    if (post_ppm_MAD < fabs(MAD_ppm_after))
    {
      LOG_INFO << "Post calibration median threshold (" << post_ppm_MAD << " ppm) not reached (median = |" << MAD_ppm_after << "| ppm). Failed to calibrate!" << std::endl;
      return false;
    }


    return true; // success
  }

}
