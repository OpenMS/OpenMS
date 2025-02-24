// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------


#include <OpenMS/PROCESSING/CALIBRATION/InternalCalibration.h>

#include <OpenMS/FORMAT/SVOutStream.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/SYSTEM/File.h>
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
        pcs[i].setMetaValue("mz_raw", pcs[i].getMZ());
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

  Size InternalCalibration::fillCalibrants(const PeakMap& exp,
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
              if (verbose)
              {
                OPENMS_LOG_INFO << "peak at [RT, m/z] " << it->getRT() << ", " << (*it)[s].getMZ() << " is NOT monoisotopic. Skipping it!\n";
              }
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
              if (verbose)
              {
                OPENMS_LOG_INFO << "peak at [RT, m/z] " << it->getRT() << ", " << (*it)[s].getMZ() << " has no +1 isotope (ppm to closest: " << Math::getPPM(mz_iso_right, (*it)[s_right].getMZ()) << ")... Skipping it!\n";
              }
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

    OPENMS_LOG_INFO << "Lock masses found across viable spectra:\n";
    for (std::map<Size, Size>::const_iterator its = stats_cal_per_spectrum.begin(); its != stats_cal_per_spectrum.end(); ++its)
    {
      OPENMS_LOG_INFO << "  " << its->first << " [of " << ref_masses.size() << "] lock masses: " << its->second << "x\n";
    }
    OPENMS_LOG_INFO << std::endl;

    // sort CalData by RT
    cal_data_.sortByRT();

    return cal_data_.size();
  }

  Size InternalCalibration::fillCalibrants(const FeatureMap& fm, double tol_ppm)
  {
    cal_data_.clear();
    CalibrantStats_ stats(tol_ppm);
    stats.cnt_total = fm.size() + fm.getUnassignedPeptideIdentifications().size();

    for (const auto& f : fm)
    {
      const std::vector<PeptideIdentification>& ids = f.getPeptideIdentifications();
      double mz_ref;
      if (ids.empty())
      {
        continue;
      }
      if (isDecalibrated_(ids[0], f.getMZ(), tol_ppm, stats, mz_ref))
      {
        continue;
      }
      cal_data_.insertCalibrationPoint(f.getRT(), f.getMZ(), f.getIntensity(), mz_ref, log(f.getIntensity()));
    }

    // unassigned peptide IDs
    fillIDs_(fm.getUnassignedPeptideIdentifications(), tol_ppm, stats);

    OPENMS_LOG_INFO << "Found " << cal_data_.size() << " calibrants (incl. unassigned) in FeatureMap." << std::endl;
    stats.print();

    // sort CalData by RT
    cal_data_.sortByRT();

    return cal_data_.size();
  }

  void InternalCalibration::fillID_(const PeptideIdentification& pep_id, const double tol_ppm, CalibrantStats_& stats)
  {
    if (pep_id.empty())
    {
      ++stats.cnt_empty;
      return;
    }
    if (!pep_id.hasMZ())
    {
      ++stats.cnt_nomz;
      return;
    }
    if (!pep_id.hasRT())
    {
      ++stats.cnt_nort;
      return;
    }
    double mz_ref;
    if (isDecalibrated_(pep_id, pep_id.getMZ(), tol_ppm, stats, mz_ref))
    {
      return;
    }

    cal_data_.insertCalibrationPoint(pep_id.getRT(), pep_id.getMZ(), 1.0, mz_ref, 1.0);
  }

  void InternalCalibration::fillIDs_( const std::vector<PeptideIdentification>& pep_ids, const double tol_ppm, CalibrantStats_& stats)
  {
    for (const auto& id : pep_ids)
    {
      fillID_(id, tol_ppm, stats);
    }
 }

  bool InternalCalibration::isDecalibrated_(const PeptideIdentification& pep_id, const double mz_obs, const double tol_ppm, CalibrantStats_& stats, double& mz_ref)
  {
    PeptideIdentification pid = pep_id;
    pid.sort();
    int q = pid.getHits()[0].getCharge();
    mz_ref = pid.getHits()[0].getSequence().getMZ(q);

    // Only use ID if precursor m/z and theoretical mass don't deviate too much.
    // as they may occur due to isotopic peak misassignments
    double delta = Math::getPPMAbs(mz_obs, mz_ref);
    if (tol_ppm < delta)
    {
      if (stats.cnt_decal < 10)
      {
        OPENMS_LOG_INFO << "Peptide " << pid.getHits()[0].getSequence().toString() << " is " << delta << " (>" << tol_ppm << ") ppm away from theoretical mass and is omitted as calibration point.\n";
      }
      else if (stats.cnt_decal == 10)
      {
        OPENMS_LOG_INFO << "More than 10 peptides are at least " << tol_ppm << " ppm away from theoretical mass and are omitted as calibration point.";
      }
      ++stats.cnt_decal;
      return true;
    }
    return false;
  }

  Size InternalCalibration::fillCalibrants(const std::vector<PeptideIdentification>& pep_ids, double tol_ppm)
  {
    cal_data_.clear();
    CalibrantStats_ stats(tol_ppm);
    stats.cnt_total = pep_ids.size();
    fillIDs_(pep_ids, tol_ppm, stats);
    OPENMS_LOG_INFO << "Found " << cal_data_.size() << " calibrants in peptide IDs." << std::endl;
    stats.print();

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
      OPENMS_LOG_INFO << "Building a global model..." << std::endl;
      tms.emplace_back();
      tms[0].train(cal_data_, model_type, use_RANSAC);
      if (MZTrafoModel::isValidModel(tms[0]))
      {
        applyTransformation(exp, target_mslvl, tms[0]);
        hasValidModels = true;
      }
    }
    else
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
        tms.emplace_back();
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
        OPENMS_LOG_INFO << "\nCalibration failed on " << invalid_models.size() << "/" << tms.size() << " [" <<  invalid_models.size() * 100 / tms.size() << " %] spectra. "
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
          if (it_right != tms.end())
          {
            dist_right = std::distance(it_center_r, it_right);
          }
          if (it_left != tms.rend())
          {
            dist_left  = std::distance(it_center_l, it_left);
          }
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
        for (Size j = 0; j < tms.size(); ++j)
        {
          if (!MZTrafoModel::isValidModel(tms[j]))
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "InternalCalibration::calibrate(): Internal error. Not all models are valid!", String(j));
          }
        }
      }
    }
    endProgress();

    // check if Rscript is available
    if (!file_models_plot.empty() || !file_residuals_plot.empty())
    {
      if (!RWrapper::findR(rscript_executable, true))
      {
        OPENMS_LOG_ERROR << "The R interpreter is required to create PNG plot files. To avoid the error, either do not request 'quality_control:*_plot' (not recommended) or fix your R installation." << std::endl;
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
          if (!MZTrafoModel::isValidModel(tms[i]))
          {
            sv << "invalid"; // this only happens if ALL models are invalid (since otherwise they would use 'neighbour')
          }
          else if (invalid_models.count(i) > 0)
          {
            sv << "neighbor";
          }
          else 
          {
            sv << "local";
          }
          sv << nl;
        }
      }
      
      // plot it
      if (!file_models_plot.empty())
      {
        if (!RWrapper::runScript("InternalCalibration_Models.R", QStringList() << out_table.toQString() << file_models_plot.toQString(), rscript_executable))
        {
          OPENMS_LOG_ERROR << "R script failed. To avoid the error, either disable the creation of 'quality_control:models_plot' (not recommended) or fix your R installation." << std::endl;
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
    if (sv != nullptr)
    {
      *sv << "# residual error after calibration" << nl
                        << "RT" << "intensity" << "mz ref" << "mz before" << "mz after" << "ppm before" << "ppm after" << nl;
    }
    Size ii(0);
    for (CalibrationData::const_iterator itc = cal_data_.begin(); itc != cal_data_.end(); ++itc, ++ii)
    {
      double rt = itc->getRT();
      // find closest model in RT
      Size idx = (global_model ? 0 : MZTrafoModel::findNearest(tms, rt));

      double mz_corrected = std::numeric_limits<double>::quiet_NaN();
      if (MZTrafoModel::isValidModel(tms[idx]))
      {
        mz_corrected = tms[idx].predict(itc->getMZ());
      }
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
        OPENMS_LOG_ERROR << "R script failed. To avoid the error, either disable the creation of 'quality_control:residuals_plot' (not recommended) or fix your R installation." << std::endl;
        return false;
      }
    }


    if (!hasValidModels)
    { // QC tables are done; quit
      OPENMS_LOG_ERROR << "Error: Could not build a single local calibration model! Check your calibrants and/or extend the search window!" << std::endl;
      if (use_RANSAC)
      {
        OPENMS_LOG_ERROR << "       Since you are using RANSAC, check the parameters as well and test different setups." << std::endl;
      }
      return false;
    }

    // use median and MAD to ignore outliers
    double median_ppm_before = Math::median(vec_ppm_before.begin(), vec_ppm_before.end());
    double MAD_ppm_before =  Math::MAD(vec_ppm_before.begin(), vec_ppm_before.end(), median_ppm_before);
    OPENMS_LOG_INFO << "\n-----\n" <<
      "ppm stats before calibration: median = " << median_ppm_before << "  MAD = " << MAD_ppm_before << "\n";
    double median_ppm_after = Math::median(vec_ppm_after.begin(), vec_ppm_after.end());
    double MAD_ppm_after =  Math::MAD(vec_ppm_after.begin(), vec_ppm_after.end(), median_ppm_after);
    OPENMS_LOG_INFO << "ppm stats after calibration: median = " << median_ppm_after << "  MAD = " << MAD_ppm_after << "\n";

    // check desired limits
    if (post_ppm_median < fabs(median_ppm_after))
    {
      OPENMS_LOG_INFO << "Post calibration median threshold (" << post_ppm_median << " ppm) not reached (median = |" << median_ppm_after << "| ppm). Failed to calibrate!" << std::endl;
      return false;
    }
    if (post_ppm_MAD < fabs(MAD_ppm_after))
    {
      OPENMS_LOG_INFO << "Post calibration MAD threshold (" << post_ppm_MAD << " ppm) not reached (MAD = |" << MAD_ppm_after << "| ppm). Failed to calibrate!" << std::endl;
      return false;
    }


    return true; // success
  }

}
