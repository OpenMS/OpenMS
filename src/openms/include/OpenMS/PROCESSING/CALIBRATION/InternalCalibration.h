// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/CalibrationData.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/PROCESSING/CALIBRATION/MZTrafoModel.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <vector>

namespace OpenMS
{
  
  class FeatureMap;
    
  /**
    @brief A mass recalibration method using linear/quadratic interpolation (robust/weighted) of given reference masses.

    @ingroup SignalProcessing
  */
 
  class OPENMS_DLLAPI InternalCalibration
    : public ProgressLogger
  {
  public:

    /// Default constructor
    InternalCalibration();

    /// Destructor
    ~InternalCalibration() override{}

    /// helper class, describing a lock mass
    struct LockMass
    {
      double mz; ///< m/z of the lock mass (incl. adducts)
      unsigned int ms_level;   ///< MS level where it occurs
      int charge;     ///< charge of the ion (to find isotopes)

      LockMass(double mz_, int lvl_, int charge_)
        : mz(mz_),
         ms_level(lvl_),
         charge(charge_)
      {}
    };

    /** 
      @brief Extract calibrants from Raw data (mzML)

      Lock masses are searched in each spectrum and added to the internal calibrant database.

      Filters can be used to exclude spurious peaks, i.e. require the calibrant peak to be monoisotopic or
      to have a +1 isotope (should not be used for very low abundant calibrants).
      If a calibrant is not found, it is added to a 'failed_lock_masses' database which is returned and not stored internally.
      The intensity of the peaks describe the reason for failed detection: 0.0 - peak not found with the given ppm tolerance;
      1.0 - peak is not monoisotopic (can only occur if 'lock_require_mono' is true)
      2.0 - peak has no +1 isotope (can only occur if 'lock_require_iso' is true)

      @param exp Peak map containing the lock masses
      @param ref_masses List of lock masses
      @param tol_ppm Search window for lock masses in 'exp'
      @param lock_require_mono Require that a lock mass is the monoisotopic peak (i.e. not an isotope peak) -- lock mass is rejected otherwise
      @param lock_require_iso Require that a lock mass has isotope peaks to its right -- lock mass is rejected otherwise
      @param failed_lock_masses Set of calibration masses which were not found, i.e. their expected m/z and RT positions;
      @param verbose Print information on 'lock_require_XXX' matches during search
      @return Number of calibration masses found

    */
    Size fillCalibrants(const PeakMap& exp,
                        const std::vector<InternalCalibration::LockMass>& ref_masses,
                        double tol_ppm,
                        bool lock_require_mono,
                        bool lock_require_iso,
                        CalibrationData& failed_lock_masses,
                        bool verbose = true);

    /** 
      @brief Extract calibrants from identifications

      Extracts only the first hit from the first peptide identification of each feature.
      Hits are sorted beforehand.
      Ambiguities should be resolved before, e.g. using IDFilter.
      RT and m/z are taken from the features, not from the identifications (for an exception see below)!

      Unassigned peptide identifications are also taken into account!
      RT and m/z are naturally taken from the IDs, since to feature is assigned.
      If you do not want these IDs, remove them from the feature map before calling this function.

      A filtering step is done in the m/z dimension using @p tol_ppm.
      Since precursor masses could be annotated wrongly (e.g. isotope peak instead of mono),
      larger outliers are removed before accepting an ID as calibrant.

      @param fm FeatureMap with peptide identifications
      @param tol_ppm Only accept ID's whose theoretical mass deviates at most this much from annotated
      @return Number of calibration masses found

    */
    Size fillCalibrants(const FeatureMap& fm, double tol_ppm);

    /** 
      @brief Extract calibrants from identifications

      Extracts only the first hit from each peptide identification.
      Hits are sorted beforehand.
      Ambiguities should be resolved before, e.g. using IDFilter.

      A filtering step is done in the m/z dimension using @p tol_ppm.
      Since precursor masses could be annotated wrongly (e.g. isotope peak instead of mono),
      larger outliers are removed before accepting an ID as calibrant.

      @param pep_ids Peptide ids (e.g. from an idXML file)
      @param tol_ppm Only accept ID's whose theoretical mass deviates at most this much from annotated
      @return Number of calibration masses found
    */
    Size fillCalibrants(const std::vector<PeptideIdentification>& pep_ids, double tol_ppm);

    /**
      @brief Get container of calibration points

      Filled using fillCalibrants() methods.

      @return Container of calibration points

    */
    const CalibrationData& getCalibrationPoints() const;

    /**
      @brief Apply calibration to data

      For each spectrum, a calibration model will be computed and applied.
      Make sure to call fillCalibrants() before, so a model can be created.

      The MSExperiment will be sorted by RT and m/z if unsorted.

      @param exp MSExperiment holding the Raw data to calibrate
      @param target_mslvl MS-levels where calibration should be applied to
      @param model_type Linear or quadratic model; select based on your instrument
      @param rt_chunk RT-window size (one-sided) of calibration points to collect around each spectrum. 
             Set to negative values, to build one global model instead.
      @param use_RANSAC Remove outliers before fitting a model?!
      @param post_ppm_median The median ppm error of the calibrants must be at least this good after calibration; otherwise this method returns false(fail)
      @param post_ppm_MAD The median absolute deviation of the calibrants must be at least this good after calibration; otherwise this method returns false(fail)
      @param file_models Output CSV filename, where model parameters are written to (pass empty string to skip)
      @param file_models_plot Output PNG image model parameters (pass empty string to skip)
      @param file_residuals Output CSV filename, where ppm errors of calibrants before and after model fitting parameters are written to (pass empty string to skip)
      @param file_residuals_plot Output PNG image of the ppm errors of calibrants (pass empty string to skip)
      @param rscript_executable Full path to the Rscript executable
      @return true upon successful calibration

    */
    bool calibrate(PeakMap& exp, 
                   const IntList& target_mslvl,
                   MZTrafoModel::MODELTYPE model_type,
                   double rt_chunk,
                   bool use_RANSAC,
                   double post_ppm_median,
                   double post_ppm_MAD,
                   const String& file_models = "",
                   const String& file_models_plot = "",
                   const String& file_residuals = "",
                   const String& file_residuals_plot = "",
                   const String& rscript_executable = "Rscript");

    /**
      @brief Transform a precursor's m/z

      Calibrate m/z of precursors.

      @param pcs Uncalibrated Precursors
      @param trafo The calibration function to apply
    */
    static void applyTransformation(std::vector<Precursor>& pcs, const MZTrafoModel& trafo);

    /**
      @brief Transform a spectrum (data+precursor)

      See applyTransformation(MSExperiment, ...) for details.

      @param spec Uncalibrated MSSpectrum
      @param target_mslvl List (can be unsorted) of MS levels to calibrate
      @param trafo The calibration function to apply
    */
    static void applyTransformation(PeakMap::SpectrumType& spec, const IntList& target_mslvl, const MZTrafoModel& trafo);

    /**
      @brief Transform spectra from a whole map (data+precursor)

      All data peaks and precursor information (if present) are calibrated in m/z.

      Only spectra whose MS-level is contained in 'target_mslvl' are calibrated.
      In addition, if a fragmentation spectrum's precursor information originates from an MS level in 'target_mslvl',
      the precursor (not the spectrum itself) is also subjected to calibration.
      E.g., If we only have MS and MS/MS spectra: for 'target_mslvl' = {1} then all MS1 spectra and MS2 precursors are calibrated.
      If 'target_mslvl' = {2}, only MS2 spectra (not their precursors) are calibrated.
      If 'target_mslvl' = {1,2} all spectra and precursors are calibrated.
            

      @param exp Uncalibrated peak map
      @param target_mslvl List (can be unsorted) of MS levels to calibrate
      @param trafo The calibration function to apply
    */
    static void applyTransformation(PeakMap& exp, const IntList& target_mslvl, const MZTrafoModel& trafo);
  
  protected:

    /// statistics when adding peptide calibrants
    struct CalibrantStats_
    {
      CalibrantStats_(const double tol_ppm)
        : tol_ppm_(tol_ppm)
      {};
      Size cnt_empty = 0; ///< cases of empty PepIDs (no hits)
      Size cnt_nomz = 0; ///< cases of no m/z value
      Size cnt_nort = 0; ///< cases of no RT value
      Size cnt_decal = 0; ///< cases of large gap (>tol_ppm) between theoretical peptide weight (from sequence) and precursor mass
      Size cnt_total = 0; ///< total number of cases

      void print() const
      {
        if (cnt_empty > 0) OPENMS_LOG_WARN << "Warning: " << cnt_empty << "/" << cnt_total << " calibrations points were skipped, since they have no peptide sequence!" << std::endl;
        if (cnt_nomz > 0) OPENMS_LOG_WARN << "Warning: " << cnt_nomz << "/" << cnt_total << " calibrations points were skipped, since they have no m/z value!" << std::endl;
        if (cnt_nort > 0) OPENMS_LOG_WARN << "Warning: " << cnt_nort << "/" << cnt_total << " calibrations points were skipped, since they have no RT value!" << std::endl;
        if (cnt_decal > 0) OPENMS_LOG_WARN << "Warning: " << cnt_decal << "/" << cnt_total << " calibrations points were skipped, since their theoretical weight is more than " << tol_ppm_ << " ppm away from their measured mass!" << std::endl;
      }

      private:
      const double tol_ppm_; ///< tolerance used for counting cnt_decal
    };

    /**
      @brief Add(no prior clear) calibrants to internal list.
      
      Extracts only the first hit from each peptide identification.
      Hits are sorted beforehand.
      Ambiguities should be resolved before, e.g. using IDFilter.

      A filtering step is done in the m/z dimension using @p tol_ppm.
      Since precursor masses could be annotated wrongly (e.g. isotope peak instead of mono),
      larger outliers are removed before accepting an ID as calibrant.

      @param pep_id A single PeptideID (e.g. from an idXML file); only the top peptide hit is used
      @param tol_ppm Only accept ID's whose theoretical mass deviates at most this much from annotated
      @param stats Update stats, if calibrant cannot be used (no RT, no MZ, no sequence, out-of tolerance)

    */
    void fillID_( const PeptideIdentification& pep_id, const double tol_ppm, CalibrantStats_& stats);

    /// calls fillID_ on all PeptideIDs
    void fillIDs_(const std::vector<PeptideIdentification>& pep_ids, const double tol_ppm, CalibrantStats_& stats);

    /// determine if sequence is within tol_ppm and update stats; fills mz_ref with the theoretical m/z of the sequence
    bool isDecalibrated_(const PeptideIdentification& pep_id, const double mz_obs, const double tol_ppm, CalibrantStats_& stats, double& mz_ref);

    /**
     @brief Calibrate m/z of a spectrum, ignoring precursors!

     This method is not exposed as public, because its easy to be misused on spectra while forgetting about the precursors of high-level spectra.
    */
    static void applyTransformation_(PeakMap::SpectrumType& spec, const MZTrafoModel& trafo);

  private:
    CalibrationData cal_data_;
  }; // class InternalCalibration
  
} // namespace OpenMS

