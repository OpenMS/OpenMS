// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <map>

namespace OpenMS
{
    /**
        @brief Common functions for DDA workflows
    
        @ingroup Analysis_ID
    */
    class OPENMS_DLLAPI DDAWorkflowCommons
    {
        public:
        /* @brief create Map between mzML file and corresponding id file
         * Checks implemented:
         * 1. Check if the number of spectra and id files match.
         * 2. If spectra and id files share common base names (without extension) 
         *    but appear in different order, throw an error.
        */
        static std::map<String, String> mapMzML2Ids(StringList & in, StringList & in_ids)
        {
            // Detect the common case that ID files have same names as spectra files
            if (!File::validateMatchingFileNames(in, in_ids, true, true, false)) // only basenames, without extension, only order
            {
            // Spectra and id files have the same set of basenames but appear in different order. -> this is most likely an error
            throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                "ID and spectra file match but order of file names seem to differ. They need to be provided in the same order.");
            }

            std::map<String, String> mzfile2idfile;
            for (Size i = 0; i != in.size(); ++i)
            {
                const String& in_abs_path = File::absolutePath(in[i]);
                const String& id_abs_path = File::absolutePath(in_ids[i]);
                mzfile2idfile[in_abs_path] = id_abs_path;      
                OPENMS_LOG_DEBUG << ("Spectra: " + in[i] + "\t Ids: " + in_ids[i],  1) << std::endl;
            }
            return mzfile2idfile;
        }
    }

/**
 * @brief Small helper to get the mapping from id files to mzML files
 *
 * Basically just reverses the mapMzML2Ids function. 
 * Potential improvement: Could be combined into a single functionexposed to the user.
 */
  std::map<String, String> mapId2MzMLs(const map<String, String>& m2i)
  {
    map<String, String> idfile2mzfile;
    for (const auto& m : m2i)
    {
      idfile2mzfile[m.second] = m.first;
    }
    return idfile2mzfile;
  }


/**
 * Estimates the median chromatographic full width at half maximum (FWHM) for a given MSExperiment.
 *
 * @param ms_centroided The centroided MSExperiment for which to estimate the FWHM.
 * @return The estimated median chromatographic FWHM based on the top 1000 intensity mass traces.
 */
  double estimateMedianChromatographicFWHM(MSExperiment & ms_centroided)
  {
    MassTraceDetection mt_ext;
    Param mtd_param = mt_ext.getParameters();

    OPENMS_LOG_DEBUG << "Parameters passed to MassTraceDetection" << mtd_param << std::endl;

    std::vector<MassTrace> m_traces;
    mt_ext.run(ms_centroided, m_traces, 1000);

    std::vector<double> fwhm_1000;
    for (auto &m : m_traces)
    {
      if (m.getSize() == 0) continue;
      m.updateMeanMZ();
      m.updateWeightedMZsd();
      double fwhm = m.estimateFWHM(false);
      fwhm_1000.push_back(fwhm);
    }

    double median_fwhm = Math::median(fwhm_1000.begin(), fwhm_1000.end());

    return median_fwhm;
  }

/**
 * @brief Recalibrates the masses of the MSExperiment using peptide identifications.
 *
 * This function recalibrates the masses of the MSExperiment by applying a mass recalibration
 * based on the theoretical masses from identification data.
 *
 * @param ms_centroided The MSExperiment object containing the centroided spectra.
 * @param peptide_ids The vector of PeptideIdentification objects containing the peptide identifications.
 * @param id_file_abs_path The absolute path of the identification file.
 */
  void recalibrateMS1(MSExperiment & ms_centroided,
    vector<PeptideIdentification>& peptide_ids,
    const String & id_file_abs_path = "")
  {
    InternalCalibration ic;
    ic.setLogType(log_type_);
    ic.fillCalibrants(peptide_ids, 25.0); // >25 ppm maximum deviation defines an outlier TODO: check if we need to adapt this
    if (ic.getCalibrationPoints().size() <= 1) return;

    // choose calibration model based on number of calibration points

    // there seem to be some problems with the QUADRATIC model that we first need to investigate
    //MZTrafoModel::MODELTYPE md = (ic.getCalibrationPoints().size() == 2) ? MZTrafoModel::LINEAR : MZTrafoModel::QUADRATIC;
    //bool use_RANSAC = (md == MZTrafoModel::LINEAR || md == MZTrafoModel::QUADRATIC);

    MZTrafoModel::MODELTYPE md = MZTrafoModel::LINEAR;
    bool use_RANSAC = true;

    Size RANSAC_initial_points = (md == MZTrafoModel::LINEAR) ? 2 : 3;
    Math::RANSACParam p(RANSAC_initial_points, 70, 10, 30, true); // TODO: check defaults (taken from tool)
    MZTrafoModel::setRANSACParams(p);
    // these limits are a little loose, but should prevent grossly wrong models without burdening the user with yet another parameter.
    MZTrafoModel::setCoefficientLimits(25.0, 25.0, 0.5); 

    IntList ms_level = {1};
    double rt_chunk = 300.0; // 5 minutes
    String qc_residual_path, qc_residual_png_path;
    if (!id_file_abs_path.empty())
    {
      const String & id_basename = File::basename(id_file_abs_path);
      qc_residual_path = id_basename + "qc_residuals.tsv";
      qc_residual_png_path = id_basename + "qc_residuals.png";
    } 

    if (!ic.calibrate(ms_centroided, 
                  ms_level, md, rt_chunk, use_RANSAC, 
                  10.0,
                  5.0, 
                  "",                      
                  "",
                  qc_residual_path,
                  qc_residual_png_path,
                  "Rscript"))
    {
      OPENMS_LOG_WARN << "\nCalibration failed. See error message above!" << std::endl;
    }
  }


/**
 * @brief Extracts seeding features from centroided MS data (e.g., for untarged extraction).
 *
 * MS1 spectra are subjected to a threshold filter to removelow-intensity peaks, 
 * and then uses the FeatureFinderMultiplex algorithm to identify potential seeding features.
 * The function also takes into account the median full width at half maximum (FWHM) of the peaks 
 * to adjust the FeatureFinderMultiplex parameters for better seed detection.
 *
 * @param[in] ms_centroided The MSExperiment object containing centroided mass spectrometry data. Only MS1 level
 *                          spectra are considered for seed calculation.
 * @param[out] seeds The FeatureMap object where the identified seeding features will be stored.
 * @param[in] median_fwhm The median FWHM of the peaks, used to adjust the FeatureFinderMultiplex parameters for
 *                        seed detection.
 *
 * @note The function currently uses hardcoded parameters for the threshold filter and FeatureFinderMultiplex
 *       algorithm, which may need to be derived from the data or provided as function arguments in future
 *       implementations.
 */
  void calculateSeeds(
    const MSExperiment & ms_centroided, 
    const double intensity_threshold,
    FeatureMap & seeds, 
    double median_fwhm
    Size charge_min = 2,
    Size charge_max = 5)
  {
    //TODO: Actually FFM provides a parameter for minimum intensity. Also it copies the full experiment again once or twice.
    MSExperiment e;
    for (const auto& s : ms_centroided)
    { 
      if (s.getMSLevel() == 1) 
      {              
        e.addSpectrum(s);
      }
    }

    ThresholdMower threshold_mower_filter;
    Param tm = threshold_mower_filter.getParameters();
    tm.setValue("threshold", intensity_threshold);;  // TODO: derive from data
    threshold_mower_filter.setParameters(tm);
    threshold_mower_filter.filterPeakMap(e);

    FeatureFinderMultiplexAlgorithm algorithm;
    Param p = algorithm.getParameters();
    p.setValue("algorithm:labels", ""); // unlabeled only
    p.setValue("algorithm:charge", String(charge_min) + ":" + String(charge_max));
    p.setValue("algorithm:rt_typical", median_fwhm * 3.0);
    p.setValue("algorithm:rt_band", 3.0); // max 3 seconds shifts between isotopic traces (not sure if needed)
    p.setValue("algorithm:rt_min", median_fwhm * 0.5);
    p.setValue("algorithm:spectrum_type", "centroid");
    algorithm.setParameters(p);
    //FIXME progress of FFM is not printed at all
    const bool progress(true);
    algorithm.run(e, progress);
    seeds = algorithm.getFeatureMap(); 
    OPENMS_LOG_INFO << "Using " << seeds.size() << " seeds from untargeted feature extraction." << endl;
  }

}

