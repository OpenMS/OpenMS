// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/QUANTITATION/DDAWorkflowCommons.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FEATUREFINDER/MassTraceDetection.h>
#include <OpenMS/PROCESSING/FILTERING/ThresholdMower.h>
#include <OpenMS/PROCESSING/CALIBRATION/InternalCalibration.h>
#include <OpenMS/PROCESSING/CALIBRATION/MZTrafoModel.h>
#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/FEATUREFINDER/FeatureFinderMultiplexAlgorithm.h>

#include <map>
#include <vector>

using namespace std;
namespace OpenMS
{
    std::map<String, String> DDAWorkflowCommons::mapId2MzMLs(const std::map<String, String>& m2i)
    {
        std::map<String, String> idfile2mzfile;
        for (const auto& m : m2i)
        {
            idfile2mzfile[m.second] = m.first;
        }
        return idfile2mzfile;
    }


    std::map<String, String> DDAWorkflowCommons::mapMzML2Ids(StringList & in, StringList & in_ids)
    {
        // validate file lists (use only basename and ignore extension)
        auto validation_result = File::validateMatchingFileNames(in, in_ids, true, true);
        // we try to fail early (without parsing files) if the input is obviously wrong
        // check for two major mistakes:
        //  1. different number of files (-> certainly wrong)
        //  2. same number of files but different order (-> certainly wrong)
        // If some files differ in names, we can't be sure at this point and skip this test for now.
        // We need to look into the ID files to infer the spectra filenames later to be sure a mistake was made.
        switch (validation_result)
        {
        case File::MatchingFileListsStatus::SET_MISMATCH:
            if (in.size() != in_ids.size())
            {
            OPENMS_LOG_FATAL_ERROR << "ID and spectra file lists differ in size. Please provide the same number of files for spectra and ID." << endl;
            throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                "ID and spectra file lists differ in size. Please provide the same number of files for spectra and ID.");          
            }
            else
            { // same number of files but filenames differ (we will try to read the spectra filenames from the id files later)
            OPENMS_LOG_DEBUG << "ID and spectra file lists differ. Please provide the same files in the same order." << std::endl;
            OPENMS_LOG_DEBUG << "File in spectra file list: " << std::endl;
            for (const auto& f : in)
            {
                OPENMS_LOG_DEBUG << f << std::endl;
            }
            OPENMS_LOG_DEBUG << "File in ID file list: " << std::endl;
            for (const auto& f : in_ids)
            {
                OPENMS_LOG_DEBUG << f << std::endl;
            }
            OPENMS_LOG_DEBUG << "Will try to infer spectra filenames from id files later." << std::endl;
            }
            break;
        case File::MatchingFileListsStatus::ORDER_MISMATCH:
            OPENMS_LOG_DEBUG << "ID and spectra file match but order of file names seem to differ. Please provide the same files in the same order." << std::endl;
            OPENMS_LOG_DEBUG << "File in spectra file list: " << std::endl;
            for (const auto& f : in)
            {
                OPENMS_LOG_DEBUG << f << std::endl;
            }
            OPENMS_LOG_WARN << "File in ID file list: " << endl;
            for (const auto& f : in_ids)
            {
                OPENMS_LOG_DEBUG << f << std::endl;
            }
            throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "ID and spectra file match but order of file names seem to differ. They need to be provided in the same order.");
            break;
        case File::MatchingFileListsStatus::MATCH:
            OPENMS_LOG_INFO << "Info: ID files have the same names as spectra files." << std::endl;
            break;      
        }

        map<String, String> mzfile2idfile;
        for (Size i = 0; i != in.size(); ++i)
        {
            const String& in_abs_path = File::absolutePath(in[i]);
            const String& id_abs_path = File::absolutePath(in_ids[i]);
            mzfile2idfile[in_abs_path] = id_abs_path;      
            OPENMS_LOG_DEBUG << "Spectra: " << in[i] << "\t Ids: " << in_ids[i] << std::endl;
        }
        return mzfile2idfile;
    }


    double DDAWorkflowCommons::estimateMedianChromatographicFWHM(MSExperiment & ms_centroided)
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
        

    void DDAWorkflowCommons::recalibrateMS1(MSExperiment & ms_centroided,
        PeptideIdentificationList& peptide_ids,
        const String & id_file_abs_path )
    {
        InternalCalibration ic;
        // ic.setLogType(log_type_);
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


    void DDAWorkflowCommons::calculateSeeds(
        const MSExperiment & ms_centroided, 
        const double intensity_threshold,
        FeatureMap & seeds, 
        double median_fwhm,
        Size charge_min,
        Size charge_max
    )
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
        OPENMS_LOG_INFO << "Using " << String(seeds.size()) << " seeds from untargeted feature extraction." << std::endl;
    }
}
