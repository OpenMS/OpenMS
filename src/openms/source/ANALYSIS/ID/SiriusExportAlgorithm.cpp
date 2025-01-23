// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka, Axel Walter $
// $Authors: Oliver Alka, Lukas Zimmermann $
// --------------------------------------------------------------------------

#include <boost/foreach.hpp> // must be first, otherwise Q_FOREACH macro will wreak havoc

#include <OpenMS/ANALYSIS/ID/SiriusExportAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/SiriusMSConverter.h>

#include <OpenMS/FORMAT/FileHandler.h>

namespace OpenMS
{
  SiriusExportAlgorithm::SiriusExportAlgorithm() :
    DefaultParamHandler("SiriusExportAlgorithm")
  {
    defaults_.setValue("filter_by_num_masstraces", 1, "Number of mass traces each feature has to have to be included. To use this parameter, setting the feature_only flag is necessary");
    defaults_.setMinInt("filter_by_num_masstraces", 1);

    defaults_.setValue("precursor_mz_tolerance", 10.0, "Tolerance window for precursor selection (Feature selection in regard to the precursor)");
    
    defaults_.setValue("precursor_mz_tolerance_unit", "ppm", "Unit of the preprocessing_precursor_mz_tolerance");
    defaults_.setValidStrings("precursor_mz_tolerance_unit", {"ppm","Da"});
    
    defaults_.setValue("precursor_rt_tolerance", 5.0, "Tolerance window (left and right) for precursor selection [seconds]");

    defaults_.setValue("isotope_pattern_iterations", 3, "Number of iterations that should be performed to extract the C13 isotope pattern. If no peak is found (C13 distance) the function will abort. Be careful with noisy data - since this can lead to wrong isotope patterns");

    defaults_.setValue("feature_only", "false", "Uses the feature information from in_featureinfo to reduce the search space to MS2 associated with a feature");
    defaults_.setValidStrings("feature_only", {"false","true"});

    defaults_.setValue("no_masstrace_info_isotope_pattern", "false", "Set to true if the masstrace information from a feature should be discarded and the isotope_pattern_iterations should be used instead");
    defaults_.setValidStrings("no_masstrace_info_isotope_pattern", {"false","true"});
    defaultsToParam_();
  }

  // ################
  // Algorithm
  // ################
  void SiriusExportAlgorithm::preprocessing(const String& featureXML_path,
                                            const MSExperiment& spectra,
                                            FeatureMapping::FeatureMappingInfo& feature_mapping_info,
                                            FeatureMapping::FeatureToMs2Indices& feature_ms2_indices) const
{
    // if fileparameter is given and should be not empty
    if (!featureXML_path.empty())
    {
      Size preprocessing_filter_by_num_masstraces = getFilterByNumMassTraces();
      if (File::exists(featureXML_path) && !File::empty(featureXML_path))
      {
        // read featureXML          
        FeatureMap feature_map;
        FileHandler().loadFeatures(featureXML_path, feature_map);

        if (preprocessing_filter_by_num_masstraces != 1 && !isFeatureOnly())
        {
          preprocessing_filter_by_num_masstraces = 1;
          OPENMS_LOG_WARN << "Parameter: preprocessing_filter_by_num_masstraces, was set to 1 to retain the adduct information for all MS2 spectra, if available. Masstrace filtering only makes sense in combination with feature_only." << std::endl;
        }

        // filter feature by number of masstraces
        auto map_it = remove_if(feature_map.begin(), feature_map.end(),
                                [&preprocessing_filter_by_num_masstraces](const Feature &feat) -> bool
                                {
                                  unsigned int n_masstraces = feat.getMetaValue(Constants::UserParam::NUM_OF_MASSTRACES);
                                  return n_masstraces < preprocessing_filter_by_num_masstraces;
                                });
        feature_map.erase(map_it, feature_map.end());

        feature_mapping_info.feature_maps.push_back(feature_map);
        feature_mapping_info.kd_tree.addMaps(feature_mapping_info.feature_maps); // KDTree references into feature_map

        // mapping of MS2 spectra to features
        feature_ms2_indices = FeatureMapping::assignMS2IndexToFeature(spectra,
                                                                  feature_mapping_info,
                                                                  getPrecursorMzTolerance(),
                                                                  getPrecursorRtTolerance(),
                                                                  precursorMzToleranceUnitIsPPM());
      }
      else
      {
        throw OpenMS::Exception::FileEmpty(__FILE__,
                                            __LINE__,
                                            __FUNCTION__,
                                            "Error: FeatureXML was empty, please provide a valid file.");
      }
    }
  }


  void SiriusExportAlgorithm::logFeatureSpectraNumber(const String& featureXML_path,
                                                        const FeatureMapping::FeatureToMs2Indices& feature_ms2_indices,
                                                        const MSExperiment& spectra) const
  {
    // number of features to be processed
    if (isFeatureOnly() && !featureXML_path.empty())
    {
      OPENMS_LOG_INFO << "Number of features to be processed: " << feature_ms2_indices.assignedMS2.size() << std::endl;
    }
    else if (!featureXML_path.empty())
    {
      OPENMS_LOG_INFO << "Number of features to be processed: " << feature_ms2_indices.assignedMS2.size() << std::endl;
      OPENMS_LOG_INFO << "Number of additional MS2 spectra to be processed: " << feature_ms2_indices.unassignedMS2.size() << std::endl;
    } 
    else
    {
      long count_ms2 = count_if(spectra.begin(), spectra.end(),
              [](const MSSpectrum &spectrum) { return spectrum.getMSLevel() == 2; });

      OPENMS_LOG_INFO << "Number of MS2 spectra to be processed: " << count_ms2 << std::endl;
    }
  }

  void SiriusExportAlgorithm::run(const StringList& mzML_files,
                                  const StringList& featureXML_files,
                                  const String& out_ms,
                                  const String& out_compoundinfo) const
  {
    // loop over all spectra in all files and write data to ofstream
    ofstream os;

    // create temporary input file (.ms)
    os.open(out_ms);
    if (!os)
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, out_ms);
    }
    os.precision(12);

    std::vector<SiriusMSFile::CompoundInfo> v_cmpinfo; // To store compound information for all files
    for (size_t i = 0; i < mzML_files.size(); ++i) 
    {
      // load experiment
      MSExperiment spectra;
      FileHandler().loadExperiment(mzML_files[i], spectra, {FileTypes::MZML});

      // run masstrace filter and feature mapping
      FeatureMapping::FeatureMappingInfo feature_mapping_info;
      FeatureMapping::FeatureToMs2Indices feature_ms2_indices;

      // check if 'featureXML_files' is empty and pass an empty string if it is
      String feature_info_to_pass = featureXML_files.empty() ? "" : featureXML_files[i];
      SiriusExportAlgorithm::preprocessing(feature_info_to_pass,
                                    spectra,
                                    feature_mapping_info,
                                    feature_ms2_indices);

      // returns Log of feature and/or spectra number
      SiriusExportAlgorithm::logFeatureSpectraNumber(feature_info_to_pass, feature_ms2_indices, spectra);

      // temporary vector to store compound information for the current file
      std::vector<SiriusMSFile::CompoundInfo> temp_cmpinfo;
      SiriusMSFile::store(spectra,
                          os,
                          feature_ms2_indices,
                          isFeatureOnly(),
                          getIsotopePatternIterations(),
                          isNoMasstraceInfoIsotopePattern(),
                          temp_cmpinfo,
                          i);
      // Append the compound information of the current file to the overall vector
      v_cmpinfo.insert(v_cmpinfo.end(), temp_cmpinfo.begin(), temp_cmpinfo.end());
    }

    os.close();

    if (!out_compoundinfo.empty()) 
    {
      SiriusMSFile::saveFeatureCompoundInfoAsTSV(v_cmpinfo, out_compoundinfo);
    }
  }
} // namespace OpenMS

/// @endcond
