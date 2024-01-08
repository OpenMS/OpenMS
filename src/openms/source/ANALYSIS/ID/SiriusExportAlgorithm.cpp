// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka, Axel Walter $
// $Authors: Oliver Alka, Lukas Zimmermann $
// --------------------------------------------------------------------------

#include <boost/foreach.hpp> // must be first, otherwise Q_FOREACH macro will wreak havoc

#include <OpenMS/ANALYSIS/ID/SiriusExportAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/SiriusMSConverter.h>

#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/SpectrumLookup.h>
#include <OpenMS/SYSTEM/File.h>

#include <QtCore/QDir>
#include <QtCore/QDirIterator>
#include <QtCore/QString>
#include <QtCore/QProcess>
#include <sstream>
#include <boost/regex.hpp>

namespace OpenMS
{
  // ###################
  // Set subtool parameters
  // ###################

    using OpenMSName = String;
    using DefaultValue = ParamValue;
    using Description = String;

    SiriusExportAlgorithm::SiriusExportAlgorithm() :
      DefaultParamHandler("SiriusExportAlgorithm"),
      preprocessing(Preprocessing(this))

    {
      // Defines the Parameters for preprocessing and SIRIUS subtools
      preprocessing.parameters();

      defaultsToParam_();
    }

    void SiriusExportAlgorithm::Preprocessing::parameters()
    {
      parameter(
                  OpenMSName("filter_by_num_masstraces"),
                  DefaultValue(1),
                  Description("Number of mass traces each feature has to have to be included. "
                              "To use this parameter, setting the feature_only flag is necessary")
                ).withMinInt(1);

      parameter(
                  OpenMSName("precursor_mz_tolerance"),
                  DefaultValue(10.0),
                  Description("Tolerance window for precursor selection (Feature selection in regard to the precursor)")
                );

      parameter(
                  OpenMSName("precursor_mz_tolerance_unit"),
                  DefaultValue("ppm"),
                  Description("Unit of the precursor_mz_tolerance")
               ).withValidStrings({"Da", "ppm"});

      parameter(
                  OpenMSName("precursor_rt_tolerance"),
                  DefaultValue(5.0),
                  Description("Tolerance window (left and right) for precursor selection [seconds]")
               );

      parameter(
                  OpenMSName("isotope_pattern_iterations"),
                  DefaultValue(3),
                  Description("Number of iterations that should be performed to extract the C13 isotope pattern. "
                              "If no peak is found (C13 distance) the function will abort. "
                              "Be careful with noisy data - since this can lead to wrong isotope patterns")
                );

      flag(
            OpenMSName("feature_only"),
            Description("Uses the feature information from in_featureinfo to reduce the search space to MS2 "
                        "associated with a feature")
          );

      flag(
            OpenMSName("no_masstrace_info_isotope_pattern"),
            Description("Use this flag if the masstrace information from a feature should be discarded "
                       "and the isotope_pattern_iterations should be used instead")
          );
    }


    void SiriusExportAlgorithm::updateExistingParameter(const OpenMS::Param &param)
    {
      for (auto it = param.begin(); it != param.end(); ++it)
      {
        const std::string name = it.getName();
        if (hasFullNameParameter(name))
        {
          vector<std::string> tags(it->tags.begin(), it->tags.end());
          param_.setValue(name, it->value, it->description, tags);
        }
      }
    }

    bool SiriusExportAlgorithm::hasFullNameParameter(const OpenMS::String &name) const
    {
      //alternative: return std::any_of(param_.begin(), param_.end(), [&](const auto& it) { return it.name == name; });
      for (auto it = param_.begin(); it != param_.end(); ++it)
      {
        if (it.getName() == name)
        {
          return true;
        }
      }
      return false;
    }


    // ################
    // Algorithm
    // ################
    void SiriusExportAlgorithm::preprocessingSirius(const String& featureinfo,
                                                     const MSExperiment& spectra,
                                                     FeatureMapping::FeatureMappingInfo& fm_info,
                                                     FeatureMapping::FeatureToMs2Indices& feature_mapping) const
    {
      // if fileparameter is given and should be not empty
      if (!featureinfo.empty())
      {
        if (File::exists(featureinfo) && !File::empty(featureinfo))
        {
          // read featureXML          
          FeatureMap feature_map;
          FileHandler().loadFeatures(featureinfo, feature_map);

          UInt num_masstrace_filter = getFilterByNumMassTraces();
          double precursor_mz_tol = getPrecursorMzTolerance();
          double precursor_rt_tol = getPrecursorRtTolerance();

          if (num_masstrace_filter != 1 && !isFeatureOnly())
          {
            num_masstrace_filter = 1;
            OPENMS_LOG_WARN << "Parameter: filter_by_num_masstraces, was set to 1 to retain the adduct information for all MS2 spectra, if available. Masstrace filtering only makes sense in combination with feature_only." << std::endl;
          }

          // filter feature by number of masstraces
          auto map_it = remove_if(feature_map.begin(), feature_map.end(),
                                  [&num_masstrace_filter](const Feature &feat) -> bool
                                  {
                                    unsigned int n_masstraces = feat.getMetaValue(Constants::UserParam::NUM_OF_MASSTRACES);
                                    return n_masstraces < num_masstrace_filter;
                                  });
          feature_map.erase(map_it, feature_map.end());
  
          fm_info.feature_maps.push_back(feature_map);
          fm_info.kd_tree.addMaps(fm_info.feature_maps); // KDTree references into feature_map
  
          // mapping of MS2 spectra to features
          feature_mapping = FeatureMapping::assignMS2IndexToFeature(spectra,
                                                                    fm_info,
                                                                    precursor_mz_tol,
                                                                    precursor_rt_tol,
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

    void SiriusExportAlgorithm::logFeatureSpectraNumber(const String& featureinfo,
                                                         const FeatureMapping::FeatureToMs2Indices& feature_mapping,
                                                         const MSExperiment& spectra) const
    {
      // number of features to be processed
      if (isFeatureOnly() && !featureinfo.empty())
      {
        OPENMS_LOG_INFO << "Number of features to be processed: " << feature_mapping.assignedMS2.size() << std::endl;
      }
      else if (!featureinfo.empty())
      {
        OPENMS_LOG_INFO << "Number of features to be processed: " << feature_mapping.assignedMS2.size() << std::endl;
        OPENMS_LOG_INFO << "Number of additional MS2 spectra to be processed: " << feature_mapping.unassignedMS2.size() << std::endl;
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
      for (size_t i = 0; i < mzML_files.size(); ++i) {
          // load experiment
          MSExperiment spectra;
          FileHandler().loadExperiment(mzML_files[i], spectra, {FileTypes::MZML});

          // run masstrace filter and feature mapping
          FeatureMapping::FeatureMappingInfo fm_info;
          FeatureMapping::FeatureToMs2Indices feature_mapping;

          // check if 'featureXML_files' is empty and pass an empty string if it is
          String feature_info_to_pass = featureXML_files.empty() ? "" : featureXML_files[i];
          SiriusExportAlgorithm::preprocessingSirius(feature_info_to_pass,
                                        spectra,
                                        fm_info,
                                        feature_mapping);

          // returns Log of feature and/or spectra number
          SiriusExportAlgorithm::logFeatureSpectraNumber(feature_info_to_pass, feature_mapping, spectra);

          // temporary vector to store compound information for the current file
          std::vector<SiriusMSFile::CompoundInfo> temp_cmpinfo;
          SiriusMSFile::store(spectra,
                              os,
                              feature_mapping,
                              SiriusExportAlgorithm::isFeatureOnly(),
                              SiriusExportAlgorithm::getIsotopePatternIterations(),
                              SiriusExportAlgorithm::isNoMasstraceInfoIsotopePattern(),
                              temp_cmpinfo,
                              i);

          // Append the compound information of the current file to the overall vector
          v_cmpinfo.insert(v_cmpinfo.end(), temp_cmpinfo.begin(), temp_cmpinfo.end());
      }

      os.close();

      if (!out_compoundinfo.empty()) 
      {
        SiriusMSFile::saveCompoundInfoAsTSV(v_cmpinfo, out_compoundinfo, SiriusExportAlgorithm::isFeatureOnly());
      }
    }

    struct OPENMS_DLLAPI SiriusWorkspaceIndex
    {
      SiriusWorkspaceIndex(int array_index, int scan_index) : array_index {array_index}, scan_index {scan_index} {}
      int array_index, scan_index;
    };
    void  SiriusExportAlgorithm::sortSiriusWorkspacePathsByScanIndex(std::vector<String>& subdirs)
    {
      std::vector<String> sorted_subdirs;
      std::vector<SiriusWorkspaceIndex> indices;

      boost::regex regexp(R"(--(?<SCAN>\d+)--)");
      for (size_t i = 0; i < subdirs.size(); i++)
      {
        indices.emplace_back(i, SpectrumLookup::extractScanNumber(subdirs[i], regexp, false));
      }

      std::sort(indices.begin(),
                indices.end(),
                [](const SiriusWorkspaceIndex& i, const SiriusWorkspaceIndex& j) { return i.scan_index < j.scan_index; } );

      sorted_subdirs.reserve(indices.size());
      for (const auto& index : indices)
      {
        sorted_subdirs.emplace_back(std::move(subdirs[index.array_index]));
      }

      sorted_subdirs.swap(subdirs);
    }

    // ################
    // Parameter handling
    // ################
    SiriusExportAlgorithm::ParameterModifier SiriusExportAlgorithm::ParameterSection::parameter(
            const String &parameter_name,
            const ParamValue &default_value,
            const String &parameter_description)
    {
      const String full_parameter = toFullParameter(parameter_name);
      openms_to_sirius[full_parameter] = parameter_name;
      enclose->defaults_.setValue(full_parameter, default_value, parameter_description);
      return ParameterModifier(full_parameter, enclose);
    }

  void SiriusExportAlgorithm::ParameterSection::flag(
            const OpenMS::String &parameter_name,
            const OpenMS::String &parameter_description)
    {
      parameter(parameter_name, DefaultValue("false"), parameter_description)
        .withValidStrings({"true", "false"});
    }
} // namespace OpenMS

/// @endcond
