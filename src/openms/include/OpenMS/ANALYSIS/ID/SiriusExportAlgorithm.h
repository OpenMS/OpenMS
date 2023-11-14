// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka, Axel Walter $
// $Authors: Oliver Alka, Lukas Zimmermann $
// --------------------------------------------------------------------------

#pragma once 

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureMapping.h>

#include <OpenMS/SYSTEM/File.h>

#include <unordered_map>
#include <QtCore/QString>
#include <QtCore/QStringList>

using namespace std;

namespace OpenMS
{
  class FeatureMap;
  class File;
  class KDTreeFeatureMaps;

  class OPENMS_DLLAPI SiriusExportAlgorithm : public DefaultParamHandler
    {
    public:
      /// default constructor
      SiriusExportAlgorithm();

      /*
       * Accessors for Preprocessing Parameters
       */

      bool isFeatureOnly() const { return preprocessing.getValue("feature_only").toBool(); }
      UInt getFilterByNumMassTraces() const { return preprocessing.getValue("filter_by_num_masstraces"); }
      double getPrecursorMzTolerance() const { return preprocessing.getValue("precursor_mz_tolerance"); }
      double getPrecursorRtTolerance() const { return preprocessing.getValue("precursor_rt_tolerance"); }
      bool precursorMzToleranceUnitIsPPM() const { return preprocessing.getValue("precursor_mz_tolerance_unit") == "ppm"; }
      bool isNoMasstraceInfoIsotopePattern() const { return preprocessing.getValue("no_masstrace_info_isotope_pattern").toBool(); }
      int getIsotopePatternIterations() const { return  preprocessing.getValue("isotope_pattern_iterations"); }

      /**
       * @brief Updates all parameters that already exist in this DefaultParamHandler
       * with the values provided by the input param object.
       *
       * @param param The Param object supplying updated parameter values. Keys that exist in the param
       *        parameter, but not in this DefaultParamHander, are ignored.
       */
      void updateExistingParameter(const Param &param);

      /**
       * @brief Checks whether this DefaultParamHandler has a ParamEntry with the provided name.
       * @param name The name of the ParamEntry that should be checked for its existence in this DefaultParamHandler
       * @return Whether this DefaultParamHandler has an ParamEntry for the provided name.
       */
      bool hasFullNameParameter(const String &name) const;

      /**
      @brief Preprocessing needed for SIRIUS

      @return FeatureToMS2Indices
            
      Filter number of masstraces and perform feature mapping.

      @param featureinfo Path to featureXML
      @param spectra Input of MSExperiment with spectra information
      @param fm_info Emtpy - stores FeatureMaps and KDTreeMaps internally 
      @param feature_mapping Empty FeatureToMs2Indices
      */
      void preprocessingSirius(const String& featureinfo,
                               const MSExperiment& spectra,
                               FeatureMapping::FeatureMappingInfo& fm_info,
                               FeatureMapping::FeatureToMs2Indices& feature_mapping) const;

      /**
      @brief logs number of features and spectra used

      Prints the number of features and spectra used (OPENMS_LOG_INFO)

      @param featureinfo Path to featureXML
      @param feature_mapping FeatureToMs2Indices with feature mapping
      @param spectra Input of MSExperiment with spectra information
      */
      void logFeatureSpectraNumber(const String& featureinfo,
                                   const FeatureMapping::FeatureToMs2Indices& feature_mapping,
                                   const MSExperiment& spectra) const;

   private:
    class ParameterModifier
    {
      const String openms_param_name;
      SiriusExportAlgorithm *enclose;

    public:
      explicit ParameterModifier(const String &param_name, SiriusExportAlgorithm *enclose) :
              openms_param_name(param_name), enclose(enclose) {}

      void withValidStrings(initializer_list<std::string> choices)
      {
        enclose->defaults_.setValidStrings(openms_param_name, choices);
      }

      void withMinInt(int value)
      {
        enclose->defaults_.setMinInt(openms_param_name, value);
      }
    };

    class ParameterSection
    {
      // Maps the OpenMS Parameter Names to the one for Sirius
      unordered_map<String, String> openms_to_sirius;

      String toFullParameter(const String &param_name) const
      {
        String result(param_name);
        result.substitute('-', '_');
        return sectionName() + ":" + result;
      }

    protected:
      ParameterModifier parameter(
              const String &parameter_name,
              const ParamValue &default_value,
              const String &parameter_description);
      void flag(
              const String &parameter_name,
              const String &parameter_description);

      explicit ParameterSection(SiriusExportAlgorithm* enclose): enclose(enclose) {}

      virtual void parameters() = 0;
      virtual String sectionName() const = 0;

      SiriusExportAlgorithm *enclose;

    public:
      virtual ~ParameterSection() = default;

      DataValue getValue(const String &param_name) const
      {
        return enclose->param_.getValue(toFullParameter(param_name));
      }

      QStringList getCommandLine() const
      {
        QStringList result;
        for (const auto &pair : openms_to_sirius)
        {
          DataValue value = enclose->param_.getValue(pair.first);
          DataValue default_value = enclose->defaults_.getValue(pair.first);
          if (!value.isEmpty() && value != default_value)
          {
           String string_value = value.toString(true);
           if (string_value == "true")
           {
             result.push_back(String("--" + pair.second).toQString());
           }
           else if (string_value != "false")
           {
             result.push_back(String("--" + pair.second + "=" + string_value).toQString());
           }
          };
        }
        return result;
      }
    };

    using SiriusSubtool = ParameterSection;

    class Preprocessing : public ParameterSection
    {
      String sectionName() const override { return "preprocessing"; }
    public:
      explicit Preprocessing(SiriusExportAlgorithm *enclose) : ParameterSection(enclose) {}
      void parameters() override;
    };

    Preprocessing preprocessing;

    };
} // namespace OpenMS
