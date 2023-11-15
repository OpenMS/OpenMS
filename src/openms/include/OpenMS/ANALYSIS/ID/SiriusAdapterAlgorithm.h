// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka $
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

  class OPENMS_DLLAPI SiriusAdapterAlgorithm : public DefaultParamHandler
    {
    public:
      /// default constructor
      SiriusAdapterAlgorithm();

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
       * @brief Accessors for Sirius Parameters
       */

      int getNumberOfSiriusCandidates() const 
      {
        int number_of_candidates = sirius.getValue("candidates");
        // default for SiriusAdapter is -1 to not pass a value to command and use SIRIUS 5 default (10)
        // therefore 10 needs to be returned in this case
        if (number_of_candidates == -1)
        {
          return 10;
        }
        else
        {
          return number_of_candidates;
        }
      }

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

      /// Struct for temporary folder structure
      class OPENMS_DLLAPI SiriusTemporaryFileSystemObjects
      {
      public:

        /// Construct temporary folder structure for SIRIUS (SiriusTemporaryFileSystemObjects)
        SiriusTemporaryFileSystemObjects(int debug_level);

        /// Destructor of SiriusTemporaryFileSystemObjects based on debug level
        ~SiriusTemporaryFileSystemObjects();

        const String& getTmpDir() const;
        const String& getTmpOutDir() const;
        const String& getTmpMsFile() const;
      
      private:
        int debug_level_;

        String tmp_dir_;
        String tmp_ms_file_;
        String tmp_out_dir_;
      };

      /**
      @brief Checks if the provided String points to a valid SIRIUS executable, otherwise tries
       to select the executable from the environment

      @return Path to SIRIUS executable

      @param executable Path to the potential executable
      */
      static String determineSiriusExecutable(String& executable);

      /**

      @brief Sort function using the extracted scan_index from the sirius workspace file path

      @return Vector of sorted sirius workspace paths based on the scan_index

      */
      static void sortSiriusWorkspacePathsByScanIndex(std::vector<String>& subdirs);


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

      /**
      @brief Log in to Sirius with personal user account (required in Sirius >= 5).

      @param email User account E-Mail.
      @param password User account password.
      */
      void logInSiriusAccount(String& executable, const String& email, const String& password) const;

      /**
      @brief Call SIRIUS with QProcess

      @param tmp_ms_file path to temporary .ms file
      @param tmp_out_dir path to temporary output folder
      @param executable path to executable
      @param out_csifingerid path to CSI:FingerID output (can be empty).

      @return Vector with paths to a compound
      */
      const vector<String> callSiriusQProcess(const String& tmp_ms_file,
                                              const String& tmp_out_dir,
                                              String& executable,
                                              const String& out_csifingerid,
                                              const bool decoy_generation) const;

    private:
    class ParameterModifier
    {
      const String openms_param_name;
      SiriusAdapterAlgorithm *enclose;

    public:
      explicit ParameterModifier(const String &param_name, SiriusAdapterAlgorithm *enclose) :
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

      explicit ParameterSection(SiriusAdapterAlgorithm* enclose): enclose(enclose) {}

      virtual void parameters() = 0;
      virtual String sectionName() const = 0;

      SiriusAdapterAlgorithm *enclose;

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
      explicit Preprocessing(SiriusAdapterAlgorithm *enclose) : ParameterSection(enclose) {}
      void parameters() override;
    };

    class Project final : public SiriusSubtool
    {
      String sectionName() const override { return "project"; }
    public:
      explicit Project(SiriusAdapterAlgorithm *enclose) : ParameterSection(enclose) {}
      void parameters() override;
    };

    class Sirius final : public SiriusSubtool
    {
      String sectionName() const override { return "sirius"; }
    public:
      explicit Sirius(SiriusAdapterAlgorithm *enclose) : ParameterSection(enclose) {}
      void parameters() override;
    };

    class Fingerid final : public SiriusSubtool
    {
      String sectionName() const override { return "fingerid"; }
    public:
      explicit Fingerid(SiriusAdapterAlgorithm *enclose) : ParameterSection(enclose) {}
      void parameters() override;
    };

    class Passatutto final : public SiriusSubtool
    {
      String sectionName() const override { return "passatutto"; }
    public:
      explicit Passatutto(SiriusAdapterAlgorithm *enclose) : ParameterSection(enclose) {}
      void parameters() override;
     };

    Preprocessing preprocessing;
    Project project;
    Sirius sirius;
    Fingerid fingerid;
    Passatutto passatutto;

    };
} // namespace OpenMS
