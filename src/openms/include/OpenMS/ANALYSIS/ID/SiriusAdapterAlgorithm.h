// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2019.
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
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#pragma once 

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/ANALYSIS/QUANTITATION/KDTreeFeatureMaps.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureMapping.h>

#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>

#include <unordered_map>
#include <QString>

using namespace std;

namespace OpenMS
{
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
      bool isNoMasstraceInfoIsotopePattern() const { return preprocessing.getValue("no_masstrace_info_isotope_pattern").toBool(); };
      int getIsotopePatternIterations() const { return  preprocessing.getValue("isotope_pattern_iterations"); }

      /*
       * Accessors for Sirius Parameters
       */
      int getNumberOfCandidates() const { return nightsky_sirius.getValue("candidates");  }

      /**
       * Updates all parameters that already exist in this DefaultParamHandler
       * with the values provided by the input param object.
       *
       * @param param The Param object supplying updated parameter values. Keys that exist in the param
       *        parameter, but not in this DefaultParamHander, are ignored.
       */
      void updateExistingParameter(const Param &param);

      /**
       * Checks whether this DefaultParamHandler has a ParamEntry with the provided name.
       * @param name The name of the ParamEntry that should be checked for its existence in this DefaultParamHandler
       * @return Whether this DefaultParamHandler has an ParamEntry for the provided name.
       */
      bool hasFullNameParameter(const String &name) const;

      /// Struct for temporary folder structure
      struct SiriusTmpStruct
      {
        String tmp_dir;
        String tmp_ms_file;
        String tmp_out_dir;
      };

      /**
        @brief Construct temporary folder structure for SIRIUS (SiriusTmpStruct)

        @return SiriusTmpStruct
      */
      static SiriusAdapterAlgorithm::SiriusTmpStruct constructSiriusTmpStruct();

      /**
      @brief Checks if the provided String points to a valid SIRIUS executable, otherwise tries
       to select the executable from the environment

      @return Path to SIRIUS executable

      @param executable Path to the potential executable
      */
      static String determineSiriusExecutable(String &executable);

      /**
      @brief Preprocessing needed for SIRIUS

      @return FeatureToMS2Indices
            
      Filter number of masstraces and perform feature mapping.

      @param featureinfo: Path to featureXML
      @param spectra: Input of MSExperiment with spectra information
      @param fp_map_kd: KDTree used for query and match spectra with features
      @param sirius_algo: Parameters needed for preprocessing
      @param feature_mapping: Empty FeatureToMs2Indices
      */
      void preprocessingSirius(const String& featureinfo,
                               const MSExperiment& spectra,
                               vector<FeatureMap>& v_fp,
                               KDTreeFeatureMaps& fp_map_kd,
                               FeatureMapping::FeatureToMs2Indices& feature_mapping);

      /**
      @brief logs number of features and spectra used

      Prints the number of features and spectra used (LOG_INFO)

      @param featureinfo: Path to featureXML
      @param feature_mapping: FeatureToMs2Indices with feature mapping
      @param spectra: Input of MSExperiment with spectra information
      */
      void logFeatureSpectraNumber(const String &featureinfo,
                                   const FeatureMapping::FeatureToMs2Indices &feature_mapping,
                                   const MSExperiment &spectra);

      /**
      @brief Call SIRIUS with QProcess

      @param tmp_ms_file: path to temporary .ms file
      @param tmp_out_dir: path to temporary output folder
      @param executable: path to executable
      @param out_csifingerid: path to CSI:FingerID output (can be empty).

      @return Vector with paths to a compound
      */
      const vector<String> callSiriusQProcess(const String& tmp_ms_file,
                                              const String& tmp_out_dir,
                                              String& executable,
                                              const String& out_csifingerid);

    private:

    class ParameterModifier
    {
      const String openms_param_name;
      SiriusAdapterAlgorithm *enclose;

    public:
      explicit ParameterModifier(const String &param_name, SiriusAdapterAlgorithm *enclose) :
              openms_param_name(param_name), enclose(enclose) {}

      void withValidStrings(initializer_list<String> choices)
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
              const DataValue &default_value,
              const String &parameter_description);
      void flag(
              const String &parameter_name,
              const String &parameter_description);

      explicit ParameterSection(SiriusAdapterAlgorithm *enclose): enclose(enclose) {};
      virtual void parameters() = 0;
      virtual String sectionName() const = 0;

      SiriusAdapterAlgorithm *enclose;

    public:
      DataValue getValue(const String &param_name) const
      {
        return enclose->param_.getValue(toFullParameter(param_name));
      }

      QStringList getCommandLine()
      {
        const DataValue omit_integer(-1);
        const DataValue omit_string("");

        QStringList result;
        for (const auto &pair : openms_to_sirius)
        {
          DataValue value = enclose->param_.getValue(pair.first);

          if (!value.isEmpty() && value != omit_integer && value != omit_string)
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

    using NightSkySubtool = ParameterSection;

    class Preprocessing : public ParameterSection
    {
      String sectionName() const override { return "preprocessing"; }
    public:
      explicit Preprocessing(SiriusAdapterAlgorithm *enclose) : ParameterSection(enclose) {}
      void parameters() override;
    };

    class Config final : public NightSkySubtool
    {
      String sectionName() const override { return "config"; }
    public:
      explicit Config(SiriusAdapterAlgorithm *enclose) : ParameterSection(enclose) {}
      void parameters() override;
    };

    class Sirius final : public NightSkySubtool
    {
      String sectionName() const override { return "sirius"; }
    public:
      explicit Sirius(SiriusAdapterAlgorithm *enclose) : ParameterSection(enclose) {}
      void parameters() override;
    };

    Preprocessing preprocessing;
    Config nightsky_config;
    Sirius nightsky_sirius;
    };
} // namespace OpenMS
