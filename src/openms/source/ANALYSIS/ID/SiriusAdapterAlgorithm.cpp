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

#include <OpenMS/ANALYSIS/ID/SiriusAdapterAlgorithm.h>

#include <OpenMS/CONCEPT/Exception.h>

#include <QtCore/QProcess>
#include <QDir>
#include <QDirIterator>

#include <fstream>
#include <include/OpenMS/DATASTRUCTURES/StringUtils.h>

namespace OpenMS
{
    using SiriusName = String;
    using OpenMSName = String;
    using DefaultValue = DataValue;
    using Description = String;

    SiriusAdapterAlgorithm::SiriusAdapterAlgorithm() :
      DefaultParamHandler("SiriusAdapterAlgorithm"),
      preprocessing(Preprocessing(this)),
      nightsky_config(Config(this)),
      nightsky_sirius(Sirius(this))
    {
      // Defines the Parameters for Preprocessing and NightSky subtools
      preprocessing.parameters();
      nightsky_config.parameters();
      nightsky_sirius.parameters();

      defaultsToParam_();
    }

    void SiriusAdapterAlgorithm::Preprocessing::parameters()
    {
      parameter(
              OpenMSName("filter_by_num_masstraces"),
              DefaultValue(1),
              Description("Number of mass traces each feature has to have to be included. "
                          "To use this parameter, setting the feature_only flag is necessary"))
          .withMinInt(1);

      parameter(
              OpenMSName("precursor_mz_tolerance"),
              DefaultValue(0.005),
              Description("Tolerance window for precursor selection (Feature selection in regard to the precursor)"));

      parameter(
              OpenMSName("precursor_mz_tolerance_unit"),
              DefaultValue("Da"),
              Description("Unit of the precursor_mz_tolerance"))
          .withValidStrings({"Da", "ppm"});

      parameter(
              OpenMSName("precursor_rt_tolerance"),
              DefaultValue(5),
              Description("Tolerance window (left and right) for precursor selection [seconds]"));

      //  defaults_.setValue".", ListUtils::create<String>("advanced"));
      parameter(
              OpenMSName("isotope_pattern_iterations"),
              DefaultValue(3),
              Description("Number of iterations that should be performed to extract the C13 isotope pattern. "
                          "If no peak is found (C13 distance) the function will abort. "
                          "Be careful with noisy data - since this can lead to wrong isotope patterns"));

      flag(OpenMSName("feature_only"),
           Description("Uses the feature information from in_featureinfo to reduce the search space to MS2 "
                          "associated with a feature"));

      flag(OpenMSName("no_masstrace_info_isotope_pattern"),
           Description("Use this flag if the masstrace information from a feature should be discarded "
                       "and the isotope_pattern_iterations should be used instead"));
    }


    void SiriusAdapterAlgorithm::Config::parameters()
    {
      parameter(
              SiriusName("NoiseIntensity"),
              DefaultValue(0),
              Description("median intensity of noise peaks"))
         .withMinInt(0);

      // TODO These parameters are currently missing
      // defaults_.setValue("sirius:isotope", "both", "how to handle isotope pattern data. Use 'score' to use them for ranking or 'filter' if you just want to remove candidates with bad isotope pattern. With 'both' you can use isotopes for filtering and scoring. Use 'omit' to ignore isotope pattern.");
      // defaults_.setValidStrings("sirius:isotope", ListUtils::create<String>("score,filter,both,omit"));
      // defaults_.setValue("sirius:top_n_hits", 10, "The number of top hits for each compound written to the CSI:FingerID output");
      // defaults_.setMinInt("sirius:top_n_hits", 1);
      // defaults_.setValue("sirius:cores", 1, "The number of cores SIRIUS is allowed to use on the system");
      // defaults_.setMinInt("sirius:cores", 1);
      // defaults_.setValue("sirius:auto_charge", "false", "Use this option if the charge of your compounds is unknown and you do not want to assume [M+H]+ as default. With the auto charge option SIRIUS will not care about charges and allow arbitrary adducts for the precursor peak.");
      // defaults_.setValidStrings("sirius:auto_charge", ListUtils::create<String>("true,false"));
      // defaults_.setValue("sirius:ion_tree", "false", "Print molecular formulas and node labels with the ion formula instead of the neutral formula", ListUtils::create<String>("advanced"));
      // defaults_.setValidStrings("sirius:ion_tree", ListUtils::create<String>("true,false"));
      // defaults_.setValue("sirius:most_intense_ms2", "false", "SIRIUS uses the fragmentation spectrum with the most intense precursor peak (for each spectrum)", ListUtils::create<String>("advanced"));
      // defaults_.setValidStrings("sirius:most_intense_ms2", ListUtils::create<String>("true,false"));
    }

    void SiriusAdapterAlgorithm::Sirius::parameters()
    {
      // --ppm-max=<ppmMax>
      // Maximum allowed mass deviation in ppm for decomposing masses.
      parameter(
              SiriusName("ppm-max"),
              DefaultValue(10),
              Description("Maximum allowed mass deviation in ppm for decomposing masses"));

      // --ppm-max-ms2=<ppmMaxMs2>
      // Maximum allowed mass deviation in ppm for decomposing masses in MS2.
      // If not specified, the same value as for the MS1 is used.
      parameter(
              SiriusName("ppm-max-ms2"),
              DefaultValue(10),
              Description("Maximum allowed mass deviation in ppm for decomposing masses in MS2."
                          "If not specified, the same value as for the MS1 is used."));

      // --tree-timeout=<treeTimeout>
      // Time out in seconds per fragmentation tree computations. 0 for an infinite amount of time.
      // Default: 0
      parameter(
              SiriusName("tree-timeout"),
              DefaultValue(0),
              Description("Time out in seconds per fragmentation tree computations. 0 for an infinite amount of time"))
         .withMinInt(0);

      //--compound-timeout=<instanceTimeout>
      // Maximal computation time in seconds for a single compound. 0 for an infinite amount of time.
      // Default: 0
      parameter(
              SiriusName("compound-timeout"),
              DefaultValue(0),
              Description(
                      "Maximal computation time in seconds for a single compound. 0 for an infinite amount of time."))
         .withMinInt(0);

      // --no-recalibration
      // Disable Recalibration of input Spectra
      flag(
              SiriusName("no-recalibration"),
              Description("Disable Recalibration of input Spectra"));

      // -p, --profile=<profile>
      // Name of the configuration profile. Some of the default profiles are: 'qtof', 'orbitrap', 'fticr'.
      parameter(
              SiriusName("profile"),
              DefaultValue("qtof"),
              Description("Name of the configuration profile"))
         .withValidStrings({"qtof", "orbitrap", "fticr"});

      // -f, --formula, --formulas=<formula>
      // Specify the neutral molecular formula of the measured compound to compute its tree or a list of candidate
      // formulas the method should discriminate. Omit this option if you want to consider all possible
      // molecular formulas
      parameter(
              SiriusName("formula"),
              DefaultValue(""),
              Description("Specify the neutral molecular formula of the measured "
                          "compound to compute its tree or a list of candidate "
                          "formulas the method should discriminate. Omit this "
                          "option if you want to consider all possible molecular formulas"));

      // -I, --ions-enforced=<ionsEnforced>
      // the iontype/adduct of the MS/MS data.
      // Example: [M+H]+, [M-H]-, [M+Cl]-, [M+Na]+, [M]+.
      // You can also provide a comma separated list of adducts.
      parameter(
              SiriusName("ions-enforced"),
              DefaultValue(-1),
              Description("the iontype/adduct of the MS/MS data. Example: [M+H]+, "
                          "[M-H]-, [M+Cl]-, [M+Na]+, [M]+. You can also provide a comma separated list of adducts"));

      //-c, --candidates=<numberOfCandidates>
      //         Number of formula candidates in the output.
      parameter(
              SiriusName("candidates"),
              DefaultValue(5),
              Description("The number of formula candidates in the SIRIUS output"))
         .withMinInt(1);

      // --candidates-per-ion=<numberOfCandidatesPerIon>
      //        Minimum number of candidates in the output for each
      // ionization. Set to force output of results for each
      // possible ionization, even if not part of highest
      // ranked results.
      parameter(
              SiriusName("candidates-per-ion"),
              DefaultValue(-1),
              Description("Minimum number of candidates in the output for each "
                          "ionization. Set to force output of results for each "
                          "possible ionization, even if not part of highest "
                          "ranked results. -1 omits parameter in Sirius."));

      // -e, --elements-considered=<detectableElements>
      //        Set the allowed elements for rare element detection.
      //        Write SBrClBSe to allow the elements S,Br,Cl,B and Se.
      parameter(
              SiriusName("elements-considered"),
              DefaultValue(""),
              Description("Set the allowed elements for rare element detection. "
                          "Write SBrClBSe to allow the elements S,Br,Cl,B and Se."));

      // -E, --elements-enforced=<enforcedElements>
      // Enforce elements for molecular formula determination.
      // Write CHNOPSCl to allow the elements C, H, N, O, P, S and Cl. Add numbers in brackets to restrict the
      // minimal and maximal allowed occurrence of these elements: CHNOP[5]S[8]Cl[1-2]. When one number is
      // given then it is interpreted as upper bound. Default is CHNOP
      parameter(
              SiriusName("elements-enforced"),
              DefaultValue(""),
              Description("Enforce elements for molecular formula determination. "
                          "Write CHNOPSCl to allow the elements C, H, N, O, P, S "
                          "and Cl. Add numbers in brackets to restrict the "
                          "minimal and maximal allowed occurrence of these "
                          "elements: CHNOP[5]S[8]Cl[1-2]. When one number is "
                          "given then it is interpreted as upper bound. Default is CHNOP"));

      // --no-isotope-score=<isotopeHandling>
      //         Disable isotope pattern score.
      flag(SiriusName("no-isotope-score"),
           Description("Disable isotope pattern score."));

      // --no-isotope-filter
      // Disable molecular formula filter. When filtering is enabled, molecular formulas are excluded if their
      // theoretical isotope pattern does not match the theoretical one, even if their MS/MS pattern has high score.
      flag(SiriusName("no-isotope-filter"),
           Description("Disable molecular formula filter. When filtering is enabled, molecular formulas are "
                      "excluded if their theoretical isotope pattern does not match the theoretical one, even if "
                      "their MS/MS pattern has high score."));

      // -i, --ions-considered=<ionsConsidered>
      //  the iontype/adduct of the MS/MS data. Example: [M+H]+,
      // [M-H]-, [M+Cl]-, [M+Na]+, [M]+. You can also provide a
      // comma separated list of adducts.
      parameter(
              SiriusName("ions-considered"),
              DefaultValue(""),
              Description("the iontype/adduct of the MS/MS data. "
                          "Example: [M+H]+, [M-H]-, [M+Cl]-, [M+Na]+, [M]+. "
                          "You can also provide a comma separated list of adducts."));
    }


    void SiriusAdapterAlgorithm::updateExistingParameter(const OpenMS::Param &param)
    {
      for (auto it = param.begin(); it != param.end(); ++it)
      {
        const String name = it.getName();
        if (hasFullNameParameter(name))
        {
          vector<String> tags(it->tags.begin(), it->tags.end());
          param_.setValue(name, it->value, it->description, tags);
        }
      }
    }

    bool SiriusAdapterAlgorithm::hasFullNameParameter(const OpenMS::String &name) const
    {
      for (auto it = param_.begin(); it != param_.end(); ++it)
      {
        if (it.getName() == name)
        {
          return true;
        }
      }
      return false;
    }

    String SiriusAdapterAlgorithm::determineSiriusExecutable(String &executable)
    { 
      // if executable was not provided
      if (executable.empty())
      {
        const String& qsiriuspathenv = QProcessEnvironment::systemEnvironment().value("SIRIUS_PATH");
        if (qsiriuspathenv.empty())
        {
          throw Exception::InvalidValue(__FILE__,
                                        __LINE__, 
                                        OPENMS_PRETTY_FUNCTION, 
                                        "FATAL: Executable of Sirius could not be found. Please either use SIRIUS_PATH env variable or provide with -executable",
                                         "");
        }
        executable = qsiriuspathenv;
      }
      const String exe = QFileInfo(executable.toQString()).canonicalFilePath().toStdString();
      LOG_WARN << "Executable is: " + exe << std::endl;
      return exe;
    }

    SiriusAdapterAlgorithm::SiriusTmpStruct SiriusAdapterAlgorithm::constructSiriusTmpStruct()
    {
      SiriusTmpStruct tmp_struct;
      QString base_dir = File::getTempDirectory().toQString();
      tmp_struct.tmp_dir = String(QDir(base_dir).filePath(File::getUniqueName().toQString()));
      tmp_struct.tmp_ms_file = QDir(base_dir).filePath((File::getUniqueName() + ".ms").toQString());
      tmp_struct.tmp_out_dir = QDir(tmp_struct.tmp_dir.toQString()).filePath("sirius_out");

      return tmp_struct;
    } 
    
    void SiriusAdapterAlgorithm::preprocessingSirius(const String& featureinfo,
                                                     const MSExperiment& spectra,
                                                     std::vector<FeatureMap>& v_fp,
                                                     KDTreeFeatureMaps& fp_map_kd,
                                                     FeatureMapping::FeatureToMs2Indices& feature_mapping)
    {
      // if fileparameter is given and should be not empty
      if (!featureinfo.empty())
      {
        if (File::exists(featureinfo) && !File::empty(featureinfo))
        {
          // read featureXML          
          FeatureXMLFile fxml;
          FeatureMap feature_map;
          fxml.load(featureinfo, feature_map);
          
          UInt num_masstrace_filter = getFilterByNumMassTraces();
          double precursor_mz_tol = getPrecursorMzTolerance();
          double precursor_rt_tol = getPrecursorRtTolerance();

          if (num_masstrace_filter != 1 && !isFeatureOnly())
          {
            num_masstrace_filter = 1;
            LOG_WARN << "Parameter: filter_by_num_masstraces, was set to 1 to retain the adduct information for all MS2 spectra, if available. Please use the masstrace filter in combination with feature_only." << std::endl;
          }

          // filter feature by number of masstraces
          auto map_it = remove_if(feature_map.begin(), feature_map.end(),
                                  [&num_masstrace_filter](const Feature &feat) -> bool
                                  {
                                    unsigned int n_masstraces = feat.getMetaValue("num_of_masstraces");
                                    return n_masstraces < num_masstrace_filter;
                                  });
          feature_map.erase(map_it, feature_map.end());
  
          v_fp.push_back(feature_map);
          fp_map_kd.addMaps(v_fp);
  
          // mapping of MS2 spectra to features
          feature_mapping = FeatureMapping::assignMS2IndexToFeature(spectra,
                                                                    fp_map_kd,
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

    void SiriusAdapterAlgorithm::logFeatureSpectraNumber(const String &featureinfo,
                                                         const FeatureMapping::FeatureToMs2Indices &feature_mapping,
                                                         const MSExperiment &spectra)
    {
      // number of features to be processed
      if (isFeatureOnly() && !featureinfo.empty())
      {
        LOG_INFO << "Number of features to be processed: " << feature_mapping.assignedMS2.size() << std::endl;
      }
      else if (!featureinfo.empty())
      {
        LOG_INFO << "Number of features to be processed: " << feature_mapping.assignedMS2.size() << std::endl;
        LOG_INFO << "Number of additional MS2 spectra to be processed: " << feature_mapping.unassignedMS2.size() << std::endl;
      } 
      else
      {
        long count_ms2 = count_if(spectra.begin(), spectra.end(),
                [](const MSSpectrum &spectrum) { return spectrum.getMSLevel() == 2; });

        LOG_INFO << "Number of MS2 spectra to be processed: " << count_ms2 << std::endl;
      }
    }

    // tmp_msfile (store), all parameters, out_dir (tmpstructure)
    const std::vector<String> SiriusAdapterAlgorithm::callSiriusQProcess(const String& tmp_ms_file,
                                                                         const String& tmp_out_dir,
                                                                               String& executable,
                                                                         const String& out_csifingerid)
    {
      // Get the command line parameters from all the subtools
      QStringList config_params = nightsky_config.getCommandLine();
      QStringList sirius_params = nightsky_sirius.getCommandLine();

      const bool run_csifingerid = ! out_csifingerid.empty();

      // Structure of the command line passed to NightSky
      QStringList command_line = QStringList({tmp_ms_file.toQString(), "config"}) + config_params + sirius_params + QStringList("sirius");

      if (run_csifingerid)
      {
        command_line << "fingerid";
      }

      LOG_INFO << "Running NightSky with the following command line parameters: " << endl;
      for (const auto &param: command_line)
      {
        LOG_INFO << param.toStdString() << " ";
      }
      LOG_INFO << endl;


        // the actual process
        QProcess qp;
        QString executable_qstring = SiriusAdapterAlgorithm::determineSiriusExecutable(executable).toQString();

        //since library paths are relative to sirius executable path
        qp.setWorkingDirectory(File::path(executable).toQString());
        qp.start(executable_qstring, command_line); // does automatic escaping etc... start
        const bool success = qp.waitForFinished(-1);

        if (!success || qp.exitStatus() != 0 || qp.exitCode() != 0)
        {
          LOG_WARN << "FATAL: External invocation of Sirius failed. Standard output and error were:" << std::endl;
          const QString sirius_stdout(qp.readAllStandardOutput());
          const QString sirius_stderr(qp.readAllStandardError());
          LOG_WARN << String(sirius_stdout) << std::endl;
          LOG_WARN << String(sirius_stderr) << std::endl;
          LOG_WARN << String(qp.exitCode()) << std::endl;
          qp.close();

          throw Exception::InvalidValue(__FILE__,
                                        __LINE__,
                                        OPENMS_PRETTY_FUNCTION,
                                        "FATAL: External invocation of Sirius failed!",
                                         "");
        }
        qp.close();

       //extract path to subfolders (sirius internal folder structure)
       std::vector<String> subdirs;
        QDirIterator it(tmp_out_dir.toQString(), QDir::Dirs | QDir::NoDotAndDotDot, QDirIterator::NoIteratorFlags);
        while (it.hasNext())
        {
          subdirs.push_back(it.next());
        }
      return subdirs;
    }


    SiriusAdapterAlgorithm::ParameterModifier SiriusAdapterAlgorithm::ParameterSection::parameter(
            const String &parameter_name,
            const DataValue &default_value,
            const String &parameter_description)
    {
      const String full_parameter = toFullParameter(parameter_name);
      openms_to_sirius[full_parameter] = parameter_name;
      enclose->defaults_.setValue(full_parameter, default_value, parameter_description);
      return ParameterModifier(full_parameter, enclose);
    }

  void SiriusAdapterAlgorithm::ParameterSection::flag(
            const OpenMS::String &parameter_name,
            const OpenMS::String &parameter_description)
    {
      parameter(parameter_name, DefaultValue("false"), parameter_description)
        .withValidStrings({"true", "false"});
    }
} // namespace OpenMS

/// @endcond
