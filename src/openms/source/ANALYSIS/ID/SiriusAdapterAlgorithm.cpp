// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka, Lukas Zimmermann $
// --------------------------------------------------------------------------

#include <boost/foreach.hpp> // must be first, otherwise Q_FOREACH macro will wreak havoc

#include <OpenMS/ANALYSIS/ID/SiriusAdapterAlgorithm.h>

#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/DATAACCESS/SiriusMzTabWriter.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/SYSTEM/File.h>

#include <QDir>
#include <QDirIterator>
#include <QString>
#include <QtCore/QProcess>
#include <sstream>

namespace OpenMS
{
  // ###################
  // Set subtool parameters
  // ###################

    using ProjectName = String;
    using SiriusName = String;
    using FingeridName = String;
    using PassatuttoName = String;
    using OpenMSName = String;
    using DefaultValue = ParamValue;
    using Description = String;

    SiriusAdapterAlgorithm::SiriusAdapterAlgorithm() :
      DefaultParamHandler("SiriusAdapterAlgorithm"),
      preprocessing(Preprocessing(this)),
      project(Project(this)),
      sirius(Sirius(this)),
      fingerid(Fingerid(this)),
      passatutto(Passatutto(this))
    {
      // Defines the Parameters for preprocessing and SIRIUS subtools
      preprocessing.parameters();
      project.parameters();
      sirius.parameters();
      fingerid.parameters();
      passatutto.parameters();

      defaults_.setValue("read_sirius_stdout", "false", "Read and print the standard output and error of the Sirius executable, even if it succeeds.", {"advanced"});
      defaults_.setValidStrings("read_sirius_stdout",{"true","false"});

      defaultsToParam_();
    }

    void SiriusAdapterAlgorithm::Preprocessing::parameters()
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

    void SiriusAdapterAlgorithm::Project::parameters()
    {
      parameter(
          ProjectName("maxmz"),
          DefaultValue(-1),
          Description("Just consider compounds with a precursor mz lower or equal\n"
                      "this maximum mz. All other compounds in the input file\n"
                      "are ignored.")
      );

      parameter(
          ProjectName("processors"),
          DefaultValue(1),
          Description("Number of cpu cores to use. If not specified SIRIUS uses all available cores.")
      );

      parameter(
          ProjectName("loglevel"),
          DefaultValue(""),
          Description("Set logging level of the Jobs SIRIUS will execute.\n"
                      "Valid values: SEVERE, WARNING, INFO, FINER, ALL\n"
                      "Default: WARNING")
      );

      flag(
          ProjectName("ignore-formula"),
          Description("Ignore given molecular formula in internal .ms format, while processing.")
      );

      flag(
          ProjectName("q"),
          Description("Suppress shell output")
      );
    }

    void SiriusAdapterAlgorithm::Sirius::parameters() // original sirius 5.6.1 parameters as defaults, except for timeouts
    {
      parameter(
                 SiriusName("ppm-max"),
                 DefaultValue(10.0),
                 Description("Maximum allowed mass deviation in ppm for decomposing masses [ppm].")
               );

      parameter(
                 SiriusName("ppm-max-ms2"),
                 DefaultValue(10.0),
                 Description("Maximum allowed mass deviation in ppm for decomposing masses in MS2 [ppm]."
                             "If not specified, the same value as for the MS1 is used. ")
                );

      parameter(
                 SiriusName("tree-timeout"),
                 DefaultValue(100), // original default = 0
                 Description("Time out in seconds per fragmentation tree computations. 0 for an infinite amount of time")
               ).withMinInt(0);

      parameter(
                 SiriusName("compound-timeout"),
                 DefaultValue(100), // original default = 0
                 Description("Maximal computation time in seconds for a single compound. 0 for an infinite amount of time.")
               ).withMinInt(0);

      flag(
            SiriusName("no-recalibration"),
            Description("Disable recalibration of input spectra")
          );

      parameter(
                 SiriusName("profile"),
                 DefaultValue("default"),
                 Description("Name of the configuration profile")
               ).withValidStrings({"default", "qtof", "orbitrap", "fticr"});

      parameter(
                 SiriusName("formulas"),
                 DefaultValue(""),
                 Description("Specify the neutral molecular formula of the measured "
                             "compound to compute its tree or a list of candidate "
                             "formulas the method should discriminate. Omit this "
                             "option if you want to consider all possible molecular formulas")
               );

      parameter(
                 SiriusName("ions-enforced"),
                 DefaultValue(""),
                 Description("the iontype/adduct of the MS/MS data. Example: [M+H]+, "
                             "[M-H]-, [M+Cl]-, [M+Na]+, [M]+. You can also provide a comma separated list of adducts")
               );

      parameter(
                 SiriusName("candidates"),
                 DefaultValue(10),
                 Description("The number of formula candidates in the SIRIUS output")
               ).withMinInt(0);

      parameter(
                 SiriusName("candidates-per-ion"),
                 DefaultValue(1),
                 Description("Minimum number of candidates in the output for each "
                             "ionization. Set to force output of results for each "
                             "possible ionization, even if not part of highest "
                             "ranked results.")
                );

      parameter(
                 SiriusName("elements-considered"),
                 DefaultValue("SBrClBSe"),
                 Description("Set the allowed elements for rare element detection. "
                             "Write SBrClBSe to allow the elements S,Br,Cl,B and Se."));

      parameter(
                 SiriusName("elements-enforced"),
                 DefaultValue("CHNOP"),
                 Description("Enforce elements for molecular formula determination. "
                             "Write CHNOPSCl to allow the elements C, H, N, O, P, S "
                             "and Cl. Add numbers in brackets to restrict the "
                             "minimal and maximal allowed occurrence of these "
                             "elements: CHNOP[5]S[8]Cl[1-2]. When one number is "
                             "given then it is interpreted as upper bound.")
                );

      flag(
            SiriusName("no-isotope-score"),
            Description("Disable isotope pattern score.")
          );

      flag(
            SiriusName("no-isotope-filter"),
            Description("Disable molecular formula filter. When filtering is enabled, molecular formulas are "
                        "excluded if their theoretical isotope pattern does not match the theoretical one, even if "
                        "their MS/MS pattern has high score.")
          );

      parameter(
                 SiriusName("ions-considered"),
                 DefaultValue("[M+H]+,[M+K]+,[M+Na]+,[M+H-H2O]+,[M+H-H4O2]+,[M+NH4]+,[M-H]-,[M+Cl]-,[M-H2O-H]-,[M+Br]-"),
                 Description("the iontype/adduct of the MS/MS data. "
                             "Example: [M+H]+, [M-H]-, [M+Cl]-, [M+Na]+, [M]+. "
                             "You can also provide a comma separated list of adducts.")
                );

      parameter(
                 SiriusName("db"),
                 DefaultValue("none"),
                 Description("Search formulas in the Union of the given "
                              "databases db-name1,db-name2,db-name3. If no database is given all possible "
                              "molecular formulas will be respected (no database "
                              "is used). "
                              "Example: possible DBs: ALL,BIO,PUBCHEM,MESH,HMDB,"
                              "KNAPSACK,CHEBI,PUBMED,KEGG,HSDB,MACONDA,METACYC,"
                              "GNPS,ZINCBIO,UNDP,YMDB,PLANTCYC,NORMAN,ADDITIONAL,"
                              "PUBCHEMANNOTATIONBIO,PUBCHEMANNOTATIONDRUG,"
                              "PUBCHEMANNOTATIONSAFETYANDTOXIC,"
                              "PUBCHEMANNOTATIONFOOD,KEGGMINE,ECOCYCMINE,"
                              "YMDBMINE")
                );

      parameter(
                 SiriusName("ions-enforced"),
                 DefaultValue(""),
                 Description("The iontype/adduct of the MS/MS data. Example: [M+H]+, \n"
                             "[M-H]-, [M+Cl]-, [M+Na]+, [M]+. You can also provide a \n"
                             "comma separated list of adducts.")
               );

      parameter(
                 SiriusName("solver"),
                 DefaultValue("CLP"),
                 Description("For GUROBI and CPLEX environment variables need to be configured. \n"
                             "(see SIRIUS manual: https://boecker-lab.github.io/docs.sirius.github.io/install/).")
               );
    }

    void SiriusAdapterAlgorithm::Fingerid::parameters()
    {
      parameter(
                 FingeridName("db"),
                 DefaultValue(""), // default = BIO
                 Description("Search structures in the Union of the given "
                              "databases db-name1,db-name2,db-name3. If no database is given all possible "
                              "molecular formulas will be respected (no database "
                              "is used). "
                              "Example: possible DBs: ALL,BIO,PUBCHEM,MESH,HMDB,"
                              "KNAPSACK,CHEBI,PUBMED,KEGG,HSDB,MACONDA,METACYC,"
                              "GNPS,ZINCBIO,UNDP,YMDB,PLANTCYC,NORMAN,ADDITIONAL,"
                              "PUBCHEMANNOTATIONBIO,PUBCHEMANNOTATIONDRUG,"
                              "PUBCHEMANNOTATIONSAFETYANDTOXIC,"
                              "PUBCHEMANNOTATIONFOOD,KEGGMINE,ECOCYCMINE,"
                              "YMDBMINE")
                );
    }

    void SiriusAdapterAlgorithm::Passatutto::parameters()
    {
    }

    void SiriusAdapterAlgorithm::updateExistingParameter(const OpenMS::Param &param)
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

    String SiriusAdapterAlgorithm::determineSiriusExecutable(String& executable)
    { 
      // if executable was not provided
      if (executable.empty())
      {
        const char* sirius_env_var = std::getenv("SIRIUS_PATH"); // returns nullptr if not found
        if (sirius_env_var == nullptr)
        {
            throw Exception::InvalidValue(__FILE__,
                                __LINE__,
                                OPENMS_PRETTY_FUNCTION,
                                "FATAL: Executable of SIRIUS could not be found. Please either use SIRIUS_PATH env variable, add the Sirius directory to our PATH or provide the executable with -sirius_executable",
                                "");
        }
        const std::string sirius_path(sirius_env_var);
        executable = sirius_path;
      }
      String exe = QFileInfo(executable.toQString()).canonicalFilePath().toStdString();
      return exe;
    }
    
    SiriusAdapterAlgorithm::SiriusTemporaryFileSystemObjects::SiriusTemporaryFileSystemObjects(int debug_level)
    {
      QString base_dir = File::getTempDirectory().toQString();
      tmp_dir_ = String(QDir(base_dir).filePath(File::getUniqueName().toQString()));
      tmp_ms_file_ = QDir(base_dir).filePath((File::getUniqueName() + ".ms").toQString());
      tmp_out_dir_ = QDir(tmp_dir_.toQString()).filePath("sirius_out");
      debug_level_ = debug_level;
    }

    SiriusAdapterAlgorithm::SiriusTemporaryFileSystemObjects::~SiriusTemporaryFileSystemObjects()
    {
      constexpr int debug_threshold = 9;

      // clean tmp directory if debug level < debug threshold
      if (debug_level_ >= debug_threshold)
      {
        OPENMS_LOG_DEBUG << "Keeping temporary files in directory " << tmp_dir_ << " and msfile at this location "<< tmp_ms_file_ << ". Set debug level lower than " << debug_threshold << " to remove them." << std::endl;
      }
      else
      {
        if (!tmp_dir_.empty())
        {
          OPENMS_LOG_DEBUG << "Deleting temporary directory " << tmp_dir_ << ". Set debug level to " << debug_threshold << " or higher to keep it." << std::endl;
          File::removeDir(tmp_dir_.toQString());
        }
        if (!tmp_ms_file_.empty())
        {
          OPENMS_LOG_DEBUG << "Deleting temporary msfile " << tmp_ms_file_ << ". Set debug level to " << debug_threshold << " or higher to keep it." << std::endl;
          File::remove(tmp_ms_file_);
        }
      }
    }

    // ################
    // Algorithm
    // ################

    const String& SiriusAdapterAlgorithm::SiriusTemporaryFileSystemObjects::getTmpDir() const
    {
      return tmp_dir_;
    }

    const String& SiriusAdapterAlgorithm::SiriusTemporaryFileSystemObjects::getTmpOutDir() const
    {
      return tmp_out_dir_;
    }

    const String& SiriusAdapterAlgorithm::SiriusTemporaryFileSystemObjects::getTmpMsFile() const
    {
      return tmp_ms_file_;
    }

    class OPENMS_DLLAPI SiriusWorkspaceIndex
    {
    public:
      int array_index, scan_index;

      SiriusWorkspaceIndex(int array_index, int scan_index) : array_index {array_index}, scan_index {scan_index} {}
    };

     void  SiriusAdapterAlgorithm::sortSiriusWorkspacePathsByScanIndex(std::vector<String>& subdirs)
    {
      std::vector<String> sorted_subdirs;
      std::vector<SiriusWorkspaceIndex> indices;

      for (size_t i = 0; i < subdirs.size(); i++)
      {
        indices.emplace_back(i, SiriusMzTabWriter::extractScanIndex(subdirs[i]));
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

    void SiriusAdapterAlgorithm::preprocessingSirius(const String& featureinfo,
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

    void SiriusAdapterAlgorithm::logFeatureSpectraNumber(const String& featureinfo,
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

    // ################
    // Algorithm
    // ################

    void SiriusAdapterAlgorithm::logInSiriusAccount(String& executable, const String& email, const String& password) const
    {
      // sirius login --email=email --password=password
      QString executable_qstring = SiriusAdapterAlgorithm::determineSiriusExecutable(executable).toQString();
      QStringList command_line{"login", String("--email="+email).toQString(), String("--password="+password).toQString()};
      
      // start QProcess with sirius login command
      QProcess qp;
      qp.start(executable_qstring, command_line);

      // print executed command as info
      std::stringstream ss;
      ss << "Executing command: " << executable_qstring.toStdString();
      for (const QString& it : command_line)
      {
        ss << " " << it.toStdString();
      }
      OPENMS_LOG_INFO << ss.str() << std::endl;

      // wait until process finished
      qp.waitForFinished(-1);

      // always print process stdout (info) and stderror (warning)
      const QString sirius_stdout(qp.readAllStandardOutput());
      const QString sirius_stderr(qp.readAllStandardError());
      OPENMS_LOG_INFO << String(sirius_stdout) << std::endl;
      OPENMS_LOG_WARN << String(sirius_stderr) << std::endl;
      
      qp.close();
    }

    // tmp_msfile (store), all parameters, out_dir (tmpstructure)
    const std::vector<String> SiriusAdapterAlgorithm::callSiriusQProcess(const String& tmp_ms_file,
                                                                         const String& tmp_out_dir,
                                                                         String& executable,
                                                                         const String& out_csifingerid,
                                                                         const bool decoy_generation) const

    {
      // get the command line parameters from all the subtools
      QStringList project_params = project.getCommandLine();
      QStringList sirius_params = sirius.getCommandLine();
      QStringList fingerid_params = fingerid.getCommandLine();

      const bool run_csifingerid = !out_csifingerid.empty();
      const bool run_passatutto = decoy_generation;

      // structure of the command line passed to NightSky (Sirius 4.X+)
      // Add noCite and instead make sure the main citations are registered in our TOPP tool.
      // Most of the time the user does not see the direct Sirius output anyway
      QStringList command_line = QStringList("--noCite") + project_params + QStringList({"--input", tmp_ms_file.toQString(), "--project", tmp_out_dir.toQString(), "--no-compression", "sirius"}) + sirius_params;

      if (run_passatutto)
      {
        command_line << "passatutto";
      }

      if (run_csifingerid)
      {
        command_line << "fingerprint" << "structure" << fingerid_params;
      }

      command_line << "write-summaries";

      // the actual process
      QProcess qp;
      QString executable_qstring = SiriusAdapterAlgorithm::determineSiriusExecutable(executable).toQString();
      qp.start(executable_qstring, command_line); // does automatic escaping etc... start
      std::stringstream ss;

      ss << "Executing command: " << executable_qstring.toStdString();
      for (const QString& it : command_line)
      {
        ss << " " << it.toStdString();
      }
      OPENMS_LOG_WARN << ss.str() << std::endl;

      const bool success = qp.waitForFinished(-1);

      if (!success || qp.exitStatus() != 0 || qp.exitCode() != 0)
      {
        OPENMS_LOG_WARN << "FATAL: External invocation of Sirius failed. Standard output and error were:" << std::endl;
        const QString sirius_stdout(qp.readAllStandardOutput());
        const QString sirius_stderr(qp.readAllStandardError());
        OPENMS_LOG_WARN << String(sirius_stdout) << std::endl;
        OPENMS_LOG_WARN << String(sirius_stderr) << std::endl;
        OPENMS_LOG_WARN << String(qp.exitCode()) << std::endl;
        qp.close();

        throw Exception::InvalidValue(__FILE__,
                                      __LINE__,
                                      OPENMS_PRETTY_FUNCTION,
                                      "FATAL: SIRIUS could not be executed!",
                                      "");
      }

      if (param_.getValue("read_sirius_stdout") == "true")
      {
        OPENMS_LOG_WARN << "Standard output and error of SIRIUS were:" << std::endl;
        const QString sirius_stdout(qp.readAllStandardOutput());
        const QString sirius_stderr(qp.readAllStandardError());
        OPENMS_LOG_WARN << String(sirius_stdout) << std::endl;
        OPENMS_LOG_WARN << String(sirius_stderr) << std::endl;
      }
      qp.close();

      //extract path to subfolders (sirius internal folder structure)
      std::vector<String> subdirs;
      QDirIterator it(tmp_out_dir.toQString(), QDir::Dirs | QDir::NoDotAndDotDot, QDirIterator::NoIteratorFlags);
      while (it.hasNext())
      {
        subdirs.emplace_back(it.next());
      }
      return subdirs;
    }

    // ################
    // Parameter handling
    // ################

    SiriusAdapterAlgorithm::ParameterModifier SiriusAdapterAlgorithm::ParameterSection::parameter(
            const String &parameter_name,
            const ParamValue &default_value,
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
