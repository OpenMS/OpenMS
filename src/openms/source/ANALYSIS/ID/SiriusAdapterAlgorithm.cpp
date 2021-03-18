// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Authors: Oliver Alka, Lukas Zimmermann $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/SiriusAdapterAlgorithm.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/DATAACCESS/SiriusMzTabWriter.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/SYSTEM/File.h>
#include <QDir>
#include <QDirIterator>
#include <QString>
#include <QtCore/QProcess>
#include <fstream>
#include <include/OpenMS/DATASTRUCTURES/StringUtils.h>

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
    using DefaultValue = DataValue;
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
          DefaultValue(),
          Description("Just consider compounds with a precursor mz lower or equal\n"
                      "this maximum mz. All other compounds in the input file\n"
                      "are ignored.")
      );

      parameter(
          ProjectName("processors"),
          DefaultValue(1),
          Description("Number of cpu cores to use. If not specified SIRIUS uses all available cores.")
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

    void SiriusAdapterAlgorithm::Sirius::parameters()
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
                 DefaultValue(0),
                 Description("Time out in seconds per fragmentation tree computations. 0 for an infinite amount of time")
               ).withMinInt(0);

      parameter(
                 SiriusName("compound-timeout"),
                 DefaultValue(100),
                 Description("Maximal computation time in seconds for a single compound. 0 for an infinite amount of time.")
               ).withMinInt(0);

      flag(
            SiriusName("no-recalibration"),
            Description("Disable recalibration of input spectra")
          );

      parameter(
                 SiriusName("profile"),
                 DefaultValue("qtof"),
                 Description("Name of the configuration profile")
               ).withValidStrings({"qtof", "orbitrap", "fticr"});

      parameter(
                 SiriusName("formula"),
                 DefaultValue(""),
                 Description("Specify the neutral molecular formula of the measured "
                             "compound to compute its tree or a list of candidate "
                             "formulas the method should discriminate. Omit this "
                             "option if you want to consider all possible molecular formulas")
               );

      parameter(
                 SiriusName("ions-enforced"),
                 DefaultValue(-1),
                 Description("the iontype/adduct of the MS/MS data. Example: [M+H]+, "
                             "[M-H]-, [M+Cl]-, [M+Na]+, [M]+. You can also provide a comma separated list of adducts")
               );

      parameter(
                 SiriusName("candidates"),
                 DefaultValue(5),
                 Description("The number of formula candidates in the SIRIUS output")
               ).withMinInt(1);

      parameter(
                 SiriusName("candidates-per-ion"),
                 DefaultValue(-1),
                 Description("Minimum number of candidates in the output for each "
                             "ionization. Set to force output of results for each "
                             "possible ionization, even if not part of highest "
                             "ranked results. -1 omits parameter in Sirius.")
                );

      parameter(
                 SiriusName("elements-considered"),
                 DefaultValue(""),
                 Description("Set the allowed elements for rare element detection. "
                             "Write SBrClBSe to allow the elements S,Br,Cl,B and Se."));

      parameter(
                 SiriusName("elements-enforced"),
                 DefaultValue(""),
                 Description("Enforce elements for molecular formula determination. "
                             "Write CHNOPSCl to allow the elements C, H, N, O, P, S "
                             "and Cl. Add numbers in brackets to restrict the "
                             "minimal and maximal allowed occurrence of these "
                             "elements: CHNOP[5]S[8]Cl[1-2]. When one number is "
                             "given then it is interpreted as upper bound. Default is CHNOP")
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
                 DefaultValue(""),
                 Description("the iontype/adduct of the MS/MS data. "
                             "Example: [M+H]+, [M-H]-, [M+Cl]-, [M+Na]+, [M]+. "
                             "You can also provide a comma separated list of adducts.")
                );

      parameter(
                 SiriusName("db"),
                 DefaultValue("all"),
                 Description("Search formulas in given database.")
               ).withValidStrings({"all",
                                   "pubchem",
                                   "mesh",
                                   "hmdb",
                                   "knapsack",
                                   "chebi",
                                   "pubmed",
                                   "bio",
                                   "kegg",
                                   "hsdb",
                                   "maconda",
                                   "metacyc",
                                   "gnps",
                                   "zincbio",
                                   "train",
                                   "undp",
                                   "pantcyc",
                                   "ymdb",
                                   "keggmine",
                                   "ecocycmine",
                                   "ymdbmine",
                                   "custom",
                                   "custom_0",
                                   "custom_1",
                                   "custom_2",
                                   "custom_3",
                                   "custom_4"});

      parameter(
                 SiriusName("ions-enforced"),
                 DefaultValue(""),
                 Description("The iontype/adduct of the MS/MS data. Example: [M+H]+, \n"
                             "[M-H]-, [M+Cl]-, [M+Na]+, [M]+. You can also provide a \n"
                             "comma separated list of adducts.")
               );
    }

    void SiriusAdapterAlgorithm::Fingerid::parameters()
    {
      parameter(
                 FingeridName("candidates"),
                 DefaultValue(10),
                 Description("Number of molecular structure candidates in the output.")
               ).withMinInt(1);

      parameter(
                 FingeridName("db"),
                 DefaultValue("all"),
                 Description("Search structure in given database. To use a custom database, please rename your database to custom_n (e.g. custom_0)")
                ).withValidStrings({"all",
                                    "pubchem",
                                    "mesh",
                                    "hmdb",
                                    "knapsack",
                                    "chebi",
                                    "pubmed",
                                    "bio",
                                    "kegg",
                                    "hsdb",
                                    "maconda",
                                    "metacyc",
                                    "gnps",
                                    "zincbio",
                                    "train",
                                    "undp",
                                    "pantcyc",
                                    "ymdb",
                                    "keggmine",
                                    "ecocycmine",
                                    "ymdbmine",
                                    "custom",
                                    "custom_0",
                                    "custom_1",
                                    "custom_2",
                                    "custom_3",
                                    "custom_4"});
    }

    void SiriusAdapterAlgorithm::Passatutto::parameters()
    {
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

    String SiriusAdapterAlgorithm::determineSiriusExecutable(String& executable)
    { 
      // if executable was not provided
      if (executable.empty())
      {
        const std::string& qsiriuspathenv(std::getenv("SIRIUS_PATH"));
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
      OPENMS_LOG_WARN << "Executable is: " + exe << std::endl;
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

      for (const auto& index : indices)
      {
        sorted_subdirs.emplace_back(std::move(subdirs[index.array_index]));
      }

      sorted_subdirs.swap(subdirs);
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
            OPENMS_LOG_WARN << "Parameter: filter_by_num_masstraces, was set to 1 to retain the adduct information for all MS2 spectra, if available. Please use the masstrace filter in combination with feature_only." << std::endl;
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

    void SiriusAdapterAlgorithm::logFeatureSpectraNumber(const String& featureinfo,
                                                         const FeatureMapping::FeatureToMs2Indices& feature_mapping,
                                                         const MSExperiment& spectra)
    {
      // number of features to be processed
      if (isFeatureOnly() && !featureinfo.empty())
      {
        OPENMS_LOG_WARN << "Number of features to be processed: " << feature_mapping.assignedMS2.size() << std::endl;
      }
      else if (!featureinfo.empty())
      {
        OPENMS_LOG_WARN << "Number of features to be processed: " << feature_mapping.assignedMS2.size() << std::endl;
        OPENMS_LOG_WARN << "Number of additional MS2 spectra to be processed: " << feature_mapping.unassignedMS2.size() << std::endl;
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

      // structure of the command line passed to NightSky
      QStringList command_line = project_params + QStringList({"--input", tmp_ms_file.toQString(), "--project", tmp_out_dir.toQString(), "sirius"}) + sirius_params;

      if (run_passatutto)
      {
        command_line << "passatutto";
      }

      if (run_csifingerid)
      {
        command_line << "fingerid" << fingerid_params;
      }

      OPENMS_LOG_INFO << "Running SIRIUS with the following command line parameters: " << endl;
      for (const auto &param: command_line)
      {
        OPENMS_LOG_INFO << param.toStdString() << " ";
      }
      OPENMS_LOG_INFO << endl;

      // the actual process
      QProcess qp;
      QString executable_qstring = SiriusAdapterAlgorithm::determineSiriusExecutable(executable).toQString();
      QString wd = File::path(executable).toQString();
      qp.setWorkingDirectory(wd); //since library paths are relative to sirius executable path
      //since library paths are relative to sirius executable path
      qp.setWorkingDirectory(File::path(executable).toQString());
      qp.start(executable_qstring, command_line); // does automatic escaping etc... start

      std::stringstream ss;
      ss << "COMMAND: " << executable_qstring.toStdString();
      for (QStringList::const_iterator it = command_line.begin(); it != command_line.end(); ++it)
      {
        ss << " " << it->toStdString();
      }
      OPENMS_LOG_WARN << ss.str() << std::endl;
      OPENMS_LOG_WARN << "Executing: " + executable_qstring.toStdString() << std::endl;

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

  // ################
  // Parameter handling
  // ################

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
