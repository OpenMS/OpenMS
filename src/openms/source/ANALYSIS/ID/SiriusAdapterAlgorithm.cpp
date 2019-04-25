// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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

namespace OpenMS
{
    SiriusAdapterAlgorithm::SiriusAdapterAlgorithm() :
      DefaultParamHandler("SiriusAdapterAlgorithm")
    {      
      // adapter parameters (preprocessing)
      defaults_.setValue("preprocessing:filter_by_num_masstraces", 1, "Features have to have at least x MassTraces. To use this parameter feature_only is neccessary");
      defaults_.setMinInt("preprocessing:filter_by_num_masstraces", 1);
      defaults_.setValue("preprocessing:precursor_mz_tolerance", 0.005, "Tolerance window for precursor selection (Feature selection in regard to the precursor)");
      defaults_.setValue("preprocessing:precursor_mz_tolerance_unit", "Da", "Unit of the precursor_mz_tolerance");
      defaults_.setValidStrings("preprocessing:precursor_mz_tolerance_unit", ListUtils::create<String>("Da,ppm"));
      defaults_.setValue("preprocessing:precursor_rt_tolerance", 5, "Tolerance window (left and right) for precursor selection [seconds]");
      defaults_.setValue("preprocessing:isotope_pattern_iterations", 3, "Number of iterations that should be performed to extract the C13 isotope pattern. If no peak is found (C13 distance) the function will abort. Be careful with noisy data - since this can lead to wrong isotope patterns.", ListUtils::create<String>("advanced"));
      // adapter flags
      defaults_.setValue("preprocessing:feature_only", "false", "Uses the feature information from in_featureinfo to reduce the search space to MS2 associated with a feature.");
      defaults_.setValidStrings("preprocessing:feature_only", ListUtils::create<String>("true,false"));
      defaults_.setValue("preprocessing:no_masstrace_info_isotope_pattern", "false", "Use this flag if the masstrace information from a feature should be discarded and the isotope_pattern_iterations should be used instead.", ListUtils::create<String>("advanced"));
      defaults_.setValidStrings("preprocessing:no_masstrace_info_isotope_pattern", ListUtils::create<String>("true,false"));
      defaults_.setSectionDescription("preprocessing", "Preprocessing");

      // parameters for SIRIUS (sirius)
      defaults_.setValue("sirius:profile", "qtof", "Specify the used analysis profile");
      defaults_.setValidStrings("sirius:profile", ListUtils::create<String>("qtof,orbitrap,fticr"));
      defaults_.setValue("sirius:candidates", 5, "The number of candidates in the SIRIUS output.");
      defaults_.setMinInt("sirius:candidates", 1);
      defaults_.setValue("sirius:database", "all", "search formulas in given database");
      defaults_.setValidStrings("sirius:database", ListUtils::create<String>("all,chebi,custom,kegg,bio,natural products,pubmed,hmdb,biocyc,hsdb,knapsack,biological,zinc bio,gnps,pubchem,mesh,maconda"));
      defaults_.setValue("sirius:noise", 0, "median intensity of noise peaks");
      defaults_.setMinInt("sirius:noise", 0);
      defaults_.setValue("sirius:ppm_max", 10, "allowed ppm for decomposing masses");
      defaults_.setValue("sirius:isotope", "both", "how to handle isotope pattern data. Use 'score' to use them for ranking or 'filter' if you just want to remove candidates with bad isotope pattern. With 'both' you can use isotopes for filtering and scoring. Use 'omit' to ignore isotope pattern.");
      defaults_.setValidStrings("sirius:isotope", ListUtils::create<String>("score,filter,both,omit"));
      defaults_.setValue("sirius:elements", "CHNOP[5]S[8]Cl[1]", "The allowed elements. Write CHNOPSCl to allow the elements C, H, N, O, P, S and Cl. Add numbers in brackets to restrict the maximal allowed occurrence of these elements: CHNOP[5]S[8]Cl[1].");
      defaults_.setValue("sirius:compound_timeout", 10, "Time out in seconds per compound. To disable the timeout set the value to 0");
      defaults_.setMinInt("sirius:compound_timeout", 0);
      defaults_.setValue("sirius:tree_timeout", 0, "Time out in seconds per fragmentation tree computation.");
      defaults_.setMinInt("sirius:tree_timeout", 0);
      defaults_.setValue("sirius:top_n_hits", 10, "The number of top hits for each compound written to the CSI:FingerID output");
      defaults_.setMinInt("sirius:top_n_hits", 1);
      defaults_.setValue("sirius:cores", 1, "The number of cores SIRIUS is allowed to use on the system");
      defaults_.setMinInt("sirius:cores", 1);
      // sirius flags
      defaults_.setValue("sirius:auto_charge", "false", "Use this option if the charge of your compounds is unknown and you do not want to assume [M+H]+ as default. With the auto charge option SIRIUS will not care about charges and allow arbitrary adducts for the precursor peak.");
      defaults_.setValidStrings("sirius:auto_charge", ListUtils::create<String>("true,false"));
      defaults_.setValue("sirius:ion_tree", "false", "Print molecular formulas and node labels with the ion formula instead of the neutral formula", ListUtils::create<String>("advanced"));
      defaults_.setValidStrings("sirius:ion_tree", ListUtils::create<String>("true,false"));
      defaults_.setValue("sirius:no_recalibration", "false", "If this option is set, SIRIUS will not recalibrate the spectrum during the analysis.", ListUtils::create<String>("advanced"));
      defaults_.setValidStrings("sirius:no_recalibration", ListUtils::create<String>("true,false"));
      defaults_.setValue("sirius:most_intense_ms2", "false", "SIRIUS uses the fragmentation spectrum with the most intense precursor peak (for each spectrum)", ListUtils::create<String>("advanced"));
      defaults_.setValidStrings("sirius:most_intense_ms2", ListUtils::create<String>("true,false"));
      defaults_.setSectionDescription("sirius", "Parameters for SIRIUS and CSI:FingerID");

      defaultsToParam_();
    }
    
    String SiriusAdapterAlgorithm::getFeatureOnly() { return feature_only_; }
    String SiriusAdapterAlgorithm::getNoMasstraceInfoIsotopePattern() { return no_masstrace_info_isotope_pattern_; }
    int SiriusAdapterAlgorithm::getIsotopePatternIterations() { return isotope_pattern_iterations_; }
    int SiriusAdapterAlgorithm::getCandidates() { return candidates_; }
    int SiriusAdapterAlgorithm::getTopNHits() { return top_n_hits_; }
     
    void SiriusAdapterAlgorithm::updateMembers_()
    {
      // adapter parameters (preprocessing)
      filter_by_num_masstraces_ = param_.getValue("preprocessing:filter_by_num_masstraces");
      precursor_mz_tolerance_ = param_.getValue("preprocessing:precursor_mz_tolerance");
      precursor_mz_tolerance_unit_ = param_.getValue("preprocessing:precursor_mz_tolerance_unit");
      precursor_rt_tolerance_ = param_.getValue("preprocessing:precursor_rt_tolerance");
      isotope_pattern_iterations_ = param_.getValue("preprocessing:isotope_pattern_iterations");
      // flags
      feature_only_ = param_.getValue("preprocessing:feature_only");
      no_masstrace_info_isotope_pattern_ = param_.getValue("preprocessing:no_masstrace_info_isotope_pattern");
      
      // parameters for SIRIUS (sirius)
      profile_ = param_.getValue("sirius:profile");
      candidates_ = param_.getValue("sirius:candidates");
      database_ = param_.getValue("sirius:database");
      noise_ = param_.getValue("sirius:noise");
      ppm_max_ = param_.getValue("sirius:ppm_max");
      isotope_ = param_.getValue("sirius:isotope");
      elements_ = param_.getValue("sirius:elements");
      compound_timeout_ = param_.getValue("sirius:compound_timeout");
      tree_timeout_ = param_.getValue("sirius:tree_timeout");
      top_n_hits_ = param_.getValue("sirius:top_n_hits");
      cores_ = param_.getValue("sirius:cores");
      // flags
      auto_charge_ = param_.getValue("sirius:auto_charge");
      ion_tree_ = param_.getValue("sirius:ion_tree");
      no_recalibration_ = param_.getValue("sirius:no_recalibration");
      most_intense_ms2_ = param_.getValue("sirius:most_intense_ms2");
    }   

    std::pair<String, String> SiriusAdapterAlgorithm::checkSiriusExecutablePath(String& executable)
    { 
      std::pair<String, String> executable_workdir;
      // if executable was not provided
      if (executable.empty())
      {
        const QProcessEnvironment env;
        const String& qsiriuspathenv = env.systemEnvironment().value("SIRIUS_PATH");
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

      // normalize file path
      QString exe = executable.toQString();
      QFileInfo file_info(exe);
      exe = file_info.canonicalFilePath();
  
      LOG_WARN << "Executable is: " + String(exe) << std::endl;
      const String path_to_executable = File::path(exe);
      executable_workdir = std::make_pair(exe.toStdString(), path_to_executable);
      
      return executable_workdir;
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
                                                     const SiriusAdapterAlgorithm& sirius_algo,
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
          
          bool feature_only = (sirius_algo.feature_only_ == "true") ? true : false;
          unsigned int num_masstrace_filter = sirius_algo.filter_by_num_masstraces_;
          double precursor_mz_tol = sirius_algo.precursor_mz_tolerance_;
          double precursor_rt_tol = sirius_algo.precursor_rt_tolerance_;
          bool ppm_prec = (sirius_algo.precursor_mz_tolerance_unit_ == "true") ? true : false; 
          
          if (num_masstrace_filter != 1 && !feature_only)
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
                                                                    ppm_prec);
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

    void SiriusAdapterAlgorithm::checkFeatureSpectraNumber(const String& featureinfo,
                                                           const FeatureMapping::FeatureToMs2Indices& feature_mapping, 
                                                           const MSExperiment& spectra,
                                                           const SiriusAdapterAlgorithm& sirius_algo)
    {
      bool feature_only = (sirius_algo.feature_only_ == "true") ? true : false;
      // number of features to be processed 
      if (feature_only && !featureinfo.empty())
      {
        LOG_WARN << "Number of features to be processed: " << feature_mapping.assignedMS2.size() << std::endl;
      }
      else if (!featureinfo.empty())
      {
        LOG_WARN << "Number of features to be processed: " << feature_mapping.assignedMS2.size() << std::endl;
        LOG_WARN << "Number of additional MS2 spectra to be processed: " << feature_mapping.unassignedMS2.size() << std::endl;
      } 
      else
      {
        int count_ms2 = 0;
        for (const auto& spec_it : spectra)
        {
          if (spec_it.getMSLevel() == 2)
          {
            count_ms2++;
          }
        }
        LOG_WARN << "Number of MS2 spectra to be processed: " << count_ms2 << std::endl;
      }
    } 

    // tmp_msfile (store), all parameters, out_dir (tmpstructure)
    const std::vector<String> SiriusAdapterAlgorithm::callSiriusQProcess(const String& tmp_ms_file,
                                                                         const String& tmp_out_dir,
                                                                         String& executable,
                                                                         const String& out_csifingerid,
                                                                         const SiriusAdapterAlgorithm& sirius_algo)
    {
      // assemble SIRIUS parameters
      QStringList process_params;
      process_params << "-p" << sirius_algo.profile_.toQString()
                     << "-e" << sirius_algo.elements_.toQString()
                     << "-d" << sirius_algo.database_.toQString()
                     << "-s" << sirius_algo.isotope_.toQString()
                     << "--noise" << QString::number(sirius_algo.noise_)
                     << "--candidates" << QString::number(sirius_algo.candidates_)
                     << "--ppm-max" << QString::number(sirius_algo.ppm_max_)
                     << "--compound-timeout" << QString::number(sirius_algo.compound_timeout_)
                     << "--tree-timeout" << QString::number(sirius_algo.tree_timeout_)
                     << "--processors" << QString::number(sirius_algo.cores_)
                     << "--quiet"
                     << "--output" << tmp_out_dir.toQString(); //internal output folder for temporary SIRIUS output file storage
  
      // add flags 
      if (sirius_algo.no_recalibration_ == "true")
      {
        process_params << "--no-recalibration";
      }
      if (!out_csifingerid.empty())
      {
        process_params << "--fingerid";
      }
      if (sirius_algo.ion_tree_ == "true")
      {
        process_params << "--iontree";
      }
      if (sirius_algo.auto_charge_ == "true")
      {
        process_params << "--auto-charge";
      }
      if (sirius_algo.most_intense_ms2_ == "true")
      {
        process_params << "--mostintense-ms2";
      }
  
      process_params << tmp_ms_file.toQString();
  
      // the actual process
      QProcess qp;
      std::pair<String, String> exe_wd = SiriusAdapterAlgorithm::checkSiriusExecutablePath(executable);
      QString exe = exe_wd.first.toQString();
      QString wd = exe_wd.second.toQString(); 
      qp.setWorkingDirectory(wd); //since library paths are relative to sirius executable path
      qp.start(exe, process_params); // does automatic escaping etc... start
      std::stringstream ss;
      ss << "COMMAND: " << String(exe);
      for (QStringList::const_iterator it = process_params.begin(); it != process_params.end(); ++it)
      {
          ss << " " << it->toStdString();
      }
      LOG_DEBUG << ss.str() << std::endl;
      LOG_WARN << "Executing: " + String(exe) << std::endl;
      LOG_WARN << "Working Dir is: " + String(wd) << std::endl;
      const bool success = qp.waitForFinished(-1); // wait till job is finished
  
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
      
      // extract path to subfolders (sirius internal folder structure)
      std::vector<String> subdirs;
      QDirIterator it(tmp_out_dir.toQString(), QDir::Dirs | QDir::NoDotAndDotDot, QDirIterator::NoIteratorFlags);
      while (it.hasNext())
      {
        subdirs.push_back(it.next());
      }
      return subdirs;
    }
} // namespace OpenMS

/// @endcond
