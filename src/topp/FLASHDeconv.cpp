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
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/MassFeatureTrace.h>
#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolutedSpectrum.h>
#include <OpenMS/FILTERING/DATAREDUCTION/MassTraceDetection.h>
#include <QDirIterator>
#include <QFileInfo>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/METADATA/SpectrumLookup.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------
/**
  @page TOPP_FLASHDeconv TOPP_FLASHDeconv
  (Need to be modified)

  @brief  @ref
  @code
  @endcode
  @verbinclude
  @htmlinclude
*/
// We do not want this class to show up in the docu:
// NEED to fill this part later


class TOPPFLASHDeconv :
    public TOPPBase
{
public:
  TOPPFLASHDeconv() :
      TOPPBase("FLASHDeconv", "Ultra-fast high-quality deconvolution enables online processing of top-down MS data",
               false)
  {
  }

protected:
  //typedef FLASHDeconvHelperStructs::Parameter Parameter;

  // this function will be used to register the tool parameters
  // it gets automatically called on tool execution
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<input_file/input_dir>", "", "Input file or directory");
    registerOutputFile_("out", "<output_file_prefix/output_dir>", "",
                        "output file prefix or output dir. "
                        "If prefix, [prefix].tsv (feature deconvolution results), [prefix]_MSn_spec.tsv (spectral deconvolution results) will be generated. "
                        "If dir, [dir]/[inputfile].tsv and [dir]/[inputfile]_MSn_spec.tsv will be generated per [inputfile].");

    registerDoubleList_("tol",
                        "<ms1_tol, ms2_tol, ...>",
                        {10.0, 10.0},
                        "ppm tolerance for MS1, 2, ...  "
                        "(e.g., -tol 10.0 15.0 to specify 10.0 and 15.0 ppm for MS1 and MS2, respectively)",
                        false);

    registerIntOption_("write_detail",
                       "<1:true 0:false>",
                       0,
                       "to write peak info per deconvoluted mass in detail or not in [prefix]_MSn_spec.tsv files. If set to 1, all peak information (m/z, intensity, charge, and isotope index) per mass is reported.",
                       false,
                       false);

    registerIntOption_("promex_out",
                       "",
                       0,
                       "if set, deconvoluted masses (for MS1 spectra) are reported in promex output format ([prefix]]_FD.ms1ft)",
                       false,
                       false);

    registerIntOption_("topfd_out",
                       "",
                       0,
                       "if set, deconvoluted masses are reported in topfd output format ([prefiX]_FD_msX.msalign, per MS level X)",
                       false,
                       false);

    registerIntOption_("mzml_out",
                       "",
                       0,
                       "if set, deconvoluted masses (for all spectra) are reported in mzml format ([preifx]_deconved.mzml)",
                       false,
                       false);

    registerIntOption_("max_MS_level", "", 2, "maximum MS level (inclusive) for deconvolution", false, true);

    registerIntOption_("use_ensamble_spectrum",
                       "",
                       0,
                       "if set to 1, all spectra will add up to a single ensamble spectrum (per MS level) for which the deconvolution is performed",
                       false,
                       false);

    registerIntOption_("use_RNA_averagine", "", 0, "if set to 1, RNA averagine model is used", false, true);

    Param fd_defaults = FLASHDeconvAlgorithm().getDefaults();
    // overwrite algorithm default so we export everything (important for copying back MSstats results)
    fd_defaults.setValue("min_charge", 1);
    fd_defaults.setValue("max_charge", 100);
    fd_defaults.setValue("min_mz", -1.0);
    fd_defaults.setValue("max_mz", -1.0);
    fd_defaults.setValue("min_RT", -1.0);
    fd_defaults.setValue("max_RT", -1.0);
    fd_defaults.setValue("min_mass", 50.0);
    fd_defaults.setValue("max_mass", 100000.0);
    fd_defaults.setValue("tol", DoubleList{10.0, 10.0}, "ppm tolerance, controlled by -tol option");
    fd_defaults.addTag("tol", "advanced"); // hide entry
    fd_defaults.setValue("min_peaks", IntList{3, 1});
    fd_defaults.addTag("min_peaks", "advanced");
    fd_defaults.setValue("min_intensity", .0, "intensity threshold");
    fd_defaults.addTag("min_intensity", "advanced");
    fd_defaults.setValue("min_isotope_cosine",
                         DoubleList{.75, .90},
                         "cosine threshold between avg. and observed isotope pattern for MS1, 2, ... (e.g., -min_isotope_cosine 0.8 0.6 to specify 0.8 and 0.6 for MS1 and MS2, respectively)");
    fd_defaults.addTag("min_isotope_cosine", "advanced");
    //fd_defaults.setValue("min_charge_cosine",
    //                     .75,
    //                     "cosine threshold between per-charge-intensity and fitted gaussian distribution (applies only to MS1)");
    //fd_defaults.addTag("min_charge_cosine", "advanced");
    fd_defaults.setValue("max_mass_count",
                         IntList{-1, -1},
                         "maximum mass count per spec for MS1, 2, ... (e.g., -max_mass_count 100 50 to specify 100 and 50 for MS1 and MS2, respectively. -1 specifies unlimited)");
    fd_defaults.addTag("max_mass_count", "advanced");
    fd_defaults.setValue("num_overlapped_scans", 10, "number of overlapped scans for MS1 deconvolution");
    fd_defaults.addTag("num_overlapped_scans", "advanced");

    Param mf_defaults = MassFeatureTrace().getDefaults();
    mf_defaults.setValue("mass_error_da",
                         1.5,
                         "da tolerance for feature tracing. Due to frequent isotope errer, 1.5 Da is recommended.");
    mf_defaults.addTag("mass_error_ppm", "advanced"); // hide entry
    mf_defaults.setValue("trace_termination_criterion", "outlier");
    mf_defaults.addTag("trace_termination_criterion", "advanced"); // hide entry
    mf_defaults.setValue("reestimate_mt_sd", "false", "");
    mf_defaults.addTag("reestimate_mt_sd", "advanced"); // hide entry
    mf_defaults.setValue("quant_method", "area", "");
    mf_defaults.addTag("quant_method", "advanced"); // hide entry
    mf_defaults.setValue("noise_threshold_int", .0, "");
    mf_defaults.addTag("noise_threshold_int", "advanced"); // hide entry
    mf_defaults.setValue("min_sample_rate", 0.01, "");
    mf_defaults.addTag("min_sample_rate", "advanced"); // hide entry
    mf_defaults.addTag("trace_termination_outliers", "advanced"); // hide entry
    //mf_defaults.addTag("min_trace_length", "advanced"); // hide entry
    //mf_defaults.setValue("trace_termination_outliers", numOverlappedScans, "");
    mf_defaults.setValue("min_trace_length", 10.0, "min feature trace length in second");//
    //mf_defaults.setValue("min_charge_cosine", .5, "controlled by -min_charge_cosine option");
    //mf_defaults.addTag("min_charge_cosine", "advanced");
    mf_defaults.setValue("min_isotope_cosine", .75, "controlled by -min_isotope_cosine option");
    mf_defaults.addTag("min_isotope_cosine", "advanced");

    Param combined;
    combined.insert("Algorithm:", fd_defaults);
    combined.insert("FeatureTracing:", mf_defaults);
    registerFullParam_(combined);
  }

  // the main_ function is called after all parameters are read
  ExitCodes main_(int, const char **) override
  {
    OPENMS_LOG_INFO << "Initializing ... " << endl;

    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String infilePath = getStringOption_("in");
    String outfilePath = getStringOption_("out");
    bool useRNAavg = getIntOption_("use_RNA_averagine") > 0;
    int maxMSLevel = getIntOption_("max_MS_level");
    bool ensamble = getIntOption_("use_ensamble_spectrum") > 0;
    bool outMzml = getIntOption_("mzml_out") > 0;
    bool outTopfd = getIntOption_("topfd_out") > 0;
    bool outPromex = getIntOption_("promex_out") > 0;
    bool writeDetail = getIntOption_("write_detail") > 0;
    //double rtWindow = getDoubleOption_("RT_window");
    //    double ms1tol = getDoubleList_("tol")[0];
    int currentMaxMSLevel = 0;

    //double maxMass = getDoubleOption_("max_mass"); // from FLASHDeconvAlgorithm

    // spectrum number per ms level per input file
    auto specCntr = new int[maxMSLevel];
    fill_n(specCntr, maxMSLevel, 0);
    // spectrum number with at least one deconvoluted mass per ms level per input file
    auto qspecCntr = new int[maxMSLevel];
    fill_n(qspecCntr, maxMSLevel, 0);
    // mass number per ms level per input file
    auto massCntr = new int[maxMSLevel];
    fill_n(massCntr, maxMSLevel, 0);
    // feature number per input file
    auto featureCntr = 0;

    // feature index written in the output file
    int featureIndex = 1;

    // total spectrum number per ms level
    auto total_specCntr = new int[maxMSLevel];
    fill_n(total_specCntr, maxMSLevel, 0);
    // total spectrum number with at least one deconvoluted mass per ms level
    auto total_qspecCntr = new int[maxMSLevel];
    fill_n(total_qspecCntr, maxMSLevel, 0);
    // total mass number per ms level
    auto total_massCntr = new int[maxMSLevel];
    // total feature number
    fill_n(total_massCntr, maxMSLevel, 0);
    auto total_featureCntr = 0;

    MSExperiment ensamble_map;
    // generate ensamble spectrum if param.ensamble is set
    if (ensamble)
    {
      for (auto i = 0; i < maxMSLevel; i++)
      {
        auto spec = MSSpectrum();
        spec.setMSLevel(i + 1);
        std::ostringstream name;
        name << "Ensamble MS" << (i + 1) << " spectrum";
        spec.setName(name.str());
        ensamble_map.addSpectrum(spec);
      }
    }

    double total_elapsed_cpu_secs = 0, total_elapsed_wall_secs = 0;
    fstream featureOut, promexOut, topfdOut_MS1, topfdOut_MS2;
    String mzmlOut;
    auto specOut = new fstream[maxMSLevel];

    //-------------------------------------------------------------
    // reading input file directory -> put that in array
    //-------------------------------------------------------------
    vector<QString> infileArray;
    QString path = QString::fromUtf8(infilePath.data(), (int) infilePath.size());
    QFileInfo check_file(path);
    if (check_file.isDir())
    {
      QDirIterator it(path, QStringList() << "*.mzml", QDir::Files, QDirIterator::Subdirectories);
      while (it.hasNext())
      {
        infileArray.push_back(it.next());
      }
    }
    else
    {
      infileArray.push_back(path);
    }
    sort(infileArray.begin(), infileArray.end());

    // check if output path is a directory
    bool isOutPathDir = (QFileInfo(QString::fromUtf8(outfilePath.data(), (int) outfilePath.size())).isDir());

    // if output path is a file name, output names should be specified accordingly here
    if (!isOutPathDir)
    {
      for (int n = 1; n <= (int) maxMSLevel; ++n)
      {
        specOut[n - 1]
            .open(outfilePath + "_MS" + n + (ensamble ? "_ensamble_spec" : "_spec") + ".tsv", fstream::out);
        DeconvolutedSpectrum::writeDeconvolutedMassesHeader(specOut[n - 1], n, writeDetail);
      }
      if (outMzml)
      {
        mzmlOut = outfilePath + "_deconved.mzml";//[preifx]_train_MSn.csv
      }

      if (!ensamble)
      {
        featureOut.open(outfilePath + ".tsv", fstream::out);
        MassFeatureTrace::writeHeader(featureOut);
      }
      if (outPromex)
      {
        promexOut.open(outfilePath + "_FD.ms1ft", fstream::out);
        MassFeatureTrace::writePromexHeader(promexOut);
      }

      if (outTopfd)
      {
        topfdOut_MS1.open(outfilePath + "_FD_ms1.msalign", fstream::out);
        topfdOut_MS2.open(outfilePath + "_FD_ms2.msalign", fstream::out);
        //writeTopFDHeader(topfdOut);
      }
    }


    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    for (auto &infile : infileArray)
    {
      // if output is a directory, initialze index and counters
      if (isOutPathDir)
      {
        fill_n(specCntr, maxMSLevel, 0);
        fill_n(qspecCntr, maxMSLevel, 0);
        fill_n(massCntr, maxMSLevel, 0);

        featureCntr = 0;
        featureIndex = 1;
      }
      MSExperiment map;
      MzMLFile mzml;

      // all for measure elapsed cup wall time
      double elapsed_cpu_secs = 0, elapsed_wall_secs = 0;
      auto elapsed_deconv_cpu_secs = new double[maxMSLevel];
      auto elapsed_deconv_wall_secs = new double[maxMSLevel];
      fill_n(elapsed_deconv_cpu_secs, maxMSLevel, 0);
      fill_n(elapsed_deconv_wall_secs, maxMSLevel, 0);

      auto begin = clock();
      auto t_start = chrono::high_resolution_clock::now();

      OPENMS_LOG_INFO << "Processing : " << infile.toStdString() << endl;

      mzml.setLogType(log_type_);
      mzml.load(infile, map);
      auto fileName = QFileInfo(infile).fileName().toStdString();

      //      double rtDuration = map[map.size() - 1].getRT() - map[0].getRT();
      auto ms1Cntr = 0;
      auto ms2Cntr = .0; // for debug...
      currentMaxMSLevel = 0;

      // read input dataset once to count spectra and generate ensamble spectrum if necessary
      for (auto &it : map)
      {
        if (it.empty())
        {
          continue;
        }
        if (it.getMSLevel() > maxMSLevel)
        {
          continue;
        }

        auto msLevel = it.getMSLevel();
        currentMaxMSLevel = currentMaxMSLevel < msLevel ? msLevel : currentMaxMSLevel;

        if (ensamble)
        {
          auto &espec = ensamble_map[it.getMSLevel() - 1];
          for (auto &p : it)
          {
            espec.push_back(p);
          }
        }

        if (it.getMSLevel() == 1)
        {
          ms1Cntr++;
        }
        if (it.getMSLevel() == 2)
        {
          ms2Cntr++;
        }
      }

      // Max MS Level is adjusted according to the input dataset
      currentMaxMSLevel = currentMaxMSLevel > maxMSLevel ? maxMSLevel : currentMaxMSLevel;
      // if an ensamble spectrum is analyzed, replace the input dataset with the ensamble one
      if (ensamble)
      {
        for (int i = 0; i < currentMaxMSLevel; ++i)
        {
          ensamble_map[i].sortByPosition();
        }
        map = ensamble_map;
      }


      std::string outfileName(fileName);

      // is output is dir, the output file name should be specified by the input file name here
      if (isOutPathDir)
      {
        std::size_t found = outfileName.find_last_of('.');
        outfileName = outfileName.substr(0, found);

        for (int n = 1; n <= (int) maxMSLevel; ++n)
        {
          specOut[n - 1]
              .open(outfilePath + outfileName + "_MS" + n + (ensamble ? "_ensamble_spec" : "_spec") + ".tsv",
                    fstream::out);
          DeconvolutedSpectrum::writeDeconvolutedMassesHeader(specOut[n - 1], n, writeDetail);
        }
        if (outMzml)
        {
          mzmlOut = outfilePath + outfileName + "_deconved.mzml";//[preifx]_train_MSn.csv
        }

        if (!ensamble)
        {
          featureOut.open(outfilePath + outfileName + ".tsv", fstream::out);
          MassFeatureTrace::writeHeader(featureOut);
        }

        if (outPromex)
        {
          promexOut.open(outfilePath + outfileName + "_FD.ms1ft", fstream::out);
          MassFeatureTrace::writePromexHeader(promexOut);
        }

        if (outTopfd)
        {
          topfdOut_MS1.open(outfilePath + outfileName + "_FD_ms1.msalign", fstream::out);
          topfdOut_MS2.open(outfilePath + outfileName + "_FD_ms2.msalign", fstream::out);
        }
      }

      // finally run FLASHDeconv here

      int scanNumber = 0;
      float prevProgress = .0;
      auto lastDeconvolutedSpectra = std::unordered_map<UInt, DeconvolutedSpectrum>();

      MSExperiment exp;


      auto fd = FLASHDeconvAlgorithm();
      Param fd_param = getParam_().copy("Algorithm:", true);
      fd_param.setValue("tol", getParam_().getValue("tol"));
      fd.setParameters(fd_param);
      fd.calculateAveragine(useRNAavg);
      auto avg = fd.getAveragine();
      auto massTracer = MassFeatureTrace();
      Param mf_param = getParam_().copy("FeatureTracing:", true);
      DoubleList isotopeCosine = fd_param.getValue("min_isotope_cosine");
      //mf_param.setValue("mass_error_ppm", ms1tol);
      mf_param.setValue("trace_termination_outliers", fd_param.getValue("num_overlapped_scans"));
      //mf_param.setValue("min_charge_cosine", fd_param.getValue("min_charge_cosine"));
      mf_param.setValue("min_isotope_cosine", isotopeCosine[0]);

      massTracer.setParameters(mf_param);

      OPENMS_LOG_INFO << "Running FLASHDeconv ... " << endl;

      for (auto it = map.begin(); it != map.end(); ++it)
      {
        scanNumber = SpectrumLookup::extractScanNumber(it->getNativeID(),
                                                       map.getSourceFiles()[0].getNativeIDTypeAccession());
        if (it->empty())
        {
          continue;
        }

        auto msLevel = it->getMSLevel();
        if (msLevel > currentMaxMSLevel)
        {
          continue;
        }

        specCntr[msLevel - 1]++;
        auto deconv_begin = clock();
        auto deconv_t_start = chrono::high_resolution_clock::now();

        auto deconvolutedSpectrum = DeconvolutedSpectrum(*it, scanNumber);
        // for MS>1 spectrum, register precursor
        if (msLevel > 1 && lastDeconvolutedSpectra.find(msLevel - 1) != lastDeconvolutedSpectra.end())
        {
          deconvolutedSpectrum.registerPrecursor(lastDeconvolutedSpectra[msLevel - 1]);
        }
        // per spec deconvolution
        fd.fillPeakGroupsInDeconvolutedSpectrum(deconvolutedSpectrum, scanNumber);

        if (outMzml)
        {
          exp.addSpectrum(deconvolutedSpectrum.toSpectrum());
        }
        elapsed_deconv_cpu_secs[msLevel - 1] += double(clock() - deconv_begin) / CLOCKS_PER_SEC;
        elapsed_deconv_wall_secs[msLevel - 1] += chrono::duration<double>(
            chrono::high_resolution_clock::now() - deconv_t_start).count();

        if (msLevel < currentMaxMSLevel)
        {
          lastDeconvolutedSpectra[msLevel] = deconvolutedSpectrum; // to register precursor in the future..
        }
        if (deconvolutedSpectrum.empty())
        {
          continue;
        }
        if (!ensamble)
        {
          massTracer.addDeconvolutedSpectrum(deconvolutedSpectrum);// add deconvoluted mass in massTracer
        }
        qspecCntr[msLevel - 1]++;
        massCntr[msLevel - 1] += deconvolutedSpectrum.size();
        deconvolutedSpectrum
            .writeDeconvolutedMasses(specOut[msLevel - 1], fileName, avg, writeDetail);

        if (outTopfd)
        {
          deconvolutedSpectrum.writeTopFD(msLevel == 1 ? topfdOut_MS1 : topfdOut_MS2, scanNumber);
        }

        deconvolutedSpectrum.clearPeakGroupsChargeInfo();
        float progress = (float) (it - map.begin()) / map.size();
        if (progress > prevProgress + .01)
        {
          printProgress(progress);
          prevProgress = progress;
        }
      }

      printProgress(1); //
      std::cout << std::endl;
      std::unordered_map<UInt, DeconvolutedSpectrum>().swap(lastDeconvolutedSpectra); // empty memory

      // massTracer run
      if (!ensamble)
      {
        massTracer
            .findFeatures(fileName, outPromex, featureCntr, featureIndex, featureOut, promexOut, fd.getAveragine());
      }
      if (outMzml)
      {
        MzMLFile mzMlFile;
        mzMlFile.store(mzmlOut, exp);
      }

      // output messages, etc.
      if (isOutPathDir)
      {
        for (int j = 0; j < (int) maxMSLevel; j++)
        {
          if (specCntr[j] == 0)
          {
            continue;
          }
          if (ensamble)
          {
            OPENMS_LOG_INFO << "In this run, FLASHDeconv found " << massCntr[j] << " masses in the ensamble MS"
                            << (j + 1) << " spectrum" << endl;

          }
          else
          {
            OPENMS_LOG_INFO << "In this run, FLASHDeconv found " << massCntr[j] << " masses in " << qspecCntr[j]
                            << " MS" << (j + 1) << " spectra out of "
                            << specCntr[j] << endl;
          }
        }

        if (featureCntr > 0)
        {
          OPENMS_LOG_INFO << "Mass tracer found " << featureCntr << " features" << endl;
        }

        for (int n = 1; n <= (int) maxMSLevel; ++n)
        {
          specOut[n - 1].close();

          if (specCntr[n - 1] == 0)// remove empty files
          {
            QString filename = QString::fromStdString(
                outfilePath + outfileName + "_MS" + n + (ensamble ? "_ensamble_spec" : "_spec") + ".tsv");
            QFile file(filename);
            file.remove();
          }
        }

        if (!ensamble)
        {
          featureOut.close();
        }
        if (outPromex)
        {
          promexOut.close();
        }

        if (outTopfd)
        {
          topfdOut_MS1.close();
          topfdOut_MS2.close();
        }

        for (int j = 0; j < (int) maxMSLevel; j++)
        {
          total_specCntr[j] += specCntr[j];
          total_qspecCntr[j] += qspecCntr[j];
          total_massCntr[j] += massCntr[j];
        }
        total_featureCntr += featureCntr;
      }
      else
      {
        for (int j = 0; j < (int) currentMaxMSLevel; j++)
        {
          if (specCntr[j] == 0)
          {
            continue;
          }

          if (ensamble)
          {
            OPENMS_LOG_INFO << "So far, FLASHDeconv found " << massCntr[j] << " masses in the ensamble MS"
                            << (j + 1) << " spectrum" << endl;

          }
          else
          {
            OPENMS_LOG_INFO << "So far, FLASHDeconv found " << massCntr[j] << " masses in " << qspecCntr[j]
                            << " MS" << (j + 1) << " spectra out of "
                            << specCntr[j] << endl;
          }
        }
        if (featureCntr > 0)
        {
          OPENMS_LOG_INFO << "Mass tracer found " << featureCntr << " features" << endl;
        }
        for (int j = 0; j < (int) maxMSLevel; j++)
        {
          total_specCntr[j] = specCntr[j];
          total_qspecCntr[j] = qspecCntr[j];
          total_massCntr[j] = massCntr[j];
        }
        total_featureCntr = featureCntr;
      }

      auto t_end = chrono::high_resolution_clock::now();
      auto end = clock();

      elapsed_cpu_secs = double(end - begin) / CLOCKS_PER_SEC;
      elapsed_wall_secs = chrono::duration<double>(t_end - t_start).count();

      OPENMS_LOG_INFO << "-- done [took " << elapsed_cpu_secs << " s (CPU), " << elapsed_wall_secs
                      << " s (Wall)] --"
                      << endl;

      auto sumCntr = 0;
      for (int j = 0; j < (int) currentMaxMSLevel; j++)
      {
        sumCntr += specCntr[j];

        OPENMS_LOG_INFO << "-- deconv per MS" << (j + 1) << " spectrum (except spec loading, feature finding) [took "
                        << 1000.0 * elapsed_deconv_cpu_secs[j] / sumCntr
                        << " ms (CPU), " << 1000.0 * elapsed_deconv_wall_secs[j] / sumCntr << " ms (Wall)] --" << endl;
      }
      delete[] elapsed_deconv_cpu_secs;
      delete[] elapsed_deconv_wall_secs;

      total_elapsed_cpu_secs += elapsed_cpu_secs;
      total_elapsed_wall_secs += elapsed_wall_secs;
    }


    OPENMS_LOG_INFO << "Total elapsed time\n-- done [took " << total_elapsed_cpu_secs << " s (CPU), "
                    << total_elapsed_wall_secs
                    << " s (Wall)] --"
                    << endl;

    if (massCntr[0] < total_massCntr[0])
    {
      for (int j = 0; j < (int) currentMaxMSLevel; j++)
      {
        if (total_specCntr[j] == 0)
        {
          continue;
        }

        if (ensamble)
        {
          OPENMS_LOG_INFO << "In total, FLASHDeconv found " << total_massCntr[j] << " masses in the ensamble MS"
                          << (j + 1) << " spectrum" << endl;
        }
        else
        {
          OPENMS_LOG_INFO << "In total, FLASHDeconv found " << total_massCntr[j] << " masses in " << total_qspecCntr[j]
                          << " MS" << (j + 1) << " spectra out of "
                          << total_specCntr[j] << endl;
        }
      }

      if (total_featureCntr > 0)
      {
        OPENMS_LOG_INFO << "Mass tracer found " << total_featureCntr << " features" << endl;
      }
    }

    if (!isOutPathDir)
    {
      for (int n = 1; n <= (int) maxMSLevel; ++n)
      {
        specOut[n - 1].close();

        if (qspecCntr[n - 1] == 0)
        {
          QString filename = QString::fromStdString(
              outfilePath + "_MS" + n + (ensamble ? "_ensamble_spec" : "_spec") + ".tsv");
          QFile file(filename);
          file.remove();
        }
      }

      if (!ensamble)
      {
        featureOut.close();
      }
      if (outPromex)
      {
        promexOut.close();
      }

      if (outTopfd)
      {
        topfdOut_MS1.close();
        topfdOut_MS2.close();
      }
    }
    return EXECUTION_OK;
  }

  static void printProgress(float progress)
  {
    float barWidth = 70;
    std::cout << "[";
    int pos = (int) (barWidth * progress);
    for (auto i = 0; i < barWidth; ++i)
    {
      if (i < pos)
      {
        std::cout << "=";
      }
      else if (i == pos)
      {
        std::cout << ">";
      }
      else
      {
        std::cout << " ";
      }
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
  }
};


// the actual main function needed to create an executable
int main(int argc, const char **argv)
{
  TOPPFLASHDeconv tool;
  return tool.main(argc, argv);
}