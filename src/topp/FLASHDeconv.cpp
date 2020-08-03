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
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/FILTERING/DATAREDUCTION/MassTraceDetection.h>
#include <QDirIterator>
#include <QFileInfo>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/MzMLFile.h>

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
  //typedef FLASHDeconvHelperStructs::PeakGroup PeakGroup;
  typedef FLASHDeconvHelperStructs::Parameter Parameter;
  //typedef FLASHDeconvHelperStructs::DeconvolutedSpectrum DeconvolutedSpectrum;

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
                        {10.0, 5.0},
                        "ppm tolerance for MS1, 2, ...  "
                        "(e.g., -tol 10.0 5.0 to specify 10.0 and 5.0 ppm for MS1 and MS2, respectively)",
                        false);

    //registerDoubleOption_("tol", "<tolerance>", 10.0, "ppm tolerance", false, false);
    registerIntOption_("min_charge", "<min_charge>", 1, "minimum charge state", false, false);
    registerIntOption_("max_charge", "<max_charge>", 100, "maximum charge state", false, false);
    registerDoubleOption_("min_mass", "<min_mass>", 50.0, "minimum mass (Da)", false, false);
    registerDoubleOption_("max_mass", "<max_mass>", 100000.0, "maximum mass (Da)", false, false);

    registerDoubleList_("min_isotope_cosine",
                        "<ms1_isotope_cos ms2_isotpe_cos, ...>",
                        {.75, .75},
                        "cosine threshold between avg. and observed isotope pattern for MS1, 2, ... "
                        "(e.g., -min_isotope_cosine 0.8 0.6 to specify 0.8 and 0.6 for MS1 and MS2, respectively)",
                        false,
                        true);

    registerDoubleOption_("min_charge_cosine",
                          "<charge_cosine>",
                          .5,
                          "cosine threshold between per-charge-intensity and fitted gaussian distribution (applies only to MS1)",
                          false,
                          true);

    registerIntList_("min_peaks",
                     "<ms1_min_supporting_peaks ms2_min_supproting_peaks, ...>",
                     {3, 1},
                     "minimum number of supporting peaks for MS1, 2, ... "
                     "(e.g., -min_peaks 3 2 to specify 3 and 2 for MS1 and MS2, respectively",
                     false,
                     true);

    registerIntList_("max_mass_count",
                     "<ms1_max_mass_count, ms2_max_mass_count, ...>",
                     {-1, -1},
                     "maximum mass count per spec for MS1, 2, ..."
                     "(e.g., -max_mass_count 100 50 to specify 100 and 50 for MS1 and MS2, respectively. -1 specifies unlimited)",
                     false,
                     true);
    //
    registerDoubleOption_("min_intensity", "<min_intensity>", 0, "intensity threshold (default 0.0)", false, true);
    registerDoubleOption_("RT_window",
                          "<seconds>",
                          0.0,
                          "RT window duration in seconds. (if 0, RT window contains 15 MS1 spectra)",
                          false,
                          true);

    registerDoubleOption_("min_RT_span", "<seconds>", 10.0, "Min feature RT span", false, true);
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
                       true);
    registerIntOption_("topfd_out",
                       "",
                       0,
                       "if set, deconvoluted masses (for MS2 spectra) are reported in topfd output format ([prefix]_FD_ms2.msalign)",
                       false,
                       true);
    registerIntOption_("mzml_out",
                       "",
                       0,
                       "if set, deconvoluted masses (for all spectra) are reported in mzml format ([preifx].mzml)",
                       false,
                       true);

    registerIntOption_("max_MS_level", "", 2, "maximum MS level (inclusive) for deconvolution", false, true);
  }

  Parameter setParameter()
  {
    Parameter param;
    param.minCharge = getIntOption_("min_charge");
    int maxCharge = getIntOption_("max_charge");

    if (maxCharge < 0 && param.minCharge < 0)
    {
      param.minCharge = -param.minCharge;
      maxCharge = -maxCharge;
      param.chargeMass = Constants::ELECTRON_MASS_U;
    }
    else
    {
      param.chargeMass = Constants::PROTON_MASS_U;
    }

    if (param.minCharge > maxCharge)
    {
      int tmp = param.minCharge;
      param.minCharge = maxCharge;
      maxCharge = tmp;
    }

    param.currentChargeRange = param.chargeRange = maxCharge - param.minCharge + 1;
    param.currentMaxMass = param.maxMass = getDoubleOption_("max_mass");
    param.minMass = getDoubleOption_("min_mass");
    param.tolerance = getDoubleList_("tol");

    for (auto j = 0; j < (int) param.tolerance.size(); j++)
    {
      param.tolerance[j] *= 1e-6;
      param.binWidth.push_back(.5 / param.tolerance[j]);
    }

    param.intensityThreshold = getDoubleOption_("min_intensity");
    param.minContinuousChargePeakCount = getIntList_("min_peaks");
    param.minIsotopeCosine = getDoubleList_("min_isotope_cosine");
    param.minChargeCosine = getDoubleOption_("min_charge_cosine");

    param.maxMassCount = getIntList_("max_mass_count");

    //param.chargeDistributionScoreThreshold = getDoubleOption_("minCDScore");
    param.RTwindow = getDoubleOption_("RT_window");
    param.minRTSpan = getDoubleOption_("min_RT_span");
    //param.threads = getIntOption_("threads");
    param.writeDetail = getIntOption_("write_detail");
    //param.jitter = getIntOption_("jitter");
    param.maxMSLevel = getIntOption_("max_MS_level");
    param.promexOut = getIntOption_("promex_out");
    param.topfdOut = getIntOption_("topfd_out");
    param.mzmlOut = getIntOption_("mzml_out");
    return param;
  }

  // the main_ function is called after all parameters are read
  ExitCodes main_(int, const char **) override
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String infilePath = getStringOption_("in");
    String outfilePath = getStringOption_("out");

    auto param = setParameter();
    auto avgine = FLASHDeconvHelperStructs::calculateAveragines(param);

    auto specCntr = new int[param.maxMSLevel];
    fill_n(specCntr, param.maxMSLevel, 0);
    auto qspecCntr = new int[param.maxMSLevel];
    fill_n(qspecCntr, param.maxMSLevel, 0);
    auto massCntr = new int[param.maxMSLevel];
    fill_n(massCntr, param.maxMSLevel, 0);
    auto featureCntr = 0;

    auto total_specCntr = new int[param.maxMSLevel];
    fill_n(total_specCntr, param.maxMSLevel, 0);
    auto total_qspecCntr = new int[param.maxMSLevel];
    fill_n(total_qspecCntr, param.maxMSLevel, 0);
    auto total_massCntr = new int[param.maxMSLevel];
    fill_n(total_massCntr, param.maxMSLevel, 0);
    auto total_featureCntr = 0;

    int specIndex = 1;
    int massIndex = 1;
    int featureIndex = 1;

    double total_elapsed_cpu_secs = 0, total_elapsed_wall_secs = 0;
    fstream featureOut, promexOut, topfdOut;
    String mzmlOut;
    auto specOut = new fstream[param.maxMSLevel];
    //auto mzmlOut = new fstream[param.maxMSLevel];
    //fstream fa;
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

    bool isOutPathDir = (QFileInfo(QString::fromUtf8(outfilePath.data(), (int) outfilePath.size())).isDir());

    OPENMS_LOG_INFO << "Initializing ... " << endl;

    if (!isOutPathDir)
    {
      for (int n = 1; n <= (int) param.maxMSLevel; ++n)
      {
        specOut[n - 1].open(outfilePath + "_MS" + n + "_spec.tsv", fstream::out);
        DeconvolutedSpectrum::writeDeconvolutedMassesHeader(specOut[n - 1], n, param.writeDetail);
      }
      if (param.mzmlOut)
      {
        mzmlOut = outfilePath + ".mzml";//[preifx]_train_MSn.csv
        //DeconvolutedSpectrum::writeAttCsvHeader(mzmlOut[n - 1]); //
      }
      //fa.open(outfilePath + "_MassList.csv", fstream::out);//[preifx]_train_MSn.csv
      //DeconvolutedSpectrum::writeAttCsvHeader(fa);
      //fsm.open(outfilePath + ".csv", fstream::out);
      //writeAttCsvHeader(fsm);

      featureOut.open(outfilePath + ".tsv", fstream::out);
      MassFeatureTrace::writeHeader(featureOut);

      if (param.promexOut)
      {
        promexOut.open(outfilePath + "_FD.ms1ft", fstream::out);
        MassFeatureTrace::writePromexHeader(promexOut);
      }

      if (param.topfdOut)
      {
        topfdOut.open(outfilePath + "_FD_ms2.msalign", fstream::out);
        //writeTopFDHeader(topfdOut);
      }
    }

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    for (auto &infile : infileArray)
    {
      if (isOutPathDir)
      {
        fill_n(specCntr, param.maxMSLevel, 0);
        fill_n(qspecCntr, param.maxMSLevel, 0);
        fill_n(massCntr, param.maxMSLevel, 0);

        featureCntr = 0;
      }
      MSExperiment map;
      MzMLFile mzml;

      double elapsed_cpu_secs = 0, elapsed_wall_secs = 0;
      //double elapsed_deconv_cpu_secs = 0, elapsed_deconv_wall_secs = 0;

      auto elapsed_deconv_cpu_secs = new double[param.maxMSLevel];
      auto elapsed_deconv_wall_secs = new double[param.maxMSLevel];
      fill_n(elapsed_deconv_cpu_secs, param.maxMSLevel, 0);
      fill_n(elapsed_deconv_wall_secs, param.maxMSLevel, 0);


      auto begin = clock();
      auto t_start = chrono::high_resolution_clock::now();

      OPENMS_LOG_INFO << "Processing : " << infile.toStdString() << endl;

      mzml.setLogType(log_type_);
      mzml.load(infile, map);

      param.fileName = QFileInfo(infile).fileName().toStdString();

      double rtDuration = map[map.size() - 1].getRT() - map[0].getRT();
      auto ms1Cntr = 0;
      auto ms2Cntr = .0; // for debug...

      for (auto &it : map)
      {
        if (it.getMSLevel() > param.maxMSLevel)
        {
          continue;
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

      //auto numMS2perMS1 = round(ms2Cntr/ms1Cntr)

      double rtDelta = rtDuration / ms1Cntr;

      auto rw = param.RTwindow;
      param.numOverlappedScans = max(param.minNumOverLappedScans, (UInt) round(rw / rtDelta));
      OPENMS_LOG_INFO << "# Overlapped MS1 scans:" << param.numOverlappedScans << " (in RT " << (rtDelta * param.numOverlappedScans) //
          << " sec)" << endl;

      std::string outfileName(param.fileName);

      if (isOutPathDir)
      {
        std::size_t found = outfileName.find_last_of('.');
        outfileName = outfileName.substr(0, found);

        for (int n = 1; n <= (int) param.maxMSLevel; ++n)
        {
          specOut[n - 1].open(outfilePath + outfileName + "_MS" + n + "_spec.tsv", fstream::out);
          DeconvolutedSpectrum::writeDeconvolutedMassesHeader(specOut[n - 1], n, param.writeDetail);
        }
        if (param.mzmlOut)
        {
          mzmlOut = outfilePath + outfileName + ".mzml";//[preifx]_train_MSn.csv
          //DeconvolutedSpectrum::writeAttCsvHeader(mzmlOut[n - 1]);
        }

        //mzmlOut.open(outfilePath + outfileName + "_MassList.csv", fstream::out);//[preifx]_train_MSn.csv
        //DeconvolutedSpectrum::writeThermoInclusionHeader(fa);

        featureOut.open(outfilePath + outfileName + ".tsv", fstream::out);
        MassFeatureTrace::writeHeader(featureOut);

        if (param.promexOut)
        {
          promexOut.open(outfilePath + outfileName + "_FD.ms1ft", fstream::out);
          MassFeatureTrace::writePromexHeader(promexOut);
        }

        if (param.topfdOut)
        {
          topfdOut.open(outfilePath + outfileName + "_FD_ms2.msalign", fstream::out);
        }
      }

      //check max ms level from the input dataset..
      param.currentMaxMSLevel = 0;
      //param.print();
      for (auto &it : map)
      {
        auto msLevel = it.getMSLevel();
        param.currentMaxMSLevel = param.currentMaxMSLevel < msLevel ? msLevel : param.currentMaxMSLevel;
      }

      param.currentMaxMSLevel = param.currentMaxMSLevel > param.maxMSLevel ? param.maxMSLevel : param.currentMaxMSLevel;

      Param common_param = getParam_().copy("algorithm:common:", true);
      writeDebug_("Common parameters passed to sub-algorithms (mtd and ffm)", common_param, 3);

      Param mtd_param = getParam_().copy("algorithm:mtd:", true);
      writeDebug_("Parameters passed to MassTraceDetection", mtd_param, 3);

      mtd_param.insert("", common_param);
      mtd_param.remove("chrom_fwhm");

      OPENMS_LOG_INFO << "Running FLASHDeconv ... " << endl;

      int scanNumber = 0;
      float prevProgress = .0;

      auto massTracer = MassFeatureTrace(param, mtd_param, avgine);
      auto lastDeconvolutedSpectra = std::unordered_map<UInt, DeconvolutedSpectrum>();

      MSExperiment exp;
      auto fd = FLASHDeconvAlgorithm(avgine, param);

      for (auto it = map.begin(); it != map.end(); ++it)
      {
        ++scanNumber;
        auto msLevel = it->getMSLevel();
        if (msLevel > param.currentMaxMSLevel)
        {
          continue;
        }

        specCntr[msLevel - 1]++;
        auto deconv_begin = clock();
        auto deconv_t_start = chrono::high_resolution_clock::now();

        // per spec deconvolution
        auto deconvolutedSpectrum = DeconvolutedSpectrum(*it, scanNumber);

        //bool proceed = true;
        param.currentChargeRange = param.chargeRange;
        param.currentMaxMass = param.maxMass;
        if (msLevel > 1 && lastDeconvolutedSpectra.find(msLevel - 1) != lastDeconvolutedSpectra.end())
        {
          deconvolutedSpectrum.registerPrecursor(lastDeconvolutedSpectra[msLevel - 1]);
        }
        //if (proceed)
        //{
        fd.getPeakGroups(deconvolutedSpectrum, specIndex, massIndex);
        //}
        if(param.mzmlOut){
          exp.addSpectrum(deconvolutedSpectrum.toSpectrum());
        }
        elapsed_deconv_cpu_secs[msLevel - 1] += double(clock() - deconv_begin) / CLOCKS_PER_SEC;
        elapsed_deconv_wall_secs[msLevel - 1] += chrono::duration<double>(
            chrono::high_resolution_clock::now() - deconv_t_start).count();

        if(msLevel<param.currentMaxMSLevel){
          lastDeconvolutedSpectra[msLevel] = deconvolutedSpectrum;
        }
        if (deconvolutedSpectrum.empty())
        {
          continue;
        }

        massTracer.addDeconvolutedSpectrum(deconvolutedSpectrum);
        qspecCntr[msLevel - 1]++;
        massCntr[msLevel - 1] += deconvolutedSpectrum.peakGroups.size();
        deconvolutedSpectrum.writeDeconvolutedMasses(specOut[msLevel - 1], param);

        if (param.topfdOut)
        {
          deconvolutedSpectrum.writeTopFD(topfdOut, specIndex - 1);
        }

        deconvolutedSpectrum.clearPeakGroupsChargeInfo();
        float progress = (float) (it - map.begin()) / map.size();
        if (progress > prevProgress + .01)
        {
          printProgress(progress); //
          prevProgress = progress;
        }
      }
      printProgress(1); //
      std::cout << std::endl;
      std::unordered_map<UInt, DeconvolutedSpectrum>().swap(lastDeconvolutedSpectra);

      massTracer.findFeatures(featureCntr, featureIndex, featureOut, promexOut);

      if(param.mzmlOut)
      {
        MzMLFile mzMlFile;
        mzMlFile.store(mzmlOut, exp);
      }

      if (isOutPathDir)
      {
        for (int j = 0; j < (int) param.maxMSLevel; j++)
        {
          if (specCntr[j] == 0)
          {
            continue;
          }
          OPENMS_LOG_INFO << "In this run, FLASHDeconv found " << massCntr[j] << " masses in " << qspecCntr[j]
                          << " MS" << (j + 1) << " spectra out of "
                          << specCntr[j] << endl;
        }

        if (featureCntr > 0)
        {
          OPENMS_LOG_INFO << "Mass tracer found " << featureCntr << " features" << endl;
        }

        for (int n = 1; n <= (int) param.maxMSLevel; ++n)
        {
          specOut[n - 1].close();

          if (specCntr[n - 1] == 0)
          {
            QString filename = QString::fromStdString(outfilePath + outfileName + "_MS" + n + "_spec.tsv");
            QFile file(filename);
            file.remove();
          }
        }
        //fa.close();
        // fsm.close();

        featureOut.close();
        if (param.promexOut)
        {
          promexOut.close();
        }

        if (param.topfdOut)
        {
          topfdOut.close();
        }

        //fsm.close();
        for (int j = 0; j < (int) param.maxMSLevel; j++)
        {
          total_specCntr[j] += specCntr[j];
          total_qspecCntr[j] += qspecCntr[j];
          total_massCntr[j] += massCntr[j];
        }
        total_featureCntr += featureCntr;
      }
      else
      {
        for (int j = 0; j < (int) param.currentMaxMSLevel; j++)
        {
          if (specCntr[j] == 0)
          {
            continue;
          }
          OPENMS_LOG_INFO << "So far, FLASHDeconv found " << massCntr[j] << " masses in " << qspecCntr[j]
                          << " MS" << (j + 1) << " spectra out of "
                          << specCntr[j] << endl;
        }
        if (featureCntr > 0)
        {
          OPENMS_LOG_INFO << "Mass tracer found " << featureCntr << " features" << endl;
        }
        for (int j = 0; j < (int) param.maxMSLevel; j++)
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
      for (int j = 0; j < (int) param.currentMaxMSLevel; j++)
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


    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    OPENMS_LOG_INFO << "Total elapsed time\n-- done [took " << total_elapsed_cpu_secs << " s (CPU), "
                    << total_elapsed_wall_secs
                    << " s (Wall)] --"
                    << endl;

    if (massCntr[0] < total_massCntr[0])
    {
      for (int j = 0; j < (int) param.currentMaxMSLevel; j++)
      {
        if (total_specCntr[j] == 0)
        {
          continue;
        }
        OPENMS_LOG_INFO << "In total, FLASHDeconv found " << total_massCntr[j] << " masses in " << total_qspecCntr[j]
                        << " MS" << (j + 1) << " spectra out of "
                        << total_specCntr[j] << endl;
      }

      if (total_featureCntr > 0)
      {
        OPENMS_LOG_INFO << "Mass tracer found " << total_featureCntr << " features" << endl;
      }
    }

    if (!isOutPathDir)
    {
      // fsm << "];";
      //  fsm.close();
      //  promexOut.close();
      for (int n = 1; n <= (int) param.maxMSLevel; ++n)
      {
        specOut[n - 1].close();

        if (specCntr[n - 1] == 0)
        {
          QString filename = QString::fromStdString(outfilePath + "_MS" + n + "_spec.tsv");
          QFile file(filename);
          file.remove();
        }
      }

      //fa.close();
      //fsm.close();
      featureOut.close();
      if (param.promexOut)
      {
        promexOut.close();
      }

      if (param.topfdOut)
      {
        topfdOut.close();
      }
    }
    return EXECUTION_OK;
  }

  static void printProgress(float progress)
  {
    //return; //
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
