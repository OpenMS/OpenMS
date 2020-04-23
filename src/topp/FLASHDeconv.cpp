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
                        "Output file prefix or output dir (if prefix, [prefix].tsv will be generated. "
                        "if dir, [dir]/[inputfile].tsv is generated per [inputfile])");

    registerDoubleList_("tol",
                        "ms1_tol ms2_tol ... (e.g., 10.0 5.0 to specify 10.0 and 5.0 ppm for MS1 and MS2, respectively",
                        {10.0, 5.0},
                        "ppm tolerance for MS1, 2, ...",
                        false);

    //registerDoubleOption_("tol", "<tolerance>", 10.0, "ppm tolerance", false, false);
    registerIntOption_("minCharge", "<min_charge>", 1, "minimum charge state", false, false);
    registerIntOption_("maxCharge", "<max_charge>", 120, "maximum charge state", false, false);
    registerDoubleOption_("minMass", "<min_mass>", 50.0, "minimum mass (Da)", false, false);
    registerDoubleOption_("maxMass", "<max_mass>", 100000.0, "maximum mass (Da)", false, false);

    registerDoubleList_("minIsotopeCosine",
                        "ms1_isotope_cos ms2_isotpe_cos ... (e.g., 0.8 0.6 to specify 0.8 and 0.6 for MS1 and MS2, respectively)",
                        {.75, .75},
                        "cosine threshold between avg. and observed isotope pattern for MS1, 2, ...",
                        false,
                        true);

    registerDoubleOption_("minChargeCosine",
                          "<charge_cosine>",
                          .5,
                          "cosine threshold between per-charge-intensity and fitted gaussian distribution (applies only to MS1)",
                          false,
                          true);

    registerIntList_("minSupportingPeaks",
                     "ms1_min_supporting_peaks ms2_min_supproting_peaks ... (e.g., 3 2 to specify 3 and 2 for MS1 and MS2, respectivly",
                     {3, 1},
                     "minimum number of supporting peaks",
                     false,
                     true);

    registerIntOption_("maxMassCount", "<max_mass_count>", -1, "maximum mass count per spec", false, true);
    //
    registerDoubleOption_("minIntensity", "<min_intensity>", 0, "intensity threshold (default 0.0)", false, true);
    registerDoubleOption_("RTwindow",
                          "<seconds>",
                          0.0,
                          "RT window (if 0, RT window contains 15 MS1 spectra)",
                          false,
                          true);
    registerDoubleOption_("minRTSpan", "<seconds>", 10.0, "Min feature RT span", false, true);
    registerIntOption_("writeDetail",
                       "<1:true 0:false>",
                       0,
                       "to write per spectrum deconvoluted masses in detail or not in [prefix]_MSn_spec.tsv files. If set, all peak information per mass is reported.",
                       false,
                       true);

    registerIntOption_("promexOutput", "", 0, "if set, promexoutput ([prefix]]_FD.ms1ft) is generated", false, true);
    registerIntOption_("topfdOutput", "", 0, "if set, topfdoutput ([prefix]_FD_ms2.msalign) is generated", false, true);
    registerIntOption_("trainOutput", "", 0, "if set, [preifx]_train_MSn.csv files are generated", false, true);

    registerIntOption_("maxMSLevel", "", 2, "maximum MS-level (inclusive) for deconvolution", false, true);
  }

  Parameter setParameter()
  {
    Parameter param;
    param.minCharge = getIntOption_("minCharge");
    param.currentChargeRange = param.chargeRange = getIntOption_("maxCharge") - param.minCharge + 1;
    param.currentMaxMass = param.maxMass = getDoubleOption_("maxMass");
    param.minMass = getDoubleOption_("minMass");
    param.tolerance = getDoubleList_("tol");

    for (auto j = 0; j < (int) param.tolerance.size(); j++)
    {
      param.tolerance[j] *= 1e-6;
      param.binWidth.push_back(.5 / param.tolerance[j]);
    }

    param.intensityThreshold = getDoubleOption_("minIntensity");
    param.minContinuousChargePeakCount = getIntList_("minSupportingPeaks");
    param.minIsotopeCosine = getDoubleList_("minIsotopeCosine");
    param.minChargeCosine = getDoubleOption_("minChargeCosine");

    param.currentMaxMassCount = param.maxMassCount = getIntOption_("maxMassCount");
    //param.chargeDistributionScoreThreshold = getDoubleOption_("minCDScore");
    param.RTwindow = getDoubleOption_("RTwindow");
    param.minRTSpan = getDoubleOption_("minRTSpan");
    //param.threads = getIntOption_("threads");
    param.writeDetail = getIntOption_("writeDetail");
    //param.jitter = getIntOption_("jitter");
    param.maxMSLevel = getIntOption_("maxMSLevel");
    param.promexOut = getIntOption_("promexOutput");
    param.topfdOut = getIntOption_("topfdOutput");
    param.trainOut = getIntOption_("trainOutput");
    return param;
  }

  static FLASHDeconvHelperStructs::PrecalcularedAveragine calculateAveragines(Parameter &param)
  {
    auto generator = new CoarseIsotopePatternGenerator();
    auto maxIso = generator->estimateFromPeptideWeight(param.maxMass);
    maxIso.trimRight(0.01 * maxIso.getMostAbundant().getIntensity());
    param.maxIsotopeCount = (int) maxIso.size() - 1;
    generator->setMaxIsotope((Size) param.maxIsotopeCount);
    return FLASHDeconvHelperStructs::PrecalcularedAveragine(50, param.maxMass, 20, generator);
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
    auto avgine = calculateAveragines(param);

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

    //int specIndex = 0, massIndex = 0;
    double total_elapsed_cpu_secs = 0, total_elapsed_wall_secs = 0;
    fstream fsf, fsp, fsfd;
    auto fs = new fstream[param.maxMSLevel];
    auto ft = new fstream[param.maxMSLevel];
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
        fs[n - 1].open(outfilePath + "_MS" + n + "_spec.tsv", fstream::out);
        DeconvolutedSpectrum::writeDeconvolutedMassesHeader(fs[n - 1], n, param.writeDetail);

        if (param.trainOut)
        {
          ft[n - 1].open(outfilePath + "_train_MS" + n + ".csv", fstream::out);//[preifx]_train_MSn.csv
          DeconvolutedSpectrum::writeAttCsvHeader(ft[n - 1]);
        }
      }

      //fsm.open(outfilePath + ".csv", fstream::out);
      //writeAttCsvHeader(fsm);

      fsf.open(outfilePath + ".tsv", fstream::out);
      MassFeatureTrace::writeHeader(fsf);

      if (param.promexOut)
      {
        fsp.open(outfilePath + "_FD.ms1ft", fstream::out);
        MassFeatureTrace::writePromexHeader(fsp);
      }

      if (param.topfdOut)
      {
        fsfd.open(outfilePath + "_FD_ms2.msalign", fstream::out);
        //writeTopFDHeader(fsfd);
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

        //specIndex = 0;
        featureCntr = 0;
        //massIndex = 0;
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
      auto msCntr = new int[1 + param.maxMSLevel];
      fill_n(msCntr, 1 + param.maxMSLevel, 0);

      for (auto &it : map)
      {
        if (it.getMSLevel() > param.maxMSLevel)
        {
          continue;
        }
        msCntr[it.getMSLevel()]++;
      }

      param.numOverlappedScans.clear();
      for (int j = 1; j <= 1; j++)
      {
        double rtDelta = rtDuration / msCntr[j];

        auto rw = param.RTwindow;
        auto count = max(param.minNumOverLappedScans, (UInt) round(rw / rtDelta));
        OPENMS_LOG_INFO << "# Overlapped MS" << j << " scans:" << count << " (in RT " << (rtDelta * count) //
                        << " sec)" << endl;

        param.numOverlappedScans.push_back(count);
      }
      for (int j = 2; j <= (int) param.maxMSLevel; j++)
      {
        param.numOverlappedScans.push_back(0);
      }
      delete[] msCntr;

      std::string outfileName(param.fileName);

      if (isOutPathDir)
      {
        std::size_t found = outfileName.find_last_of('.');
        outfileName = outfileName.substr(0, found);

        for (int n = 1; n <= (int) param.maxMSLevel; ++n)
        {
          fs[n - 1].open(outfilePath + outfileName + "_MS" + n + "_spec.tsv", fstream::out);
          DeconvolutedSpectrum::writeDeconvolutedMassesHeader(fs[n - 1], n, param.writeDetail);

          if (param.trainOut)
          {
            ft[n - 1].open(outfilePath + outfileName + "_train_MS" + n + ".csv", fstream::out);//[preifx]_train_MSn.csv
            DeconvolutedSpectrum::writeAttCsvHeader(ft[n - 1]);
          }

        }

        //fsm.open(outfilePath + outfileName + "PerSpecMasses.csv", fstream::out);
        //writeAttCsvHeader(fsm);

        fsf.open(outfilePath + outfileName + ".tsv", fstream::out);
        MassFeatureTrace::writeHeader(fsf);

        if (param.promexOut)
        {
          fsp.open(outfilePath + outfileName + "_FD.ms1ft", fstream::out);
          MassFeatureTrace::writePromexHeader(fsp);
        }

        if (param.topfdOut)
        {
          fsfd.open(outfilePath + outfileName + "_FD_ms2.msalign", fstream::out);
        }
      }

      //check max ms level from the input dataset..
      param.currentMaxMSLevel = 0;

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

      int scanNumber = 0, id = 0, specIndex = 0, massIndex = 0;
      float prevProgress = .0;
      auto massTracer = MassFeatureTrace(param, mtd_param, avgine);
      auto lastDeconvolutedSpectra = std::unordered_map<UInt, DeconvolutedSpectrum>();

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
        auto fd = FLASHDeconvAlgorithm(specIndex, massIndex, avgine, param);
        auto deconvolutedSpectrum = DeconvolutedSpectrum(*it);

        bool proceed = true;
        param.currentChargeRange = param.chargeRange;
        param.currentMaxMass = param.maxMass;

        if (msLevel > 1 && lastDeconvolutedSpectra.find(msLevel - 1) != lastDeconvolutedSpectra.end())
        {
          proceed = deconvolutedSpectrum.registerPrecursor(lastDeconvolutedSpectra[msLevel - 1]);
          if (proceed)
          {
            param.currentChargeRange = deconvolutedSpectrum.precursorPeak->charge;
            param.currentMaxMass = deconvolutedSpectrum.precursorPeakGroup->monoisotopicMass;
          }
        }

        if (proceed)
        {
          fd.Deconvolution(deconvolutedSpectrum, scanNumber);
        }

        elapsed_deconv_cpu_secs[msLevel-1] += double(clock() - deconv_begin) / CLOCKS_PER_SEC;
        elapsed_deconv_wall_secs[msLevel-1] += chrono::duration<double>(chrono::high_resolution_clock::now() - deconv_t_start)
            .count();

        if (deconvolutedSpectrum.empty())
        {
          continue;
        }

        lastDeconvolutedSpectra[msLevel] = deconvolutedSpectrum;
        massTracer.addDeconvolutedSpectrum(deconvolutedSpectrum);
        qspecCntr[msLevel - 1]++;
        massCntr[msLevel - 1] += deconvolutedSpectrum.peakGroups.size();

        deconvolutedSpectrum.writeDeconvolutedMasses(fs[msLevel - 1], param);
        if (param.topfdOut)
        {
          deconvolutedSpectrum.writeTopFD(fsfd, id++);
        }
        if(param.trainOut){
          deconvolutedSpectrum.writeAttCsv(ft[msLevel-1], msLevel);
        }

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

      massTracer.findFeatures(featureCntr, fsf, fsp);

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
          fs[n - 1].close();
          if(param.trainOut){
            ft[n-1].close();
          }
          if (specCntr[n - 1] == 0)
          {
            QString filename = QString::fromStdString(outfilePath + outfileName + "_MS" + n + "_spec.tsv");
            QFile file(filename);
            file.remove();
          }
        }
        // fsm.close();

        fsf.close();
        if (param.promexOut)
        {
          fsp.close();
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
        for (int j = 0; j < (int) param.maxMSLevel; j++)
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
      for (int j = 0; j < (int) param.maxMSLevel; j++)
      {
        sumCntr += specCntr[j];

        OPENMS_LOG_INFO << "-- deconv per MS"<< (j+1) <<" spectrum (except spec loading, feature finding) [took "
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
      for (int j = 0; j < (int) param.maxMSLevel; j++)
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
      //  fsp.close();
      for (int n = 1; n <= (int) param.maxMSLevel; ++n)
      {
        fs[n - 1].close();
        if(param.trainOut){
          ft[n-1].close();
        }
        if (specCntr[n - 1] == 0)
        {
          QString filename = QString::fromStdString(outfilePath + "_MS" + n + "_spec.tsv");
          QFile file(filename);
          file.remove();
        }
      }

      //fsm.close();
      fsf.close();
      if (param.promexOut)
      {
        fsp.close();
      }

      if (param.topfdOut)
      {
        fsfd.close();
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
