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
  typedef FLASHDeconvHelperStructs::PeakGroup PeakGroup;
  typedef FLASHDeconvHelperStructs::PrecalcularedAveragine PrecalcularedAveragine;
  typedef FLASHDeconvHelperStructs::Parameter Parameter;

  // this function will be used to register the tool parameters
  // it gets automatically called on tool execution
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<input file>", "", "Input file");
//    setValidFormats_("in", ListUtils::create<String>("mzML"),);
    registerOutputFile_("out", "<output file prefix/output dir>", "",
                        "Output file prefix or output dir (if prefix, [file prefix].tsv will be generated. "
                        "if dir, [dir]/[inputfile].tsv is generated per [inputfile])");

    registerDoubleOption_("tol", "<tolerance>", 10.0, "ppm tolerance", false, false);
    registerIntOption_("minC", "<min charge>", 2, "minimum charge state", false, false);
    registerIntOption_("maxC", "<max charge>", 100, "maximum charge state", false, false);
    registerDoubleOption_("minM", "<min mass>", 1000.0, "minimum mass (Da)", false, false);
    registerDoubleOption_("maxM", "<max mass>", 100000.0, "maximum mass (Da)", false, false);

    registerDoubleOption_("minIC", "<cosine threshold 0 - 1>", .6, "cosine threshold between avg. and observed isotope pattern", false, false);
    registerDoubleOption_("minCC", "<cosine threshold 0 - 1>", .6, "cosine threshold between per-charge-intensity and fitted gaussian distribution", false, false);
    registerDoubleOption_("minICS", "<cosine threshold 0 - 1>", .4, "cosine threshold between avg. and observed isotope pattern (spectrum level)", false, true);
    registerDoubleOption_("minCCS", "<cosine threshold 0 - 1>", .4, "cosine threshold between per-charge-intensity and fitted gaussian distribution (spectrum level)", false, true);

    registerIntOption_("minCP", "<min continuous charge peak count>", 3, "minimum number of peaks of continuous charges per mass", false, true);

    registerIntOption_("maxMC", "<max mass count>", -1, "maximum mass count per spec", false, true);
    //
    registerDoubleOption_("minIT", "<min intensity>", 0.0, "intensity threshold (default 0.0)", false, true);
    registerDoubleOption_("RTwindow", "<seconds>", 0.0, "RT window (if 0, 1% total gradient time)", false, true);
    registerDoubleOption_("minRTspan", "<seconds>", 10.0, "Min feature RT span", false, true);
    registerIntOption_("writeSpecDeconv", "<1:true 0:false>", 0, "to write per spectrum deconvoluted masses or not. If set, [prefix]PerSpecMasses.tsv is generated", false, true);

    registerIntOption_("maxMSL", "", 1, "maximum MS-level (inclusive) for deconvolution", false, true);

    registerIntOption_("jitter", "<1:true 0:false>", 0, "jitter universal pattern to generate decoy features (output file will end with *Decoy.tsv)", false, true);
  }

  Parameter setParameter()
  {
    Parameter param;
    param.minCharge = getIntOption_("minC");
    param.chargeRange = getIntOption_("maxC") - param.minCharge + 1;
    param.maxMass = getDoubleOption_("maxM");
    param.minMass = getDoubleOption_("minM");
    param.tolerance = getDoubleOption_("tol") * 1e-6;
    param.binWidth = .5 / param.tolerance;
    param.intensityThreshold = getDoubleOption_("minIT");
    param.minContinuousChargePeakCount = getIntOption_("minCP");
    param.minIsotopeCosine = getDoubleOption_("minIC");
    param.minChargeCosine = getDoubleOption_("minCC");

    param.minIsotopeCosineSpec = getDoubleOption_("minICS");
    param.minChargeCosineSpec = getDoubleOption_("minCCS");

    //  param.maxIsotopeCosine = getDoubleOption_("minIC1");
    //param.maxIsotopeCount = getIntOption_("maxIC");
    param.maxMassCount = getIntOption_("maxMC");
    //param.chargeDistributionScoreThreshold = getDoubleOption_("minCDScore");
    param.RTwindow = getDoubleOption_("RTwindow");
    param.minRTSpan = getDoubleOption_("minRTspan");
    param.threads = getIntOption_("threads");
    param.writeSpecTsv = getIntOption_("writeSpecDeconv");
    param.jitter = getIntOption_("jitter");
    param.maxMSLevel = getIntOption_("maxMSL");
    return param;
  }


  // the main_ function is called after all parameters are read
  ExitCodes main_(int, const char **) override
  {
    FLASHDeconvAlgorithm fd;

    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String infilePath = getStringOption_("in");
    String outfilePath = getStringOption_("out");

    auto param = setParameter();
    auto averagines = getPrecalculatedAveragines(param);
    int specCntr = 0, qspecCntr = 0, massCntr = 0, featureCntr = 0;
    int total_specCntr = 0, total_qspecCntr = 0, total_massCntr = 0, total_featureCntr = 0;
    double total_elapsed_cpu_secs = 0, total_elapsed_wall_secs = 0;
    fstream fs, fsf, fsm;

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
      if (param.writeSpecTsv > 0)
      {
        fs.open(outfilePath + "PerSpecMasses.tsv", fstream::out);
      }
      // if (param.RTwindow > 0)
      // {

      if (param.jitter == 0)
      {
        fsf.open(outfilePath + ".tsv", fstream::out);
      }
      else
      {
        fsf.open(outfilePath + "Decoy.tsv", fstream::out);
      }
      //  }

      writeHeader(fs, fsf, true);
      //  fsm.open(outfilePath + ".m", fstream::out);
      //  fsm << "m=[";

      // fsp.open(outfilePath + "peak.m", fstream::out);
    }

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    for (auto &infile : infileArray)
    {
      if (isOutPathDir)
      {
        specCntr = qspecCntr = massCntr = featureCntr = 0;
      }
      MSExperiment map;
      MzMLFile mzml;

      double elapsed_cpu_secs = 0, elapsed_wall_secs = 0;
      double elapsed_deconv_cpu_secs = 0, elapsed_deconv_wall_secs = 0;

      auto begin = clock();
      auto t_start = chrono::high_resolution_clock::now();

      OPENMS_LOG_INFO << "Processing : " << infile.toStdString() << endl;

      mzml.setLogType(log_type_);
      mzml.load(infile, map);

      param.fileName = QFileInfo(infile).fileName().toStdString();

      int ms1Cntr = 0;
      for (auto it = map.begin(); it != map.end(); ++it)
      {
        //cout<<it->getMSLevel()<<endl;
        if ((int) it->getMSLevel() > param.maxMSLevel)
        {
          continue;
        }
        ms1Cntr++;
      }

      double rtDuration = map[map.size() - 1].getRT() - map[0].getRT();
      double rtDelta = rtDuration / ms1Cntr;
      if (param.RTwindow <= 0)
      {
        param.RTwindow = rtDuration * .01;
      }

      param.numOverlappedScans = max(param.minNumOverLappedScans, (int) (.5 + param.RTwindow / rtDelta));
      OPENMS_LOG_INFO << "# Overlapped MS1 scans:" << param.numOverlappedScans << " (in RT " << param.RTwindow << " sec)" << endl;
      if (isOutPathDir)
      {
        std::string outfileName(param.fileName);
        std::size_t found = outfileName.find_last_of(".");
        outfileName = outfileName.substr(0, found);

        if (param.writeSpecTsv > 0)
        {
          fs.open(outfilePath + outfileName + "PerSpecMasses.tsv", fstream::out);
        }
        if (param.jitter == 0)
        {
          fsf.open(outfilePath + outfileName + ".tsv", fstream::out);
        }
        else
        {
          fsf.open(outfilePath + outfileName + "Decoy.tsv", fstream::out);
        }

        // fsm.open(outfilePath + outfileName + "Annotated.m", fstream::out); //
        writeHeader(fs, fsf, true);

        //outfileName.erase(std::remove(outfileName.begin(), outfileName.end(), '_'), outfileName.end());
        //outfileName.erase(std::remove(outfileName.begin(), outfileName.end(), '-'), outfileName.end());
        //   fsm.open(outfilePath + "m" + outfileName + ".m", fstream::out);
        //   fsm << "m=[";
        //   fsp.open(outfilePath + "m" + outfileName + "peak.m", fstream::out);
      }

      OPENMS_LOG_INFO << "Running FLASHDeconv ... " << endl;
      auto deconv_begin = clock();
      auto deconv_t_start = chrono::high_resolution_clock::now();
      //continue;
      auto peakGroups = FLASHDeconvAlgorithm::Deconvolution(map, param, averagines, specCntr, qspecCntr, massCntr);

      auto deconv_t_end = chrono::high_resolution_clock::now();
      auto deconv_end = clock();

      //writeAnnotatedSpectra(peakGroups,map,fsm);//

      if (!peakGroups.empty() && specCntr > 0 && map.size() > 1)
      {
        findFeatures(peakGroups, map, featureCntr, fsf, averagines, param);
      }

      cout<< "after running" << endl;

      if (param.writeSpecTsv)
      {
        OPENMS_LOG_INFO << endl << "writing per spec deconvolution results ...";
        OPENMS_LOG_INFO.flush();

        for (auto &pg : peakGroups)
        {
          writePeakGroup(pg, param, fs);
        }

        OPENMS_LOG_INFO << "done" << endl;

      }

      // cout<<4.5<<endl;
      if (isOutPathDir)
      {
        OPENMS_LOG_INFO << "In this run, FLASHDeconv found " << massCntr << " masses in " << qspecCntr
             << " MS1 spectra out of "
             << specCntr << endl;
        if (featureCntr > 0)
        {
          OPENMS_LOG_INFO << "Mass tracer found " << featureCntr << " features" << endl;
        }

        //   fsm << "];";
        //   fsm.close();
        //   fsp.close();
        if (param.writeSpecTsv > 0)
        {
          fs.close();
        }

        fsf.close();

        //fsm.close();
        total_specCntr += specCntr;
        total_qspecCntr += qspecCntr;
        total_massCntr += massCntr;
        total_featureCntr += featureCntr;
      }
      else
      {
        OPENMS_LOG_INFO << "So far, FLASHDeconv found " << massCntr << " masses in " << qspecCntr
             << " MS1 spectra out of "
             << specCntr << endl;
        if (featureCntr > 0)
        {
          OPENMS_LOG_INFO << "Mass tracer found " << featureCntr << " features" << endl;
        }

        total_specCntr = specCntr;
        total_qspecCntr = qspecCntr;
        total_massCntr = massCntr;
        total_featureCntr = featureCntr;

      }

      //cout<<5<<endl;
      auto t_end = chrono::high_resolution_clock::now();
      auto end = clock();

      elapsed_deconv_cpu_secs = double(deconv_end - deconv_begin) / CLOCKS_PER_SEC;
      elapsed_deconv_wall_secs = chrono::duration<double>(deconv_t_end - deconv_t_start).count();

      elapsed_cpu_secs = double(end - begin) / CLOCKS_PER_SEC;
      elapsed_wall_secs = chrono::duration<double>(t_end - t_start).count();

      OPENMS_LOG_INFO << "-- done [took " << elapsed_cpu_secs << " s (CPU), " << elapsed_wall_secs
           << " s (Wall)] --"
           << endl;
      OPENMS_LOG_INFO << "-- deconv per spectrum (except spec loading, feature finding) [took "
           << 1000.0 * elapsed_deconv_cpu_secs / specCntr
           << " ms (CPU), " << 1000.0 * elapsed_deconv_wall_secs / specCntr << " ms (Wall)] --" << endl;

      total_elapsed_cpu_secs += elapsed_cpu_secs;
      total_elapsed_wall_secs += elapsed_wall_secs;
    }


    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    OPENMS_LOG_INFO << "Total elapsed time\n-- done [took " << total_elapsed_cpu_secs << " s (CPU), " << total_elapsed_wall_secs
         << " s (Wall)] --"
         << endl;

    if (massCntr < total_massCntr)
    {
      OPENMS_LOG_INFO << "In total, FLASHDeconv found " << total_massCntr << " masses in " << total_qspecCntr
           << " MS1 spectra out of "
           << total_specCntr << endl;
      if (featureCntr > 0)
      {
        OPENMS_LOG_INFO << "Mass tracer found " << total_featureCntr << " features" << endl;
      }
    }

    if (!isOutPathDir)
    {
      // fsm << "];";
      //  fsm.close();
      //  fsp.close();
      if (param.writeSpecTsv > 0)
      {
        fs.close();
      }
      fsf.close();
    }

    return EXECUTION_OK;
  }

  static void writeAnnotatedSpectra(vector<PeakGroup> &pgs,
                                    MSExperiment &map,
                                    fstream &fs)//, fstream &fsm, fstream &fsp)
  {

    boost::unordered_map<double, vector<PeakGroup>> pgmap;

    for (auto &pg : pgs)
    {
      pgmap[pg.spec->getRT()].push_back(pg);
    }

    int index = 1;
    for (auto it = map.begin(); it != map.end(); ++it)
    {
      double rt = it->getRT();
      /*if (rt < 3250)
      {
        continue;
      }
      if (rt > 3350)
      {
        break;
      }*/

      if (pgmap.find(rt) == pgmap.end())
      {
        continue;
      }
      auto t = pgmap[rt];
      vector<double> masses;
      int mi = 1;
      for (auto &p : t)
      {
        fs << "apeaks" << index << "{" << mi++ << "}=[";
        for (auto &lp:p.peaks)
        {
          auto &op = lp.orgPeak;
          fs << op->getMZ() << "," << op->getIntensity() << ";";
        }

        fs << "];\n";
      }
      fs << "aspec{" << index << "}=[";
      for (auto &p : t)
      {
        fs << p.monoisotopicMass << "," << p.intensity << ";";
      }
      fs << "];\n";

      fs << "spec{" << index << "}=[";
      for (auto &p : *it)
      {
        fs << p.getMZ() << "," << p.getIntensity() << ";";
      }
      fs << "];\n";


      index++;
      //if (my_hash_map.find(non-existent key) == my_hash_map.end())

    }


  }

  static void writePeakGroup(PeakGroup &pg, Parameter &param, fstream &fs)//, fstream &fsm, fstream &fsp)
  {
    //return;//
    if (pg.peaks.empty())
    {
      return;
    }
    double &m = pg.monoisotopicMass;
    double &am = pg.avgMass;
    double &intensity = pg.intensity;
    //int nm = getNominalMass(m);
    sort(pg.peaks.begin(), pg.peaks.end());
    int minCharge = param.chargeRange + param.minCharge;
    int maxCharge = -1;
    for (auto &p : pg.peaks)
    {
      minCharge = minCharge < p.charge ? minCharge : p.charge;
      maxCharge = maxCharge > p.charge ? maxCharge : p.charge;
    }
    //cout<<1<<endl;
    fs << pg.massIndex << "\t" << pg.specIndex << "\t" << param.fileName << "\t" << pg.spec->getNativeID() << "\t"
       << pg.massCntr << "\t"
       << fixed << setprecision(3) << m << "\t" << am << "\t"
       << (maxCharge - minCharge + 1) << "\t" << minCharge << "\t" << maxCharge << "\t"
       << fixed << setprecision(1) << intensity << "\t" << pg.spec->getRT()
       << "\t" << pg.peaks.size() << "\t";

    fs << fixed << setprecision(2);
    for (auto &p : pg.peaks)
    {
      fs << p.orgPeak->getMZ() << ";";
    }
    fs << "\t";
    for (auto &p : pg.peaks)
    {
      fs << p.charge << ";";
    }
    fs << "\t";
    for (auto &p : pg.peaks)
    {
      fs << p.getMass() << ";";
    }
    fs << "\t";
    for (auto &p : pg.peaks)
    {
      fs << p.isotopeIndex << ";";
    }
    fs << "\t";
    fs << fixed << setprecision(1);
    for (auto &p : pg.peaks)
    {
      fs << p.orgPeak->getIntensity() << ";";
    }
    fs << fixed << setprecision(3);
    fs << "\t" << pg.isotopeCosineScore
       << "\t" << pg.chargeCosineScore
       << "\n";

    /*
    //cout<<1<<endl;

            fsp << "pg" << (int) (pg.monoisotopicMass * 10) << "rt" << (int) (pg.spec->getRT())
                << "=[";

            for (auto &p : pg.peaks) {
                fsp << p.charge << "," << p.isotopeIndex << "," << p.orgPeak->getIntensity() << ";";
            }

            fsp << "];\n";
            //cout<<3<<endl;
 */
    //fsm << m << "," << nm << "," << intensity << "," << pg.spec->getRT() << "\n";
    //cout<<4<<endl;


  }

  void findFeatures(vector<PeakGroup> &peakGroups,
                    MSExperiment &map,
                    int &featureCntr,
                    fstream &fsf,
                    PrecalcularedAveragine &averagines,
                    Parameter &param)
  {
    boost::unordered_map<float, PeakGroup> *peakGroupMap;
    boost::unordered_map<float, int> rtSpecMap;

    for (auto it = map.begin(); it != map.end(); ++it)
    {
      it->clear(false);
    }

    int maxSpecIndex = 0;
    for (auto &pg : peakGroups)
    {
      auto &spec = pg.spec;

      Peak1D tp(pg.monoisotopicMass, (float) pg.intensity);

      rtSpecMap[spec->getRT()] = pg.specIndex;
      maxSpecIndex = max(maxSpecIndex, pg.specIndex);

      spec->push_back(tp);
    }
    peakGroupMap = new boost::unordered_map<float, PeakGroup>[maxSpecIndex + 1];

    for (auto &pg : peakGroups)
    {
      //      auto &spec = pg.spec;

      auto &pgMap = peakGroupMap[pg.specIndex];

      pgMap[pg.monoisotopicMass] = pg;
    }

    for (auto it = map.begin(); it != map.end(); ++it)
    {
      it->sortByPosition();
    }

    Param common_param = getParam_().copy("algorithm:common:", true);
    writeDebug_("Common parameters passed to sub-algorithms (mtd and ffm)", common_param, 3);

    Param mtd_param = getParam_().copy("algorithm:mtd:", true);
    writeDebug_("Parameters passed to MassTraceDetection", mtd_param, 3);

    MassTraceDetection mtdet;
    mtd_param.insert("", common_param);
    mtd_param.remove("chrom_fwhm");

    //mtd_param.setValue("mass_error_da", .3,// * (param.chargeRange+ param.minCharge),
    //                   "Allowed mass deviation (in da).");
    mtd_param.setValue("mass_error_ppm", param.tolerance * 1e6 * 2, "");
    mtd_param.setValue("trace_termination_criterion", "outlier", "");

    mtd_param.setValue("reestimate_mt_sd", "true", "");
    mtd_param.setValue("quant_method", "area", "");
    mtd_param.setValue("noise_threshold_int", .0, "");

    //double rtDuration = (map[map.size() - 1].getRT() - map[0].getRT()) / ms1Cntr;
    mtd_param.setValue("min_sample_rate", 0.01, "");
    mtd_param.setValue("trace_termination_outliers", param.numOverlappedScans, "");
    mtd_param.setValue("min_trace_length", param.minRTSpan, "");
    //mtd_param.setValue("max_trace_length", 1000.0, "");
    mtdet.setParameters(mtd_param);

    vector<MassTrace> m_traces;
    mtdet.run(map, m_traces);  // m_traces : output of this function

    // cout<<1<<endl;
    double *perChargeIntensity = new double[param.chargeRange + param.minCharge + 1];
    double *perChargeMaxIntensity = new double[param.chargeRange + param.minCharge + 1];
    double *perChargeMz = new double[param.chargeRange + param.minCharge + 1];
    double *perIsotopeIntensity = new double[param.maxIsotopeCount];

    for (auto &mt : m_traces)
    {
      //if (mt.getSize() < 3)
      //{
      //continue;
      // }
      int minCharge = param.chargeRange + param.minCharge + 1;
      int maxCharge = 0;
      boost::dynamic_bitset<> charges(param.chargeRange + param.minCharge + 1);
      cout << "where? 1" << endl;
      fill_n(perChargeIntensity, param.chargeRange + param.minCharge + 1, 0);
      fill_n(perChargeMaxIntensity, param.chargeRange + param.minCharge + 1, 0);
      fill_n(perChargeMz, param.chargeRange + param.minCharge + 1, 0);
      fill_n(perIsotopeIntensity, param.maxIsotopeCount, 0);
      //mt.getIntensity()
      //double sum_intensity = .0;
      //double mass = 0;//mt.getCentroidMZ();
      //double avgMass = 0;
      double massDiff = 0;
      double max_intensity = -1;

      for (auto &p2 : mt)
      {
        cout << "where? 2" << endl;
        int specIndex = rtSpecMap[(float) p2.getRT()];
        auto &pgMap = peakGroupMap[specIndex];
        auto &pg = pgMap[(float) p2.getMZ()];
        minCharge = min(minCharge, pg.minCharge);
        maxCharge = max(maxCharge, pg.maxCharge);
        //sum_intensity += pg.intensity;
        if (pg.intensity > max_intensity)
        {
          max_intensity = pg.intensity;
          massDiff = pg.avgMass - pg.monoisotopicMass;
        }
        cout << "where? 22" << endl;
        for (auto &p : pg.peaks)
        {
          if (p.isotopeIndex < 0 || p.isotopeIndex >= param.maxIsotopeCount || p.charge < 0 ||
              p.charge >= param.chargeRange + param.minCharge + 1)
          {
            continue;
          }
          cout << "where? 222" << endl;
          charges[p.charge] = true;
          cout << "p.charge:" << p.charge << "," << p.orgPeak->getIntensity() << endl;
          perChargeIntensity[p.charge] += p.orgPeak->getIntensity();
          cout << "where? 222-2" << endl;
          perIsotopeIntensity[p.isotopeIndex] += p.orgPeak->getIntensity();
          cout << "where? 222-3" << endl;
          if (perChargeMaxIntensity[p.charge] > p.orgPeak->getIntensity())
          {
            continue;
          }
          cout << "where? 2222" << endl;
          perChargeMaxIntensity[p.charge] = p.orgPeak->getIntensity();
          perChargeMz[p.charge] = p.orgPeak->getMZ();
          cout << "where? 3" << endl;
        }
        cout << "where? 4" << endl;

        /*if (max_intensity > pg.intensity)
        {
          continue;
        }
        max_intensity = pg.intensity;
        mass = pg.monoisotopicMass;*/
      }
      cout << "where? 5" << endl;

      // cout<<2<<endl;
      if (massDiff <= 0)
      {
        continue;
      }

      double chargeScore = FLASHDeconvAlgorithm::getChargeFitScore(perChargeIntensity, param.minCharge + param.chargeRange + 1);
      if (chargeScore < param.minChargeCosine) //
      {
        continue;
      }

      int offset = 0;
      double mass = mt.getCentroidMZ();
      double isoScore = FLASHDeconvAlgorithm::getIsotopeCosineAndDetermineIsotopeIndex(mass,
                                                                 perIsotopeIntensity,
                                                                 param.maxIsotopeCount,
                                                                 averagines, offset);
      cout << "where? 6" << endl;
      if (isoScore < param.minIsotopeCosine)
      {
        continue;
      }

      if (offset != 0)
      {
        mass += offset * Constants::C13C12_MASSDIFF_U;
        //avgMass += offset * Constants::C13C12_MASSDIFF_U;
        //p.isotopeIndex -= offset;
      }
      cout << "where? 7" << endl;
      //auto mass = mt.getCentroidMZ();
      fsf << ++featureCntr << "\t" << param.fileName << "\t" << to_string(mass) << "\t"
          << to_string(mass + massDiff) << "\t"
          << mt.getSize() << "\t"
          //fsf << ++featureCntr << "\t" << param.fileName << "\t" << mass << "\t"
          //<< getNominalMass(mass) << "\t"
          << mt.begin()->getRT() << "\t"
          << mt.rbegin()->getRT() << "\t"
          << mt.getTraceLength() << "\t"
          << mt[mt.findMaxByIntPeak()].getRT() << "\t"
          << mt.getMaxIntensity(false) << "\t"
          // << mt.computePeakArea() << "\t"
          << minCharge << "\t"
          << maxCharge << "\t"
          << charges.count() << "\t"
          << isoScore << "\t"
          << chargeScore << "\n";

      // cout<<3<<endl;

      /*for (int charge = param.minCharge; charge < param.chargeRange + param.minCharge; charge++)
      {
        if (perChargeIntensity[charge] <= 0)
        {
          continue;
        }
        fsf << charge << "," << perChargeMz[charge] << "," << perChargeIntensity[charge] << ";";
      }
      fsf << param.chargeRange + param.minCharge << "," <<
          perChargeMz[param.chargeRange + param.minCharge] << "," <<
          perChargeIntensity[param.chargeRange + param.minCharge] << "\n";
  */
      cout << "where? 8" << endl;
    }
    // cout<<4<<endl;
    delete[] perIsotopeIntensity;
    delete[] perChargeMz;
    delete[] perChargeMaxIntensity;
    delete[] perChargeIntensity;

    delete[] peakGroupMap;
    // cout<<4.1<<endl;
    cout << "where? 9" << endl;
  }

  static PrecalcularedAveragine getPrecalculatedAveragines(Parameter &param)
  {
    auto generator = new CoarseIsotopePatternGenerator();
    auto maxIso = generator->estimateFromPeptideWeight(param.maxMass);
    maxIso.trimRight(0.01 * maxIso.getMostAbundant().getIntensity());
    param.maxIsotopeCount = (int) maxIso.size() - 1;
    generator->setMaxIsotope((Size) param.maxIsotopeCount);
    return PrecalcularedAveragine(100, param.maxMass, 50, generator);
  }

  static void writeHeader(fstream &fs, fstream &fsf, bool featureOut = false)
  {
    fs
        << "MassIndex\tSpecIndex\tFileName\tSpecID\tMassCountInSpec\tMonoisotopicMass\tAvgMass\t"
           "PeakChargeRange\tPeakMinCharge\tPeakMaxCharge\t"
           "AggregatedIntensity\tRetentionTime\tPeakCount\tPeakMZs\tPeakCharges\tPeakMasses\tPeakIsotopeIndices\t"
           "PeakIntensities\tIsotopeCosineScore\tChargeIntensityCosineScore\n";
    if (!featureOut)
    {
      return;
    }
    fsf << "ID\tFileName\tMonoisotopicMass\tAverageMass\tMassCount\tStartRetentionTime"
           "\tEndRetentionTime\tRetentionTimeDuration\tApexRetentionTime"
           "\tMaxIntensity\tMinCharge\tMaxCharge\tChargeCount\tIsotopeCosineScore\tChargeIntensityCosineScore\n";

    return;
  }




};



// the actual main function needed to create an executable
int main(int argc, const char ** argv)
{
  TOPPFLASHDeconv tool;
  return tool.main(argc, argv);
}
