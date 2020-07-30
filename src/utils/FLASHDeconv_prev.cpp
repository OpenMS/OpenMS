//
// Created by JihyungKim on 07.12.18.
//

#include <OpenMS/config.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include "boost/dynamic_bitset.hpp"
#include <iostream>
#include <iomanip>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <QDirIterator>
#include <QFileInfo>
#include <chrono>
#include <OpenMS/FILTERING/DATAREDUCTION/MassTraceDetection.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
//#include <OpenMS/MATH/STATISTICS/CumulativeBinomial.h>
#include <Eigen/Dense>

using namespace OpenMS;
using namespace std;

class FLASHDeconv :
    public TOPPBase
{
public:
  FLASHDeconv() :
      TOPPBase("FLASHDeconv",
               "Ultra-fast high-quality deconvolution enables online processing of top-down MS data",
               false)
  {
  }

  //static const bool jitter = false;

  struct Parameter
  {
    int minCharge;
    double minMass;
    double maxMass;
    double tolerance;
    String fileName;// up to here: ordinary user accessible parameters

    double intensityThreshold;// advanced parameters
    double minIsotopeCosine;
    double minChargeCosine;

    double minIsotopeCosineSpec;
    double minChargeCosineSpec;

    int minContinuousChargePeakCount;
    int maxIsotopeCount;
    int maxMassCount;
    int maxMSLevel = 1;//maxMSL;
    //double charg = 1eDistributionScoreThreshold;
    double RTwindow;
    double minRTSpan;
    vector<int> hCharges{2, 3, 5,}; // automated or fixed parameters
    int chargeRange;
    double binWidth;
    int minNumOverLappedScans = 1;
    int numOverlappedScans = minNumOverLappedScans;
    int threads = 1;
    int writeSpecTsv = 0;
    int jitter = 0;
  };

  struct PrecalcularedAveragine
  {
    vector<IsotopeDistribution> isotopes;
    vector<double> norms;
    vector<Size> mostAbundantIndices;

    double massInterval;
    double minMass;

    PrecalcularedAveragine(double m, double M, double delta, CoarseIsotopePatternGenerator *generator)
        :
        massInterval(delta), minMass(m)
    {
      int i = 0;
      while (true)
      {
        double a = i * massInterval;
        i++;
        if (a < m)
        {
          continue;
        }
        if (a > M)
        {
          break;
        }
        auto iso = generator->estimateFromPeptideWeight(a);
        iso.trimRight(0.01 * iso.getMostAbundant().getIntensity());
        isotopes.push_back(iso);
        double norm = .0;
        Size mostAbundantIndex = 0;
        double mostAbundantInt = 0;

        for (Size k = 0; k < iso.size(); k++)
        {
          norm += iso[k].getIntensity() * iso[k].getIntensity();
          if (mostAbundantInt >= iso[k].getIntensity())
          {
            continue;
          }
          mostAbundantInt = iso[k].getIntensity();
          mostAbundantIndex = k;
        }
        mostAbundantIndices.push_back(mostAbundantIndex);
        norms.push_back(norm);


        //auto mostAbundant = iso.getMostAbundant();


      }
    }

    IsotopeDistribution get(double mass)
    {
      Size i = (Size) (.5 + (mass - minMass) / massInterval);
      i = i >= isotopes.size() ? isotopes.size() - 1 : i;
      return isotopes[i];
    }

    double getNorm(double mass)
    {
      Size i = (Size) (.5 + (mass - minMass) / massInterval);
      i = i >= isotopes.size() ? isotopes.size() - 1 : i;
      return norms[i];
    }

    Size getMostAbundantIndex(double mass)
    {
      Size i = (Size) (.5 + (mass - minMass) / massInterval);
      i = i >= isotopes.size() ? isotopes.size() - 1 : i;
      return mostAbundantIndices[i];
    }

    void print(double mass)
    {
      auto iso = get(mass);
      cout << "avgiso=[";
      for (Size k = 0; k < iso.size(); k++)
      {
        cout << iso[k].getIntensity() << " ";
      }
      cout << "];sqrdsumint=" << getNorm(mass) << ";" << endl;
    }

  };

  struct LogMzPeak
  {
    Peak1D *orgPeak;
    double logMz = 0;
    double mass = .0;
    int charge = 0;
    int isotopeIndex = -1;
    //double score;

    LogMzPeak() :
        orgPeak(nullptr), logMz(-1000), charge(0), isotopeIndex(0)
    {
    }

    explicit LogMzPeak(Peak1D &peak) :
        orgPeak(&peak), logMz(getLogMz(peak.getMZ())), charge(0), isotopeIndex(0)
    {
    }

    LogMzPeak(Peak1D &peak, int c, int i) :
        orgPeak(&peak), logMz(getLogMz(peak.getMZ())), charge(c),
        isotopeIndex(i)
    {
    }

    ~LogMzPeak()
    {

    }

    double getMass()
    {
      if (mass <= 0)
      {
        mass = exp(logMz) * charge;
      }
      return mass;
    }

    bool operator<(const LogMzPeak &a) const
    {
      return this->logMz < a.logMz;
    }

    bool operator>(const LogMzPeak &a) const
    {
      return this->logMz >= a.logMz;
    }

  };

  struct PeakGroup
  {
    vector<LogMzPeak> peaks;
    double monoisotopicMass = .0;
    double avgMass = .0;
    double intensity = .0;
    Size massBinIndex = 0;
    // double chargeDistributionScore = .0;

    double isotopeCosineScore = .0;
    double chargeCosineScore = .0;
    int massIndex, specIndex, massCntr;
    int maxCharge, minCharge;

    MSSpectrum *spec;

    ~PeakGroup()
    {
      vector<LogMzPeak>().swap(peaks);
    }

    void push_back(LogMzPeak &p)
    {
      peaks.push_back(p);
    }

    void reserve(Size n)
    {
      peaks.reserve(n);
    }

    bool operator<(const PeakGroup &a) const
    {
      return this->spec->getRT() < a.spec->getRT();
    }

    bool operator>(const PeakGroup &a) const
    {
      return this->spec->getRT() >= a.spec->getRT();
    }


    bool operator==(const PeakGroup &a) const
    {
      return this->spec->getRT() == a.spec->getRT() && this->monoisotopicMass == a.monoisotopicMass
             && this->intensity == a.intensity;
    }

    void updateMassesAndIntensity(PrecalcularedAveragine &averagines, int offset = 0, int maxIsoIndex = 0)
    {
      double maxIntensityForMonoIsotopeMass = -1;

      if (offset != 0)
      {
        vector<LogMzPeak> tmpPeaks;
        tmpPeaks.swap(peaks);
        peaks.reserve(tmpPeaks.size());

        for (auto &p : tmpPeaks)
        {
          p.isotopeIndex -= offset;
          if (p.isotopeIndex < 0 || p.isotopeIndex >= maxIsoIndex)
          {
            continue;
          }
          peaks.push_back(p);
        }
      }

      intensity = .0;
      for (auto &p : peaks)
      {
        double pi = p.orgPeak->getIntensity();
        intensity += pi;
        if (maxIntensityForMonoIsotopeMass > pi)
        {
          continue;
        }
        maxIntensityForMonoIsotopeMass = pi;
        monoisotopicMass = p.getMass() - p.isotopeIndex * Constants::C13C12_MASSDIFF_U;
        //if (updateAvgMass)
        //{
        int mostAbundantIndex = averagines.getMostAbundantIndex(monoisotopicMass);
        avgMass = p.getMass() + (mostAbundantIndex - p.isotopeIndex) * Constants::C13C12_MASSDIFF_U;
        //}
      }
    }
  };


protected:

  static double getLogMz(double mz)
  {
    return log(mz - Constants::PROTON_MASS_U);
  }

  static int getNominalMass(double &m)
  {
    return (int) (m * 0.999497 + .5);
  }

  static double getBinValue(Size bin, double minV, double binWidth)
  {
    return minV + bin / binWidth;
  }

  static Size getBinNumber(double v, double minV, double binWidth)
  {
    if (v < minV)
    {
      return 0;
    }
    return (Size) (((v - minV) * binWidth) + .5);

  }


  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<input file>", "", "Input file");
    //setValidFormats_("in", ListUtils::create<String>("mzML"),);
    registerOutputFile_("out",
                        "<output file prefix/output dir>",
                        "",
                        "Output file prefix or output dir (if prefix, [file prefix].tsv will be generated. "
                        "if dir, [dir]/[inputfile].tsv is generated per [inputfile])");
    registerDoubleOption_("tol", "<tolerance>", 10.0, "ppm tolerance", false, false);
    registerIntOption_("minC", "<min charge>", 2, "minimum charge state", false, false);
    registerIntOption_("maxC", "<max charge>", 100, "maximum charge state", false, false);
    registerDoubleOption_("minM", "<min mass>", 1000.0, "minimum mass (Da)", false, false);
    registerDoubleOption_("maxM", "<max mass>", 100000.0, "maximum mass (Da)", false, false);

    registerDoubleOption_("minIC",
                          "<cosine threshold 0 - 1>",
                          .6,
                          "cosine threshold between avg. and observed isotope pattern",
                          false,
                          false);
    registerDoubleOption_("minCC",
                          "<cosine threshold 0 - 1>",
                          .6,
                          "cosine threshold between per-charge-intensity and fitted gaussian distribution",
                          false,
                          false);


    registerDoubleOption_("minICS",
                          "<cosine threshold 0 - 1>",
                          .4,
                          "cosine threshold between avg. and observed isotope pattern (spectrum level)",
                          false,
                          true);
    registerDoubleOption_("minCCS",
                          "<cosine threshold 0 - 1>",
                          .4,
                          "cosine threshold between per-charge-intensity and fitted gaussian distribution (spectrum level)",
                          false,
                          true);

    registerIntOption_("minCP", "<min continuous charge peak count>", 3,
                       "minimum number of peaks of continuous charges per mass", false, true);


    registerIntOption_("maxMC", "<max mass count>", -1, "maximum mass count per spec", false, true);
    //
    registerDoubleOption_("minIT", "<min intensity>", 0.0, "intensity threshold (default 0.0)", false, true);
    registerDoubleOption_("RTwindow", "<seconds>", 0.0,
                          "RT window (if 0, 1% total gradient time)", false, true);
    registerDoubleOption_("minRTspan", "<seconds>", 10.0,
                          "Min feature RT span", false, true);
    registerIntOption_("writeSpecDeconv",
                       "<1:true 0:false>",
                       0,
                       "to write per spectrum deconvoluted masses or not. If set, [prefix]PerSpecMasses.tsv is generated",
                       false,
                       true);

    registerIntOption_("maxMSL",
                       "",
                       1,
                       "maximum MS-level (inclusive) for deconvolution",
                       false,
                       true);

    registerIntOption_("jitter",
                       "<1:true 0:false>",
                       0,
                       "jitter universal pattern to generate decoy features (output file will end with *Decoy.tsv)",
                       false,
                       true);


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

  ExitCodes main_(int, const char **) override
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String infilePath = getStringOption_("in");
    String outfilePath = getStringOption_("out");

    //-------------------------------------------------------------
    // input file path --> put in array
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

    cout << "Initializing ... " << endl;

    auto param = setParameter();
    auto averagines = getPrecalculatedAveragines(param);
    int specCntr = 0, qspecCntr = 0, massCntr = 0, featureCntr = 0;
    int total_specCntr = 0, total_qspecCntr = 0, total_massCntr = 0, total_featureCntr = 0;
    double total_elapsed_cpu_secs = 0, total_elapsed_wall_secs = 0;
    fstream fs, fsf, fsm;

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

      cout << "Processing : " << infile.toStdString() << endl;

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
      cout << "# Overlapped MS1 scans:" << param.numOverlappedScans << " (in RT " << param.RTwindow << " sec)" << endl;
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

      cout << "Running FLASHDeconv ... " << endl;
      auto deconv_begin = clock();
      auto deconv_t_start = chrono::high_resolution_clock::now();
      //continue;
      auto peakGroups = Deconvolution(map, param, averagines, specCntr, qspecCntr, massCntr);

      auto deconv_t_end = chrono::high_resolution_clock::now();
      auto deconv_end = clock();

      //writeAnnotatedSpectra(peakGroups,map,fsm);//

      if (!peakGroups.empty() && specCntr > 0 && map.size() > 1)
      {
        findFeatures(peakGroups, map, featureCntr, fsf, averagines, param);
      }

      if (param.writeSpecTsv)
      {
        cout << endl << "writing per spec deconvolution results ...";
        cout.flush();

        for (auto &pg : peakGroups)
        {
          writePeakGroup(pg, param, fs);
        }

        cout << "done" << endl;

      }

      // cout<<4.5<<endl;
      if (isOutPathDir)
      {
        cout << "In this run, FLASHDeconv found " << massCntr << " masses in " << qspecCntr
             << " MS1 spectra out of "
             << specCntr << endl;
        if (featureCntr > 0)
        {
          cout << "Mass tracer found " << featureCntr << " features" << endl;
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
        cout << "So far, FLASHDeconv found " << massCntr << " masses in " << qspecCntr
             << " MS1 spectra out of "
             << specCntr << endl;
        if (featureCntr > 0)
        {
          cout << "Mass tracer found " << featureCntr << " features" << endl;
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

      cout << "-- done [took " << elapsed_cpu_secs << " s (CPU), " << elapsed_wall_secs
           << " s (Wall)] --"
           << endl;
      cout << "-- deconv per spectrum (except spec loading, feature finding) [took "
           << 1000.0 * elapsed_deconv_cpu_secs / specCntr
           << " ms (CPU), " << 1000.0 * elapsed_deconv_wall_secs / specCntr << " ms (Wall)] --" << endl;

      total_elapsed_cpu_secs += elapsed_cpu_secs;
      total_elapsed_wall_secs += elapsed_wall_secs;
    }

    cout << "Total elapsed time\n-- done [took " << total_elapsed_cpu_secs << " s (CPU), " << total_elapsed_wall_secs
         << " s (Wall)] --"
         << endl;

    if (massCntr < total_massCntr)
    {
      cout << "In total, FLASHDeconv found " << total_massCntr << " masses in " << total_qspecCntr
           << " MS1 spectra out of "
           << total_specCntr << endl;
      if (featureCntr > 0)
      {
        cout << "Mass tracer found " << total_featureCntr << " features" << endl;
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

  void
  findFeatures(vector<PeakGroup> &peakGroups,
               MSExperiment &mp,
               int &featureCntr,
               fstream &fsf,
               PrecalcularedAveragine &averagines,
               Parameter &param)
  {

    MSExperiment map;
    boost::unordered_map<float, PeakGroup> *peakGroupMap;
    // boost::unordered_map<float, MSSpectrum> rtOrignalSpecMap;
    boost::unordered_map<float, int> rtSpecMap;
    std::map<int, MSSpectrum> indexSpecMap;

    int maxSpecIndex = 0;

    for (auto &pg : peakGroups)
    {
      auto &spec = pg.spec;

      if (indexSpecMap.find(pg.specIndex) == indexSpecMap.end()){
        indexSpecMap[pg.specIndex] = MSSpectrum();
      }
      auto &deconvSpec = indexSpecMap[pg.specIndex];
      rtSpecMap[spec->getRT()] = pg.specIndex;
      maxSpecIndex = maxSpecIndex > pg.specIndex ? maxSpecIndex : pg.specIndex ;

      deconvSpec.setRT(spec->getRT());
      Peak1D tp(pg.monoisotopicMass, (float) pg.intensity);
      deconvSpec.push_back(tp);
    }

    //int tmp = 0;
    for(auto iter = indexSpecMap.begin(); iter != indexSpecMap.end(); ++iter)
    {
      //tmp+=(iter->second).size();
      map.addSpectrum(iter->second);
    }

    //std::cout<<map.size()<< " " <<tmp<< std::endl;
    peakGroupMap = new boost::unordered_map<float, PeakGroup>[maxSpecIndex + 1];

    for (auto &pg : peakGroups)
    {
      //      auto &spec = pg.spec;
      auto &pgMap = peakGroupMap[pg.specIndex];

      pgMap[pg.monoisotopicMass] = pg;
      //std::cout<<pg.monoisotopicMass<< " " << pg.specIndex << std::endl;

    }

    for (auto it = map.begin(); it != map.end(); ++it)
    {
      it->sortByPosition();
      // cout<<it->size()<<endl;
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

        for (auto &p : pg.peaks)
        {
          if (p.isotopeIndex < 0 || p.isotopeIndex >= param.maxIsotopeCount || p.charge < 0 ||
              p.charge >= param.chargeRange + param.minCharge + 1)
          {
            continue;
          }
          charges[p.charge] = true;
          perChargeIntensity[p.charge] += p.orgPeak->getIntensity();
          perIsotopeIntensity[p.isotopeIndex] += p.orgPeak->getIntensity();
          if (perChargeMaxIntensity[p.charge] > p.orgPeak->getIntensity())
          {
            continue;
          }
          perChargeMaxIntensity[p.charge] = p.orgPeak->getIntensity();
          perChargeMz[p.charge] = p.orgPeak->getMZ();
        }
        /*if (max_intensity > pg.intensity)
        {
          continue;
        }
        max_intensity = pg.intensity;
        mass = pg.monoisotopicMass;*/
      }
      // cout<<2<<endl;
      if (massDiff <= 0)
      {
        continue;
      }

      double chargeScore = getChargeFitScore(perChargeIntensity, param.minCharge + param.chargeRange + 1);
      if (chargeScore < param.minChargeCosine) //
      {
        continue;
      }

      int offset = 0;
      double mass = mt.getCentroidMZ();
      double isoScore = getIsotopeCosineAndDetermineIsotopeIndex(mass,
                                                                 perIsotopeIntensity,
                                                                 param.maxIsotopeCount,
                                                                 averagines, offset);

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
    }
    // cout<<4<<endl;
    delete[] perIsotopeIntensity;
    delete[] perChargeMz;
    delete[] perChargeMaxIntensity;
    delete[] perChargeIntensity;

    delete[] peakGroupMap;
    // cout<<4.1<<endl;
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

  static PrecalcularedAveragine getPrecalculatedAveragines(Parameter &param)
  {
    auto generator = new CoarseIsotopePatternGenerator();
    auto maxIso = generator->estimateFromPeptideWeight(param.maxMass);
    maxIso.trimRight(0.01 * maxIso.getMostAbundant().getIntensity());
    param.maxIsotopeCount = (int) maxIso.size() - 1;
    generator->setMaxIsotope((Size) param.maxIsotopeCount);
    return PrecalcularedAveragine(100, param.maxMass, 50, generator);
  }

  static vector<PeakGroup> Deconvolution(MSExperiment &map, Parameter &param,
                                         PrecalcularedAveragine &averagines, int &specCntr, int &qspecCntr,
                                         int &massCntr)
  {
    double *filter = new double[param.chargeRange];
    //long **hBinOffsets = new long *[param.chargeRange];
    auto harmonicFilter = new double *[param.hCharges.size()];
    int *random = new int[param.chargeRange];
    std::srand(std::time(nullptr));

    for (int i = 0; i < param.chargeRange; i++)
    {
      //auto r = ((float) rand()) / (float) RAND_MAX - .5;
      filter[i] = log(
          1.0 / (i + param.minCharge)); //.124 * (1 + i%2) + should be descending, and negative!  - 0.2 * (1+(i%3)/2.0))

      //cout << i<< " " << filter[i]<<endl;
      /*hBinOffsets[i] = new long[param.hCharges.size()];
      for (Size k = 0; k < param.hCharges.size(); k++)
      {
        auto &hc = param.hCharges[k];
        float n = (float) (hc / 2);
        auto harmonicFilter = log(1.0 / (i + n / hc + param.minCharge));
        hBinOffsets[i][k] = (long) floor((filter[i] - harmonicFilter) * param.binWidth);
      }*/
    }

    for (Size k = 0; k < param.hCharges.size(); k++)
    {
      harmonicFilter[k] = new double[param.chargeRange];
      auto &hc = param.hCharges[k];
      float n = (float) (hc / 2);

      for (int i = 0; i < param.chargeRange; i++)
      {
        //auto r = ((float) rand()) / (float) RAND_MAX - .5;log(1.0 / (i + n / hc + param.minCharge));
        harmonicFilter[k][i] = log(
            1.0 / (i - n / hc +
                   param.minCharge)); //.124 * (1 + i%2) + should be descending, and negative!  - 0.2 * (1+(i%3)/2.0))
      }
    }

    if (param.jitter != 0)
    {
      double *tfilter = new double[param.chargeRange];
      auto m = filter[0];
      auto M = filter[param.chargeRange - 1];
      for (int i = 0; i < param.chargeRange; i++)
      {
        tfilter[i] = -filter[param.chargeRange - i - 1] + M + m;
      }
      filter = tfilter;
    }

    float prevProgress = .0;
    vector<PeakGroup> allPeakGroups;
    allPeakGroups.reserve(100000);
    //to overlap previous mass bins.
    vector<vector<Size>> prevMassBinVector;
    vector<double> prevMinBinLogMassVector;
    //vector <MSSpectrum> prevSpectra;

    for (auto it = map.begin(); it != map.end(); ++it)
    {
      if ((int) it->getMSLevel() > param.maxMSLevel)
      {
        continue;
      }

      //  if (it->getRT() < 3290 || it->getRT() > 3330)
      // {
      //  continue;
      // }
      /*double is = 0;
      for(auto&p : *it){
        is+=p.getIntensity();
      }

      cout<< it->getRT() <<" "<<is<<endl;
      continue;
 */
      float progress = (float) (it - map.begin()) / map.size();
      if (progress > prevProgress + .01)
      {
        printProgress(progress); //
        prevProgress = progress;
      }

      specCntr++;

      auto logMzPeaks = getLogMzPeaks(*it, param);
      if (logMzPeaks.empty())
      {
        continue;
      }
      auto peakGroups = getPeakGroupsFromSpectrum(logMzPeaks, filter, harmonicFilter,
                                                  prevMassBinVector, prevMinBinLogMassVector,
                                                  averagines,
                                                  param, specCntr);


      if (peakGroups.empty())
      {
        continue;
      }
      //vector<PeakGroup>().swap(peakGroups);
      //prevPgs = peakGroups;
      qspecCntr++;

      //allPeakGroups.reserve(allPeakGroups.size() + peakGroups.size());

      for (auto &pg : peakGroups)
      {
        massCntr++;
        pg.spec = &(*it);
        pg.massIndex = massCntr;
        pg.specIndex = qspecCntr;
        pg.massCntr = (int) peakGroups.size();
        allPeakGroups.push_back(pg);
      }
    }

    delete[] filter;

    for (Size k = 0; k < param.hCharges.size(); k++)
    {
      delete[] harmonicFilter[k];
    }
    delete[] harmonicFilter;
    delete[] random;
    printProgress(1); //
    //allPeakGroups.shrink_to_fit();
    return allPeakGroups; //
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

  static void printProgress(float progress)
  {
    //return;
    int barWidth = 70;
    cout << "[";
    int pos = (int) (barWidth * progress);
    for (int i = 0; i < barWidth; ++i)
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
    cout << "] " << int(progress * 100.0) << " %\r";
    cout.flush();
  }

  static vector<LogMzPeak> getLogMzPeaks(MSSpectrum &spec, const Parameter &param)
  {
    vector<LogMzPeak> logMzPeaks;
    logMzPeaks.reserve(spec.size());
    for (auto &peak: spec)
    {
      if (peak.getIntensity() <= param.intensityThreshold)
      {
        continue;
      }
      LogMzPeak logMzPeak(peak);
      logMzPeaks.push_back(logMzPeak);
    }

    //logMzPeaks.shrink_to_fit();

    return logMzPeaks;
  }


  static vector<PeakGroup>
  getPeakGroupsFromSpectrum(vector<LogMzPeak> &logMzPeaks, double *filter, double **harmonicFilter,
                            vector<vector<Size>> &prevMassBinVector,
                            vector<double> &prevMinBinLogMassVector,
                            PrecalcularedAveragine &averagines,
                            const Parameter &param, int &specCntr)
  {
    int sn = 1;
    double massDelta = (param.maxMass - param.minMass) / sn;

    double minMass = param.minMass + massDelta * (specCntr % sn);
    double maxMass = minMass + massDelta;

    double massBinMaxValue = min(
        logMzPeaks[logMzPeaks.size() - 1].logMz -
        filter[param.chargeRange - param.minContinuousChargePeakCount - 1],
        log(param.maxMass));

    double massBinMinValue = logMzPeaks[0].logMz - filter[param.minContinuousChargePeakCount];
    double mzBinMinValue = logMzPeaks[0].logMz;
    double mzBinMaxValue = logMzPeaks[logMzPeaks.size() - 1].logMz;
    Size massBinNumber = getBinNumber(massBinMaxValue, massBinMinValue, param.binWidth) + 1;

    long *binOffsets = new long[param.chargeRange];

    for (int i = 0; i < param.chargeRange; i++)
    {
      binOffsets[i] = (long) round((mzBinMinValue - filter[i] - massBinMinValue) *
                                   param.binWidth);//rand() %10000 + 100;
    }

    auto hBinOffsets = new long *[param.hCharges.size()];
    for (Size k = 0; k < param.hCharges.size(); k++)
    {
      hBinOffsets[k] = new long[param.chargeRange];
      for (int i = 0; i < param.chargeRange; i++)
      {
        hBinOffsets[k][i] = (long) round((mzBinMinValue - harmonicFilter[k][i] - massBinMinValue) *
                                         param.binWidth);
      }
    }


    if (param.jitter > 0)
    {

      for (int i = 0; i < param.chargeRange - 1; i++)
      {
        int diff = binOffsets[i + 1] - binOffsets[i];

        binOffsets[i] += rand() % diff;
      }

      binOffsets[param.chargeRange - 1] += rand() % 50;
      //random[i] = rand() % 400 - 200;
      //
      //param.jitter >0 ? random[i] : 0)
    }

    Size mzBinNumber = getBinNumber(mzBinMaxValue, mzBinMinValue, param.binWidth) + 1;
    float *intensities = new float[mzBinNumber];

    auto mzBins = getMzBins(logMzPeaks, mzBinMinValue, mzBinNumber, param.binWidth, intensities);

    boost::dynamic_bitset<> massBins(massBinNumber);
    float *sumLogIntensities = new float[massBinNumber];

    auto unionMassBins = getUnionMassBin(massBins, massBinMinValue, prevMassBinVector, prevMinBinLogMassVector,
                                         param);
    auto perMassChargeRanges = getMassBins(massBins, mzBins, massBinMinValue,
                                           sumLogIntensities,
                                           binOffsets,
                                           hBinOffsets,
                                           unionMassBins,
                                           intensities, param,
                                           minMass, maxMass);

    //cout<<1<<endl;
    auto peakGroups = getPeakGroupsWithMassBins(unionMassBins,
                                                logMzPeaks,
                                                mzBinMinValue, massBinMinValue, sumLogIntensities,
                                                binOffsets,
                                                perMassChargeRanges,
                                                param);

    //cout<<2<<endl;
    // cout<<peakGroups.size() << " "; //
    auto filteredPeakGroups = scoreAndFilterPeakGroups(peakGroups, averagines, param);
    if (prevMassBinVector.size() > 0 && prevMassBinVector.size() >= (Size) param.numOverlappedScans)
    {
      prevMassBinVector.erase(prevMassBinVector.begin());
      prevMinBinLogMassVector.erase(prevMinBinLogMassVector.begin());
    }
    //cout<<3<<endl;
    //cout<<filteredPeakGroups.size() << endl; //

    vector<Size> mb;
    mb.reserve(filteredPeakGroups.size());
    for (auto &pg : filteredPeakGroups)//filteredPeakGroups
    {
      pg.peaks.shrink_to_fit();
      if (massBins[pg.massBinIndex])
      {
        mb.push_back(pg.massBinIndex);
      }
    }

    prevMassBinVector.push_back(mb); //
    prevMinBinLogMassVector.push_back(massBinMinValue);

    prevMassBinVector.shrink_to_fit();
    prevMinBinLogMassVector.shrink_to_fit();

    //cout<<4<<endl;

    for (Size k = 0; k < param.hCharges.size(); k++)
    {
      delete[] hBinOffsets[k];
    }
    delete[] hBinOffsets;

    delete[] binOffsets;
    for (int i = 0; i < 3; i++)
    {
      delete[] perMassChargeRanges[i]; // delete array within matrix
    }// delete actual matrix
    delete[] perMassChargeRanges;
    delete[] intensities;
    delete[] sumLogIntensities;
    return filteredPeakGroups;
  }

  static boost::dynamic_bitset<> getUnionMassBin(boost::dynamic_bitset<> &massBins, double &massBinMinValue,
                                                 vector<vector<Size>> &prevMassBinVector,
                                                 vector<double> &prevMassBinMinValue,
                                                 const Parameter &param)
  {
    boost::dynamic_bitset<> u(massBins.size());
    if (u.size() == 0)
    {
      return u;
    }
    for (Size i = 0; i < prevMassBinVector.size(); i++)
    {
      auto &pmb = prevMassBinVector[i];
      if (pmb.empty())
      {
        continue;
      }
      long shift = (long) (round((massBinMinValue - prevMassBinMinValue[i]) * param.binWidth));

      for (auto &index : pmb)
      {
        long j = index - shift;
        if (j < 0)
        {
          continue;
        }
        if ((Size) j >= u.size())
        {
          break;
        }
        u[j] = true;
      }
    }
    return u;
  }

  static vector<PeakGroup> getPeakGroupsWithMassBins(boost::dynamic_bitset<> &unionedMassBins,
      //boost::dynamic_bitset<> &massBins,
                                                     vector<LogMzPeak> &logMzPeaks,
                                                     double &mzBinMinValue,
                                                     double &massBinMinValue,
                                                     float *sumLogIntensities,
                                                     long *binOffsets,
                                                     Byte **chargeRanges,
                                                     const Parameter &param)
  {
    double binWidth = param.binWidth;
    double tol = param.tolerance * 2;
    int minCharge = param.minCharge;
    int chargeRange = param.chargeRange;
    int maxIsotopeCount = param.maxIsotopeCount;

    int logMzPeakSize = (int) logMzPeaks.size();
    Size massBinSize = unionedMassBins.size();
    int *currentPeakIndex = new int[param.chargeRange];
    fill_n(currentPeakIndex, param.chargeRange, 0); //

    vector<PeakGroup> peakGroups;
    peakGroups.reserve(unionedMassBins.count());
    auto &minChargeRanges = chargeRanges[0];
    auto &maxChargeRanges = chargeRanges[1];
    auto &mzChargeRanges = chargeRanges[2];

    //vector<Size> toRemove;
    //toRemove.reserve(1000);
    //boost::dynamic_bitset<> tmp(logMzPeaks.size());

    auto massBinIndex = unionedMassBins.find_first();
    // Size lastSetMassBinIndex = unionedMassBins.size();
    Size *peakBinNumbers = new Size[logMzPeakSize];
    for (int i = 0; i < logMzPeakSize; i++)
    {
      peakBinNumbers[i] = getBinNumber(logMzPeaks[i].logMz, mzBinMinValue, binWidth);
    }

    while (massBinIndex != unionedMassBins.npos)
    {
      double logM = getBinValue(massBinIndex, massBinMinValue, binWidth);
      double diff = Constants::C13C12_MASSDIFF_U / exp(logM);
      double isoLogM1 = logM - diff;
      double isoLogM2 = logM + diff;

      auto b1 = getBinNumber(isoLogM1, massBinMinValue, binWidth);
      if (b1 > 0)
      {
        if (sumLogIntensities[massBinIndex] < sumLogIntensities[b1])
        {
          massBinIndex = unionedMassBins.find_next(massBinIndex);
          continue;
        }
      }

      auto b2 = getBinNumber(isoLogM2, massBinMinValue, binWidth);
      if (b2 < unionedMassBins.size())
      {
        if (sumLogIntensities[massBinIndex] < sumLogIntensities[b2])
        {
          massBinIndex = unionedMassBins.find_next(massBinIndex);
          continue;
        }
      }

      if (sumLogIntensities[b1] == 0 && sumLogIntensities[b2] == 0)
      {
        massBinIndex = unionedMassBins.find_next(massBinIndex);
        continue;
      }

      int isoOff = 0;
      PeakGroup pg;
      pg.reserve(chargeRange * 30);

      for (auto j = minChargeRanges[massBinIndex]; j <= maxChargeRanges[massBinIndex]; j++)
      {
        long &binOffset = binOffsets[j];
        auto bi = massBinIndex - binOffset;
        if (mzChargeRanges[bi] < chargeRange && mzChargeRanges[bi] != j)
        {
          continue;
        }

        int charge = j + minCharge;
        auto &cpi = currentPeakIndex[j];
        double maxIntensity = 0.0;
        int maxIntensityPeakIndex = -1;

        while (cpi < logMzPeakSize - 1)
        {
          //auto bi = peakBinNumbers[cpi] + binOffset;
          if (peakBinNumbers[cpi] == bi)
          {
            auto intensity = logMzPeaks[cpi].orgPeak->getIntensity();
            if (intensity > maxIntensity)
            {
              maxIntensity = intensity;
              maxIntensityPeakIndex = cpi;
            }
          }
          else if (peakBinNumbers[cpi] > bi)
          {
            break;
          }
          cpi++;
        }

        if (maxIntensityPeakIndex >= 0)
        {
          const double mz = logMzPeaks[maxIntensityPeakIndex].orgPeak->getMZ();
          //double &logMz = logMzPeaks[maxIntensityPeakIndex].logMz;
          const double isof = Constants::C13C12_MASSDIFF_U / charge;
          double mzDelta = tol * mz;
          //cout<<mzDelta<<endl;
          int maxI = 0;
          for (int d = -1; d <= 1; d += 2)
          { // negative then positive direction.
            int peakIndex = maxIntensityPeakIndex + (d < 0 ? d : 0);
            int lastPeakIndex = -100;
            for (int i = 0; i < maxIsotopeCount && (peakIndex >= 0 && peakIndex < logMzPeakSize); i++)
            {
              maxI = max(maxI, i);
              const double centerMz = mz + isof * i * d;
              const double centerMzMin = centerMz - mzDelta;
              const double centerMzMax = centerMz + mzDelta;
              bool isotopePeakPresent = false;
              if (lastPeakIndex >= 0)
              {
                peakIndex = lastPeakIndex;
              }//maxIntensityPeakIndex + (d < 0 ? d : 0);
              for (; peakIndex >= 0 && peakIndex < logMzPeakSize; peakIndex += d)
              {
                const double observedMz = logMzPeaks[peakIndex].orgPeak->getMZ();
                if (observedMz < centerMzMin)
                {
                  if (d < 0)
                  {
                    break;
                  }
                  else
                  {
                    continue;
                  }
                }
                if (observedMz > centerMzMax)
                {
                  if (d < 0)
                  {
                    continue;
                  }
                  else
                  {
                    break;
                  }
                }

                isotopePeakPresent = true;
                if (peakIndex != lastPeakIndex)
                {
                  // auto &mzBin = ;
                  const auto bin = peakBinNumbers[peakIndex] + binOffset;

                  // if (mzChargeRanges[peakBinNumbers[peakIndex]] < chargeRange && mzChargeRanges[peakBinNumbers[peakIndex]] != j)
                  // {
                  //  continue;
                  //}

                  if (bin < massBinSize) //
                  {
                    LogMzPeak p(*logMzPeaks[peakIndex].orgPeak, charge, i * d);
                    pg.peaks.push_back(p);
                    //isoOff = min(isoOff, p.isotopeIndex);
                    lastPeakIndex = peakIndex;
                    //}
                    // if (massBinIndex != bin && unionedMassBins[bin])
                    //{
                    //unionedMassBins[bin] = false;
                    //if (bin>0) unionedMassBins[bin-1] = false;
                    //if (bin < massBinSize -1)unionedMassBins[bin+1] = false;
                    // toRemove.push_back(bin);
                    //   }
                  }
                }
              }
              if (!isotopePeakPresent)
              {
                break;
              }
            }
          }


          //int minIsoIndex = maxI;
          for (auto &p : pg.peaks)
          {// assign the nearest isotope index..
            if (p.charge != charge)
            {
              continue;
            }
            for (int d = -1; d <= 1; d += 2)
            { // negative then positive direction.
              int maxId = 0;
              double minMzDelta = maxI;

              for (int i = 0; i <= maxI; i++)
              {
                double centerMz = mz + isof * (p.isotopeIndex + i * d);
                double delta = abs(centerMz - p.orgPeak->getMZ());

                if (delta > minMzDelta)
                {
                  break;
                }
                maxId = i * d;
                minMzDelta = delta;
              }
              //if (maxId != 0 )cout<< maxId <<endl;
              p.isotopeIndex += maxId;
            }
            isoOff = min(isoOff, p.isotopeIndex);
          }

        }
      }
      //cout<<pg.peaks.size()<<endl;
      //pg.peaks.shrink_to_fit();
      if (!pg.peaks.empty())
      {
        int minIi = 10000;
        int maxIi = -10000;
        for (auto &p : pg.peaks)
        {
          minIi = min(minIi, p.isotopeIndex);
          maxIi = max(maxIi, p.isotopeIndex);
          p.isotopeIndex -= isoOff;
        }
        if (minIi != maxIi)
        {

          //pg.updateMassesAndIntensity();
          pg.massBinIndex = massBinIndex;
          /*if(lastSetMassBinIndex == massBinIndex-1){// remove duplicate
              auto prevIntensity = peakGroups[peakGroups.size()-1].intensity;
              auto currentIntensity = pg.intensity;
              if(prevIntensity < currentIntensity){
                  peakGroups[peakGroups.size()-1] = pg;
              }
          }else*/
          peakGroups.push_back(pg);
        }
        //lastSetMassBinIndex = massBinIndex;
      }

      //if (massBinIndex < unionedMassBins.size() - 1 && !unionedMassBins[massBinIndex + 1])
      //{
      // for (auto &b : toRemove)
      // {
      //unionedMassBins[b] = false;
      //}
      //toRemove.clear();
      //}

      massBinIndex = unionedMassBins.find_next(massBinIndex);
    }
    delete[] currentPeakIndex;
    delete[] peakBinNumbers;
    return peakGroups;
  }

  static boost::dynamic_bitset<>
  getMzBins(vector<LogMzPeak> &logMzPeaks, double &mzBinMinValue, Size &binNumber, double binWidth,
            float *intensities
  )
  {
    boost::dynamic_bitset<> mzBins(binNumber);
    // double *intensities = new double[binNumber];
    fill_n(intensities, binNumber, .0);

    for (auto &p : logMzPeaks)
    {
      Size bi = getBinNumber(p.logMz, mzBinMinValue, binWidth);
      if (bi >= binNumber)
      {
        continue;
      }
      mzBins.set(bi);
      intensities[bi] += p.orgPeak->getIntensity();

      auto delta = (p.logMz - getBinValue(bi, mzBinMinValue, binWidth));

      if (delta > 0)
      {
        //add bin + 1
        if (bi < binNumber - 1)
        {
          mzBins.set(bi + 1);
          intensities[bi + 1] += p.orgPeak->getIntensity();
        }
      }
      else if (delta < 0)
      {
        if (bi > 0)
        {
          mzBins.set(bi - 1);
          intensities[bi - 1] += p.orgPeak->getIntensity();
        }
      }
    }


    //  delete[] intensities;
    return mzBins;
  }

  /*
  double calculateBinScoreThreshold(boost::dynamic_bitset<> &massBins, Byte *massBinScores, const Parameter &param){
      int *binScoreDist = new int[param.chargeRange + 1];
      fill_n(binScoreDist, param.chargeRange + 1, 0);

      auto setBinIndex = massBins.find_first();
      while (setBinIndex != massBins.npos) {
          binScoreDist[massBinScores[setBinIndex]]++;
          setBinIndex = massBins.find_next(setBinIndex);
      }
      int massCounter = 0;
      int binScoreThreshold = param.minContinuousChargePeakCount;
      for (int i = param.chargeRange; i >= 0; i--) {
          massCounter += binScoreDist[i];
          if (massCounter >= param.maxMassCount * 2) {
              binScoreThreshold = i > binScoreThreshold ? i : binScoreThreshold;
              break;
          }
      }
      return binScoreThreshold;
  }
 */
  /*
      int getMassBins(boost::dynamic_bitset<> &massBins, boost::dynamic_bitset<> &mzBins,
                      double &massBinMinValue, double &mzBinMinValue,
                      double *filter, const Parameter &param) {

          long *binOffsets = new long[param.chargeRange];
          for (int i = 0; i < param.chargeRange; i++) {
              binOffsets[i] = (long) (((mzBinMinValue - filter[i] - massBinMinValue) / param.tolerance) + .5);
          }

          //auto mzBinIndex = mzBins.find_first();
          int chargeRange = -1;

          //int minN2 = 3;
          int *minChargeCount = new int[param.chargeRange];
          // int* minContinuousChargePeakCount = new int[param.chargeRange];
          Size binNumber = massBins.size();
          Size binThresholdMinMass = getBinNumber(log(param.minMass), massBinMinValue, param.tolerance);
          double tol = param.tolerance;
          //

          Byte *numCharges = new Byte[binNumber];
          fill_n(numCharges, binNumber, 0);

          Byte *prevCharges = new Byte[binNumber];
          fill_n(prevCharges, binNumber, -1);

          Byte *chargeDiffSum = new Byte[binNumber];
          fill_n(chargeDiffSum, binNumber, -1);

          boost::dynamic_bitset<> massBins_(massBins.size());

          vector <boost::dynamic_bitset<>> massBinVector;
          massBinVector.reserve((Size) param.chargeRange);

          for (Byte j = 0; j < param.chargeRange; j++) {
              boost::dynamic_bitset<> mbin(binNumber);
              long shift = binOffsets[j];
              if(shift>0) mbin |= mzBins << (Size)shift;
              else mbin |= mzBins >> (Size)(-shift);

              mbin.resize(binNumber);
              if(mbin.count() > 0){
                  chargeRange = j;
              }

              massBinVector.push_back(mbin);
          }
          //cout<<chargeRange<<endl;

          int chargeRange_ = chargeRange;
          for (Byte j = 0; j <= chargeRange_; j++) {
              auto &mbin = massBinVector[j];
              auto mb = mbin.find_first();
              bool set = false;
              while (mb != mbin.npos) {
                  auto &nc = numCharges[mb];
                  nc++;
                  auto &p = prevCharges[mb];
                  chargeDiffSum[mb] += (j - p);
                  p = j;
                  if (nc >= param.minChargeCount) {
                      massBins_[mb] = true;
                      set = true;
                  }
                  mb = mbin.find_next(mb);
              }
              if(set) chargeRange = j;
          }


          auto mzBin = mzBins.find_first();
          while(mzBin != mzBins.npos) {
              long maxIndex = 0;
              int max = param.minChargeCount - 1;

              long minIndex = 0;
              float min = 100;

              for (Byte j = 0; j <= chargeRange; j++) {
                  auto mb = mzBin + binOffsets[j];
                  if (mb >= binNumber) break;
                  if(mb < binThresholdMinMass) continue;
                  //if (numCharges[mb] < 2) continue;
                  if (!massBins_[mb]) continue;

                  if (max < numCharges[mb]) {
                      max = numCharges[mb];
                      maxIndex = mb;
                  }

                  float cd = abs(((float) chargeDiffSum[mb] / (numCharges[mb]) - 1) - 1);

                  if (min > cd) {
                      min = cd;
                      minIndex = mb;
                  }
              }
              if (maxIndex > 0 && maxIndex == minIndex) {
                  massBins[maxIndex] = true;
              }
              mzBin = mzBins.find_next(mzBin);
          }

          return chargeRange + 1;// getFinalMassBinsAndDeleteHarmonicArtifacts(massBins, mzBins, mzBinMinValue, binOffsets,
          // harmonicBinOffsetVector, param);
      }*/


  static Byte **getMassBins(boost::dynamic_bitset<> &massBins, boost::dynamic_bitset<> &mzBins,
                            double &massBinMinValue,
                            float *sumLogIntensities,
                            long *binOffsets,
                            long **hBinOffsets,
                            boost::dynamic_bitset<> &unionMassBins,
                            float *intensities,
                            const Parameter &param, double &minMass, double &maxMass)
  {
    long binThresholdMinMass = (long) getBinNumber(log(minMass), massBinMinValue, param.binWidth);
    long binThresholdMaxMass = (long) min(massBins.size(),
                                          1 + getBinNumber(log(maxMass), massBinMinValue, param.binWidth));
    boost::dynamic_bitset<> isQualified(massBins.size());

    //cout<<minMass << " " << maxMass<<endl;
    //fill_n(sumLogIntensities, massBins.size(), 0);
    //cout<<0<<endl;
    getInitialMassBins(massBins, mzBins, isQualified,
                       sumLogIntensities,
        //continuousChargePeakPairCount,
                       hBinOffsets, binOffsets,
                       intensities,
                       param);
    //cout<<1<<endl;
    //printMasses(isQualified, massBinMinValue, continuousChargePeakPairCount, param);
    //auto toSkip = (isQualified | isQualified).flip();
    auto perMassChargeRanges = getFinalMassBins(massBins, mzBins, isQualified, unionMassBins,
                                                sumLogIntensities,
        //     massBinMinValue,
        // continuousChargePeakPairCount,
                                                binOffsets,
                                                param,
                                                binThresholdMinMass,
                                                binThresholdMaxMass);
    //cout<<2<<endl;
    //printMasses(massBins, massBinMinValue, continuousChargePeakPairCount, param);

    //cout<<isQualified.count()<<" " << massBins.count()<<" ";//
    //delete[] continuousChargePeakPairCount;
    //delete[] sumLogIntensities;
    return perMassChargeRanges;
  }

  static void
  printMasses(boost::dynamic_bitset<> &massBins, double &massBinMinValue, Byte *continuousChargePeakPairCount,
              const Parameter &param)
  {
    auto index = massBins.find_first();
    while (index != massBins.npos)
    {
      auto m = exp(getBinValue(index, massBinMinValue, param.binWidth));
      if (m < 50500 && m > 50400)
      {

        cout << m <<
             " " << (int) continuousChargePeakPairCount[index] <<
             //" " << (int) noneContinuousChargePeakPairCount[index] <<
             endl;
      }
      index = massBins.find_next(index);
    }
    //cout << endl;
  }


  static void getInitialMassBins(boost::dynamic_bitset<> &massBins,
                                 boost::dynamic_bitset<> &mzBins,
                                 boost::dynamic_bitset<> &isQualified,
                                 float *signal,
                                 long **hBinOffsets,
                                 long *binOffsets,
                                 float *intensities,
                                 const Parameter &param
  )
  {
    int chargeRange = param.chargeRange;
    int hChargeSize = (int) param.hCharges.size();
    int minContinuousChargePeakCount = param.minContinuousChargePeakCount;
    //    long mzBinSize = (long) mzBins.size();
    long binEnd = (long) massBins.size();

    Byte *continuousChargePeakPairCount = new Byte[massBins.size()];
    fill_n(continuousChargePeakPairCount, massBins.size(), 0);

    Byte *prevCharges = new Byte[massBins.size()];
    fill_n(prevCharges, massBins.size(), (Byte) (chargeRange + 2));

    float *prevIntensities = new float[massBins.size()];
    fill_n(prevIntensities, massBins.size(), 1);

    //auto prevHcharges = new Byte*[hChargeSize];
    //for(auto k=0;k<hChargeSize;k++){
    //  prevHcharges[k] = new Byte[massBins.size()];
    //  fill_n(prevHcharges[k], massBins.size(), (Byte) (chargeRange + 2));
    // }

    //float* signal = new float[massBins.size()];
    fill_n(signal, massBins.size(), 0);
    auto noise = new float *[hChargeSize + 1];
    for (auto k = 0; k < hChargeSize + 1; k++)
    {
      noise[k] = new float[massBins.size()];
      fill_n(noise[k], massBins.size(), 0);
    }

    //double th = 1;//
    float factor = 4.0;
    auto mzBinIndex = mzBins.find_first();
    while (mzBinIndex != mzBins.npos)
    {
      auto &intensity = intensities[mzBinIndex];
      for (Byte j = 0; j < chargeRange; j++)
      {
        long massBinIndex = mzBinIndex + binOffsets[j];
        if (massBinIndex < 0)
        {
          continue;
        }
        if (massBinIndex >= binEnd)
        {
          break;
        }
        auto cd = prevCharges[massBinIndex] - j;

        auto &prevIntensity = prevIntensities[massBinIndex];
        float minInt = min(intensity, prevIntensity);
        float maxInt = max(intensity, prevIntensity);

        float id = maxInt / minInt;

        if (prevCharges[massBinIndex] < chargeRange && cd != 1 && id<factor)
        {
          noise[hChargeSize][massBinIndex] += minInt;
        }

        if (cd != 1 || id > factor)
        {
          continuousChargePeakPairCount[massBinIndex] = 0;
        }
        else
        {
          int maxHcharge = -1;
          float maxHint = .0;
          for (auto k = 0; k < hChargeSize; k++)
          {
            auto hmzBinIndex = massBinIndex - hBinOffsets[k][j];
            if (hmzBinIndex > 0 && hmzBinIndex < (long) mzBins.size() && mzBins[hmzBinIndex])
            {
              auto &hintensity = intensities[hmzBinIndex];
              if (hintensity > minInt
                  &&  hintensity < factor * maxInt
                  //&&  hintensity > maxInt/factor
                  )
              {
                noise[k][massBinIndex] += hintensity;

                if (hintensity < maxHint)
                {
                  continue;
                }
                maxHint = hintensity;
                maxHcharge = k;
              }
            }
          }
          if (maxHcharge >= 0)
          {
            //  cout<<1<<endl;
            continuousChargePeakPairCount[massBinIndex] = 0;
          }
          else
          {
            signal[massBinIndex] += intensity;
            if (!isQualified[massBinIndex])
            {
              isQualified[massBinIndex] =
                  ++continuousChargePeakPairCount[massBinIndex] >= minContinuousChargePeakCount;
            }
          }
        }
        prevIntensity = intensity;
        prevCharges[massBinIndex] = j;
      }
      mzBinIndex = mzBins.find_next(mzBinIndex);
    }

    auto mindex = isQualified.find_first();
    while (mindex != isQualified.npos)
    {
      auto &s = signal[mindex];
      // auto msnr = s / (noise[0][mindex]);
      float maxNoise = .0;
      for (auto k = 0; k < hChargeSize + 1; k++)
      {
        maxNoise = max(maxNoise, noise[k][mindex]);
        // msnr = min(snr, msnr);
      }
      //s -= maxNoise;
      s -= maxNoise;
      mindex = isQualified.find_next(mindex);
    }

    delete[] prevIntensities;
    delete[] continuousChargePeakPairCount;
    delete[] prevCharges;
    for (auto k = 0; k < hChargeSize + 1; k++)
    {
      delete[] noise[k];
    }
    //for(auto k=0;k<hChargeSize; k++)
    //{
    //  delete[] prevHcharges[k];
    //}
    //delete[] prevHcharges;
    delete[] noise;
    //delete[] signal;

  }

  static Byte **getFinalMassBins(boost::dynamic_bitset<> &massBins, boost::dynamic_bitset<> &mzBins,
                                 boost::dynamic_bitset<> &isQualified,
                                 boost::dynamic_bitset<> &unionMassBins,
                                 float *sumLogIntensities,
      // double &massBinMinValue,
                                 long *binOffsets,
                                 const Parameter &param,
                                 long &binStart, long &binEnd
  )
  {
    int chargeRange = param.chargeRange;
    // auto binWidth = param.binWidth;
    Byte *maxChargeRanges = new Byte[massBins.size()];
    fill_n(maxChargeRanges, massBins.size(), 0);

    Byte *minChargeRanges = new Byte[massBins.size()];
    fill_n(minChargeRanges, massBins.size(), chargeRange + 1);

    Byte *mzChargeRanges = new Byte[mzBins.size()];
    fill_n(mzChargeRanges, mzBins.size(), chargeRange + 1);
    //Byte *selected = new Byte[massBins.size()];
    //fill_n(selected, massBins.size(), 0);

    auto mzBinIndex = mzBins.find_first();
    long binSize = (long) massBins.size();
    //int minContinuousChargePeakCount = param.minContinuousChargePeakCount;

    auto toSkip = (isQualified | unionMassBins).flip();
    unionMassBins.reset();


    while (mzBinIndex != mzBins.npos)
    {
      long maxIndex = -1;
      float maxCount = -1e11;
      Byte charge = 0;

      //vector<Byte> setCharges;
      //setCharges.reserve(chargeRange);
      for (Byte j = 0; j < chargeRange; j++)
      {
        long massBinIndex = mzBinIndex + binOffsets[j];
        if (massBinIndex < 0)
        {
          continue;
        }
        if (massBinIndex >= binSize)
        {
          break;
        }
        if (toSkip[massBinIndex])
        {
          continue;
        }

        //if (continuousChargePeakPairCount[massBinIndex] < round(log2(j+param.minCharge)))
        //  continue;

        auto &t = sumLogIntensities[massBinIndex];// + noneContinuousChargePeakPairCount[massBinIndex];//
        if (t == 0) { // no signal
          continue;
        }
        if (maxCount < t)
        {
          maxCount = t;
          maxIndex = massBinIndex;
          charge = j;
        }
      }

      // if (maxIndex >= 0)
      //{
      if (maxIndex > binStart && maxIndex < binEnd)
      {
        //if (sumLogIntensities[maxIndex] >  sumLogIntensities[maxIndex-1] &&
        //    sumLogIntensities[maxIndex] >  sumLogIntensities[maxIndex+1])
        {
          maxChargeRanges[maxIndex] = max(maxChargeRanges[maxIndex], charge);
          minChargeRanges[maxIndex] = min(minChargeRanges[maxIndex], charge);
          //  ++selected[maxIndex];
          //int t = max(continuousChargePeakPairCount[maxIndex],minContinuousChargePeakCount);
          // if (++selected[maxIndex] >= 3)
          //&& selected[maxIndex] >= continuousChargePeakPairCount[maxIndex])
          //{
          //massBins[maxIndex] = isQualified[maxIndex];

          massBins[maxIndex] = isQualified[maxIndex];

          mzChargeRanges[mzBinIndex] = charge;//minChargeRanges[maxIndex];//...

          unionMassBins[maxIndex] = true;
        }
        // }
      }
      //}
      mzBinIndex = mzBins.find_next(mzBinIndex);
    }


    // bool isChargeWellDistributed = checkChargeDistribution(perChargeIntensity,
    //                                                       param.chargeRange,
    //                                                      3);

    //cout<<" " <<massBins.count()<<endl; //
    Byte **chargeRanges = new Byte *[3];
    chargeRanges[0] = minChargeRanges;
    chargeRanges[1] = maxChargeRanges;
    chargeRanges[2] = mzChargeRanges;
    // delete[] selected;
    return chargeRanges;
  }

  /*
  static double
  getIsotopeCosineThreshold(double mass, double minCosine, double maxCosine, double minMass, double maxMass)
  {
    double isotopeCosineThreshold;
    if (mass < minMass)
    {
      isotopeCosineThreshold = maxCosine;
    }
    else if (mass > maxMass)
    {
      isotopeCosineThreshold = minCosine;
    }
    else
    {
      isotopeCosineThreshold = minCosine + (maxCosine - minCosine) / (maxMass - minMass) * (maxMass - mass);
    }
    //cout<<mass<<" "<<minCosine<< " " << maxCosine<< " " <<  isotopeCosineThreshold<<endl;
    return isotopeCosineThreshold;
  }
*/
  static vector<PeakGroup> scoreAndFilterPeakGroups(vector<PeakGroup> &peakGroups,
                                                    PrecalcularedAveragine &averagines,
                                                    const Parameter &param)
  {
    vector<PeakGroup> filteredPeakGroups;
    filteredPeakGroups.reserve(peakGroups.size());
    double threshold = .0;

    Size mc = (Size) param.maxMassCount;
    if (mc > 0)
    {
      vector<double> intensities;
      intensities.reserve(peakGroups.size());

      for (auto &pg : peakGroups)
      {
        pg.updateMassesAndIntensity(averagines);
        intensities.push_back(pg.intensity);
      }

      if (intensities.size() > mc)
      {
        sort(intensities.begin(), intensities.end());
        threshold = intensities[intensities.size() - mc];
      }
      vector<double>().swap(intensities);
    }
    //filterPeakGroupsByIntensity(peakGroups, intensities, param);


    //auto perIsotopeMaxCharge = new int[param.maxIsotopeCount];
    //auto perIsotopeMinCharge = new int[param.maxIsotopeCount];
    //auto perChargeMaxIsotope = new int[param.chargeRange];
    //auto perChargeMinIsotope = new int[param.chargeRange];// polish later

    auto perIsotopeIntensity = new double[param.maxIsotopeCount];
    auto perChargeIntensity = new double[param.chargeRange];

    for (auto &pg : peakGroups)
    {
      if (pg.intensity < threshold)
      {
        continue;
      }
      updatePerChargeIsotopeIntensity(//perIsotopeMinCharge, perIsotopeMaxCharge,
          //perChargeMinIsotope, perChargeMaxIsotope,
          perIsotopeIntensity, perChargeIntensity,
          pg, param);

      pg.chargeCosineScore = getChargeFitScore(perChargeIntensity, param.chargeRange);

      if (pg.peaks.empty() || (pg.chargeCosineScore < param.minChargeCosineSpec))
      {
        continue;
      }

      /*bool isIsotopeSpanGood = checkSpanDistribution(perChargeMinIsotope, perChargeMaxIsotope,
                                                     param.chargeRange, 3);

      if (!isIsotopeSpanGood)
      {
        continue;
      }


      bool isChargeSpanGood = checkSpanDistribution(perIsotopeMinCharge, perIsotopeMaxCharge,
                                                    param.maxIsotopeCount, 3);

      if (!isChargeSpanGood)
      {
        continue;
      }
*/

      bool isChargeWellDistributed = checkChargeDistribution(perChargeIntensity,
                                                             param.chargeRange,
                                                             param.minContinuousChargePeakCount);

      if (!isChargeWellDistributed)
      {
        continue;
      }

      int offset = 0;
      pg.isotopeCosineScore = getIsotopeCosineAndDetermineIsotopeIndex(pg.peaks[0].getMass(),
                                                                       perIsotopeIntensity,
                                                                       param.maxIsotopeCount,
                                                                       averagines, offset);

      //double isotopeCosineThreshold = param.minIsotopeCosineSpec;// getIsotopeCosineThreshold(pg.peaks[0].getMass(),
      //   param.minIsotopeCosine, param.maxIsotopeCosine,
      //  0, 100000);

      if (pg.peaks.empty() || (pg.isotopeCosineScore < param.minIsotopeCosineSpec))
      {
        continue;
      }


      pg.updateMassesAndIntensity(averagines, offset, param.maxIsotopeCount);
      //cout<<"mass "<< pg.monoisotopicMass <<endl;


      /*updatePerChargeIsotopeIntensity(perIsotopeMinCharge, perIsotopeMaxCharge,
                                      perChargeMinIsotope, perChargeMaxIsotope,
                                      perIsotopeIntensity, perChargeIntensity,
                                      pg, param);
      if (pg.peaks[0].mass>60000){
        cout<<"iso=[";
        for(int i=0;i<param.maxIsotopeCount;i++){
          //if (perIsotopeIntensity[i]<=0){
          //  break;
          //}
          cout<<perIsotopeIntensity[i]<<" ";
        }
        cout<<"];";

        averagines.print(pg.peaks[0].mass);
      }
 */
      filteredPeakGroups.push_back(pg);

    }

    vector<PeakGroup>().swap(peakGroups);

    removeOverlappingPeakGroups(filteredPeakGroups, param.tolerance);
    //removeHarmonicPeakGroups(filteredPeakGroups, param);

    delete[] perIsotopeIntensity;
    delete[] perChargeIntensity;
    //delete[] perIsotopeMinCharge;
    //delete[] perIsotopeMaxCharge;
    //delete[] perChargeMinIsotope;
    //delete[] perChargeMaxIsotope;


    //filteredPeakGroups.shrink_to_fit();
    return filteredPeakGroups;
  }


  static void removeOverlappingPeakGroups(vector<PeakGroup> &pgs, double tol)
  { // pgs are sorted
    vector<PeakGroup> merged;
    merged.reserve(pgs.size());

    for (Size i = 0; i < pgs.size(); i++)
    {
      bool select = true;
      auto &pg = pgs[i];
      double massTol = pg.monoisotopicMass * tol * 2;

      for (Size j = i + 1; j < pgs.size(); j++)
      {
        auto &pgo = pgs[j];
        if (!select || pgo.monoisotopicMass - pg.monoisotopicMass > massTol)
        {
          break;
        }
        select &= pg.intensity > pgo.intensity;
      }

      if (!select)
      {
        continue;
      }
      for (int j = i - 1; j >= 0; j--)
      {
        auto &pgo = pgs[j];
        if (!select || pg.monoisotopicMass - pgo.monoisotopicMass > massTol)
        {
          break;
        }
        select &= pg.intensity > pgo.intensity;
      }
      if (!select)
      {
        continue;
      }
      merged.push_back(pg);
    }
    //pgs = merged;
    //merged.shrink_to_fit();
    vector<PeakGroup>().swap(pgs);
    merged.swap(pgs);
    //vector<PeakGroup>().swap(merged);
  }

  static void updatePerChargeIsotopeIntensity(//int *perIsotopeMinCharge, int *perIsotopeMaxCharge,
      //int *perChargeMinIsotope, int *perChargeMaxIsotope,
      double *perIsotopeIntensity,
      double *perChargeIntensity,
      PeakGroup &pg,
      const Parameter &param)
  {


    //fill_n(perIsotopeMinCharge, param.maxIsotopeCount, param.chargeRange);
    //fill_n(perIsotopeMaxCharge, param.maxIsotopeCount, -1);

    //fill_n(perChargeMinIsotope, param.chargeRange, param.maxIsotopeCount);
    //fill_n(perChargeMaxIsotope, param.chargeRange, -1);

    fill_n(perIsotopeIntensity, param.maxIsotopeCount, 0);
    fill_n(perChargeIntensity, param.chargeRange, 0);

    // bool *tmp = new bool[param.chargeRange * param.maxIsotopeCount];
    int minCharge = param.chargeRange + param.minCharge + 1;
    int maxCharge = 0;

    for (auto &p : pg.peaks)
    {
      if (p.isotopeIndex < 0 || p.isotopeIndex >= param.maxIsotopeCount)
      {
        continue;
      }
      minCharge = min(minCharge, p.charge);
      maxCharge = max(maxCharge, p.charge);

      int index = p.charge - param.minCharge;
      perIsotopeIntensity[p.isotopeIndex] += p.orgPeak->getIntensity();
      perChargeIntensity[index] += p.orgPeak->getIntensity();

      //perChargeMinIsotope[index] = min(perChargeMinIsotope[index], p.isotopeIndex);
      //perChargeMaxIsotope[index] = max(perChargeMaxIsotope[index], p.isotopeIndex);

      //perIsotopeMinCharge[p.isotopeIndex] = min(perIsotopeMinCharge[p.isotopeIndex], index);
      //perIsotopeMaxCharge[p.isotopeIndex] = max(perIsotopeMaxCharge[p.isotopeIndex], index);
      //  tmp[p.isotopeIndex + param.maxIsotopeCount * index]=true;
    }
    pg.maxCharge = maxCharge;
    pg.minCharge = minCharge;

  }


  static double getIsotopeCosineAndDetermineIsotopeIndex(double mass,
                                                         double *perIsotopeIntensities,
                                                         int perIsotopeIntensitiesSize,
                                                         PrecalcularedAveragine &averagines,
                                                         int &offset)
  {
    auto iso = averagines.get(mass);
    auto isoNorm = averagines.getNorm(mass);

    int isoSize = (int) iso.size();

    //double isoDiff = Constants::C13C12_MASSDIFF_U;// OpenMS::Math::mean(diffs.begin(), diffs.end());

    offset = 0;
    double maxCosine = -1;
    int maxIsotopeIndex = 0, minIsotopeIndex = -1;

    for (int i = 0; i < perIsotopeIntensitiesSize; i++)
    {
      if (perIsotopeIntensities[i] <= 0)
      {
        continue;
      }
      maxIsotopeIndex = i;
      if (minIsotopeIndex < 0)
      {
        minIsotopeIndex = i;
      }
    }

    for (int f = -isoSize + minIsotopeIndex; f <= maxIsotopeIndex; f++)
    {
      auto cos = getCosine(perIsotopeIntensities, minIsotopeIndex, maxIsotopeIndex, iso, isoSize, isoNorm, f);

      if (maxCosine <= cos)
      {
        maxCosine = cos;
        offset = f;
      }
    }

    return maxCosine;
  }


  static bool checkSpanDistribution(int *mins, int *maxs, int range, int threshold)
  {
    int nonZeroStart = -1, nonZeroEnd = 0;
    int maxSpan = 0;

    for (int i = 0; i < range; i++)
    {
      if (maxs[i] >= 0)
      {
        if (nonZeroStart < 0)
        {
          nonZeroStart = i;
        }
        nonZeroEnd = i;
        maxSpan = max(maxSpan, maxs[i] - mins[i]);
      }
    }
    if (maxSpan <= 0)
    {
      return false;
    }

    int prevCharge = nonZeroStart;
    int n_r = 0;

    double spanThreshold = maxSpan / 1.5;//

    for (int k = nonZeroStart + 1; k <= nonZeroEnd; k++)
    {
      if (maxs[k] < 0)
      {
        continue;
      }


      if (k - prevCharge == 1)
      {
        int intersectSpan = min(maxs[prevCharge], maxs[k])
                            - max(mins[prevCharge], mins[k]);

        if (spanThreshold <= intersectSpan)
        { //
          n_r++;
        }
      }
      prevCharge = k;
    }

    if (n_r < threshold)
    {
      return -100.0;
    }

    for (int i = 2; i < min(12, range); i++)
    {
      for (int l = 0; l < i; l++)
      {
        int t = 0;
        prevCharge = nonZeroStart + l;
        for (int k = prevCharge + i; k <= nonZeroEnd; k += i)
        {
          if (maxs[k] < 0)
          {
            continue;
          }


          if (k - prevCharge == i)
          {
            int intersectSpan = min(maxs[prevCharge], maxs[k])
                                - max(mins[prevCharge], mins[k]);

            if (spanThreshold <= intersectSpan)
            {
              t++;
            }
          }
          prevCharge = k;
        }
        if (n_r <= t)
        {
          return false;
        }
      }
    }

    return true;
  }

  static double getChargeFitScore(double *perChargeIntensity,
                                  int range)
  {
    double maxPerChargeIntensity = .0;
    vector<double> xs;
    vector<double> ys;

    xs.reserve(range + 2);
    ys.reserve(range + 2);

    for (int i = 0; i < range; i++)
    {
      maxPerChargeIntensity = max(maxPerChargeIntensity, perChargeIntensity[i]);
    }

    double th = maxPerChargeIntensity * .02;//2%
    int first = -1, last = 0;
    for (int i = 0; i < range; i++)
    {
      if (perChargeIntensity[i] <= th)
      {
        continue;
      }
      if (first < 0)
      {
        first = i;
        // xs.push_back(first-1000);
        //ys.push_back(minInt/10);
      }

      last = i;
      //cout<<log(perChargeIntensity[i])<<endl;
      //  xs.push_back(i);
      // ys.push_back((1+perChargeIntensity[i]));
    }
    if (last - first < 2)
    {
      return 0;
    }
    //xs.push_back(first-2);
    //ys.push_back(1);
    //xs.push_back(first-1);
    //ys.push_back(1);

    for (int i = first; i <= last; i++)
    {
      //if (perChargeIntensity[i]<=0)
      //{
      //  continue;
      //}
      xs.push_back(i);
      ys.push_back((1 + perChargeIntensity[i]));
    }
    //xs.push_back(last+1);
    //ys.push_back(1);
    //xs.push_back(last+2);
    //ys.push_back(1);


    // xs.push_back(last + 1000);
    // ys.push_back((minInt/10));

    Eigen::Matrix3d m;
    Eigen::Vector3d v;

    double s0 = 0, s1 = 0, s2 = 0, s3 = 0, s4 = 0;
    double t0 = 0, t1 = 0, t2 = 0;

    for (Size i = 0; i < xs.size(); i++)
    {
      auto &x = xs[i];
      auto y = log(ys[i]);
      s0++;
      s1 += x;
      s2 += x * x;
      s3 += x * x * x;
      s4 += x * x * x * x;
      t0 += y;
      t1 += y * x;
      t2 += y * x * x;
    }
    m(0, 0) = s0;
    m(1, 0) = m(0, 1) = s1;
    m(2, 0) = m(1, 1) = m(0, 2) = s2;
    m(2, 1) = m(1, 2) = s3;
    m(2, 2) = s4;

    //cout<<m<<endl;
    //cout<<v<<endl;
    auto im = m.inverse();
    //cout<<im<<endl;
    // cout<<m<<endl;
    v(0) = t0;
    v(1) = t1;
    v(2) = t2;
    //cout<<v<<endl;
    auto abc = im * v;
    //cout<<abc<<endl;
    double mu = -abc(1) / abc(2) / 2;
    double omega = -1 / abc(2) / 2;
    //double a = (abc(0)-abc(1)*abc(1))/4/abc(2);

    if (omega <= 0)
    {
      //omega = -omega;
      //mu = -mu;
      return 0;
    }
    //cout<<"**"<<endl;
    //cout<<mu <<" * " << omega<<endl;
    vector<double> tys;

    //cout<<endl;
    for (Size i = 0; i < ys.size(); i++)
    {
      double ty = exp(-(xs[i] - mu) * (xs[i] - mu) / 2 / omega);
      tys.push_back(ty);
      //  cout<< xs[i]<<" "<<ys[i] << " " << ty <<endl;
    }
    // cout<<endl;
    // cout<<"cos " << getCosine(ys, tys) <<endl;
    return getCosine(ys, tys);

  }


  static bool checkChargeDistribution(double *perChargeIntensity,
                                      int range,
                                      int threshold)
  {
    double maxPerChargeIntensity = .0;
    int nonZeroStart = -1, nonZeroEnd = 0;
    for (int i = 0; i < range; i++)
    {
      if (perChargeIntensity[i] > 0)
      {
        // intensities.push_back(intensity);
        maxPerChargeIntensity = max(maxPerChargeIntensity, perChargeIntensity[i]);
        if (nonZeroStart < 0)
        {
          nonZeroStart = i;
        }
        nonZeroEnd = i;
      }
    }

    int prevCharge = nonZeroStart;

    int n_r = .0;

    double intensityThreshold = maxPerChargeIntensity / 4.0;//intensities[intensities.size()*95/100] / 5.0;
    for (int k = prevCharge + 1; k <= nonZeroEnd; k++)
    {
      if (perChargeIntensity[k] <= intensityThreshold)
      {
        continue;
      }

      if (k - prevCharge == 1)
      {
        n_r++;
      }
      prevCharge = k;
    }

    if (n_r < threshold)
    {
      return false;
    }


    for (int i = 2; i < min(7, range); i++)
    {
      for (int l = 0; l < i; l++)
      {
        int t = 0;
        prevCharge = nonZeroStart + l;
        for (int k = prevCharge + i; k <= nonZeroEnd; k += i)
        {
          if (perChargeIntensity[k] <= intensityThreshold)
          { //
            continue;
          }
          if (k - prevCharge == i)
          {
            t++;
          }
          prevCharge = k;
        }
        if (n_r <= t)
        {
          return false;
        }
      }
    }

    return true;
  }

  static double
  getCosine(double *a, int &aStart, int &aEnd, IsotopeDistribution &b, int &bSize, double &bNorm, int offset = 1)
  {
    double n = .0, d1 = .0;

    for (int j = aStart; j < aEnd; j++)
    {
      d1 += a[j] * a[j];
      int i = j - offset;
      if (i < 0 || i >= bSize)
      {
        continue;
      }
      n += a[j] * b[i].getIntensity(); //
    }

    double d = (d1 * bNorm);
    if (d <= 0)
    {
      return 0;
    }
    return n / sqrt(d);
  }

  static double
  getCosine(vector<double> &a, vector<double> &b, int off = 0)
  {
    double n = .0, d1 = .0, d2 = .0;
    Size size = a.size();
    //int overlapCntr = 0;
    for (Size j = off; j < size - off; j++)
    {
      d1 += a[j] * a[j];
      d2 += b[j] * b[j];
      n += a[j] * b[j];
      //  if(a[j] > 0 && b[j] > 0) overlapCntr++;
    }

    //if(overlapCntr < 2) return 0; //
    double d = (d1 * d2);
    if (d <= 0)
    {
      return 0;
    }
    return n / sqrt(d);
  }

  static void
  filterPeakGroupsByIntensity(vector<PeakGroup> &peakGroups, vector<double> &intensities, const Parameter &param)
  {
    if (param.maxMassCount < 0 || intensities.size() <= (Size) param.maxMassCount)
    {
      return;
    }
    Size mc = (Size) param.maxMassCount;
    sort(intensities.begin(), intensities.end());
    auto threshold = intensities[intensities.size() - mc];
    for (auto pg = peakGroups.begin(); pg != peakGroups.end();)
    {
      if (peakGroups.size() <= mc)
      {
        break;
      }
      if (pg->intensity < threshold)
      {
        pg = peakGroups.erase(pg);
        continue;
      }
      ++pg;
    }
  }


  /*
      void findNominalMassFeatures(vector<PeakGroup> &peakGroups, int &id, fstream& fs, const Parameter & param){

          map<int, vector<double>> massMap;
          map<int, vector<String>> writeMap;

          double delta = param.RTwindow;

          //sort(peakGroups.begin(), peakGroups.end());
          for(auto& pg : peakGroups){
              int nm = getNominalMass(pg.monoisotopicMass);
              double rt = pg.spec->getRT();
              double &intensity = pg.intensity;
              auto it = massMap.find(nm);
              vector<double> v; // start rt, end rt, apex rt, intensity, maxIntensity, abundance
              if (it == massMap.end()){
                  v.push_back(rt); // start rt 0
                  v.push_back(rt); // end rt 1
                  v.push_back(rt); // apex rt 2
                  v.push_back(intensity); // max intensity 3
                  v.push_back(0.0); // abundance 4
                  massMap[nm] = v;
                  continue;
              }
              v = it->second;
              double prt = v[1];
              double pMaxIntensity = v[3];

              if(rt - prt <delta){ // keep feature
                  v[1] = rt;
                  v[4] += (pMaxIntensity + intensity)/2.0*(rt - prt);
                  if(pMaxIntensity < intensity) {
                      v[2] = rt;
                      v[3] = intensity;
                  }
                  massMap[nm] = v;
              }else{ // write feature and clear key
                  if(v[1]-v[0]>0) {
                      stringstream s;
                      s << "\t" << param.fileName << "\t" << nm << "\t" << fixed << setprecision(5) << v[0] << "\t"
                        << v[1] << "\t" << v[1] - v[0] << "\t" << v[2] << "\t" << v[3] << "\t" << v[4];
                      auto sit = writeMap.find(nm);
                      if (sit == writeMap.end()) {
                          vector <String> vs;
                          writeMap[nm] = vs;
                      }
                      writeMap[nm].push_back(s.str());
                  }
                  massMap.erase(it);
              }
          }

          for (auto it=massMap.begin(); it!=massMap.end(); ++it){
              int nm = it->first;
              auto v = it->second;
              if(v[1] - v[0] <= 0) continue;
              stringstream s;
              s<<"\t"<< param.fileName << "\t" <<nm<<"\t"<<fixed<<setprecision(5)<<v[0]<<"\t"<<v[1]<<"\t"<<v[1]-v[0]<<"\t"
               <<v[2]<<"\t"<<v[3]<< "\t" << v[4];
              auto sit = writeMap.find(nm);
              if (sit == writeMap.end()) {
                  vector<String> vs;
                  writeMap[nm] = vs;
              }
              writeMap[nm].push_back(s.str());
          }

          for (auto it=writeMap.begin(); it!=writeMap.end(); ++it) {
              for(auto& s : it->second){
                  fs<<++id<<s<<"\n";
              }
          }
      }*/
};

int main(int argc, const char **argv)
{
  FLASHDeconv tool;
  return tool.main(argc, argv);
}

/*
   double getLogLikelihood(double int1, double int2, bool isH0){
       double tmp = 1e-1;
       double ret;
       if(int1<=0){
           if(int2<=0){
               if(isH0) ret = .8;
               else
                   ret = .9;
           }else{
               if(isH0) ret = .2;
               else
                   ret = .1;
           }
       }else{
           if(int2<=0){
               if(isH0) ret = .8;
               else ret = .4;
           }else if(int1<int2){
               if(isH0) ret = .1;
               else ret = .1;
           }else {
               if (isH0) ret = .1;
               else ret = .5;
           }
       }
       return log10(ret);
   }

*/

/*
vector<Peak1D> filterSpectrum(MSSpectrum &spectrum, double window, double factor, double th){
   deque<Peak1D> peaksInWindow;
   vector<Peak1D> filtered;
   vector<double> intensityHeap;
   intensityHeap.reserve(spectrum.size());
   filtered.reserve(spectrum.size());
   int wsIndex = 0, weIndex = 0;
   double w = window/2;
   double prevMedian = 0;
   vector<Peak1D> initFiltered;
   initFiltered.reserve(spectrum.size());
   for(auto &p : spectrum) if(p.getIntensity() > th) initFiltered.push_back(p);

   for (int i=0;i<initFiltered.size();i++) {
       auto p = initFiltered[i];
       auto mz = p.getMZ();
       double median = prevMedian;
       while(initFiltered[wsIndex].getMZ() < mz - w){
           auto firstp = peaksInWindow.front();
           auto j = lower_bound(intensityHeap.begin(), intensityHeap.end(), firstp.getIntensity());
           intensityHeap.erase(j);
           // find firstp in heap and remove it using binary search
           peaksInWindow.pop_front();
           wsIndex++;
       }
       while(weIndex< initFiltered.size() && initFiltered[weIndex].getMZ() < mz + w){
           auto lastp = spectrum[weIndex];
           peaksInWindow.push_back(lastp);
           auto j = lower_bound(intensityHeap.begin(), intensityHeap.end(), lastp.getIntensity());
           intensityHeap.insert(j, lastp.getIntensity());
           median = intensityHeap[intensityHeap.size()/2];
           weIndex++;
       }
       if(p.getIntensity() >= median * factor)
           filtered.push_back(p);

       prevMedian = median;
   }
   return filtered;
}*/



/*
   void sortMatrix(vector<vector<LogMzPeak>> &matrix, vector<LogMzPeak> &result){
       priority_queue< ppi, vector<ppi>, greater<ppi> > pq;

       for (Size i=0; i< matrix.size(); i++){
           pq.push({matrix[i][0].logMz, {i, 0}});
       }

       while (!pq.empty()) {
           ppi curr = pq.top();
           pq.pop();

           // i ==> Array Number (filter index)
           // j ==> Index in the array number (peak index)
           Size i = curr.second.first;
           Size j = curr.second.second;

           result.push_back(matrix[i][j]);

           // The next element belongs to same array as current.
           if (j + 1 < matrix[i].size())
               pq.push({ matrix[i][j + 1].logMz, { i, j + 1 } });
       }
   }
*/

