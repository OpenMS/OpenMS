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
#include <OpenMS/MATH/STATISTICS/CumulativeBinomial.h>


using namespace OpenMS;
using namespace std;

class FLASHDeconv :
    public TOPPBase
{

public:
  FLASHDeconv() :
      TOPPBase("FLASHDeconv & FeatureFinderIntact",
               "Ultra-fast high-quality deconvolution enables online processing of top-down MS data",
               false)
  {
  }

  struct Parameter
  {
    int minCharge;
    double minMass;
    double maxMass;
    double tolerance;
    String fileName;// up to here: ordinary user accessible parameters

    double intensityThreshold;// advanced parameters
    int chargeRange;
    int minContinuousChargePeakCount;
    int maxIsotopeCount;
    int maxMassCount;
    //double chargeDistributionScoreThreshold;
    double RTwindow;

    vector<int> hCharges{2, 3, 5, 7}; // automated or fixed parameters
    double binWidth;
    int minNumOverLappedScans = 10;
    int numOverlappedScans = minNumOverLappedScans;
    int threads = 1;
  };

  struct PrecalcularedAveragine
  {
    vector<IsotopeDistribution> isotopes;
    vector<double> norms;

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
        for (Size k = 0; k < iso.size(); k++)
        {
          norm += iso[k].getIntensity() * iso[k].getIntensity();
        }
        norms.push_back(norm);
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
  };

  struct LogMzPeak
  {
    Peak1D *orgPeak;
    double logMz;
    double mass = .0;
    int charge;
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
    int massIndex, specIndex, massCntr;
    int maxCharge, minCharge;

    MSSpectrum *spec;

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

    void updateMassesAndIntensity(int mostAbundantIndex = -1)
    {
      double maxIntensityForMonoIsotopeMass = -1;
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
        if (mostAbundantIndex >= 0)
        {
          avgMass = p.getMass() + (mostAbundantIndex - p.isotopeIndex) * Constants::C13C12_MASSDIFF_U;
        }
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
                        "Output file prefix or output dir (if prefix, [file prefix].tsv , [file prefix]feature.tsv, and [file prefix].m will be generated. "
                        "if dir, [dir]/[inputfile].tsv, [dir]/[inputfile]feature.tsv, and [dir]/[inputfile].m are generated per [inputfile])");
    registerIntOption_("minC", "<min charge>", 2, "minimum charge state", false, false);
    registerIntOption_("maxC", "<max charge>", 150, "maximum charge state", false, false);
    registerDoubleOption_("minM", "<min mass>", 1000.0, "minimum mass (Da)", false, false);
    registerIntOption_("minCC", "<min continuous charge peak count>", 3,
                       "minimum number of peaks of continuous charges per mass", false, true);
    //registerIntOption_("maxIC", "<max isotope count>", 300, "maximum isotope count", false, true);
    registerIntOption_("maxMC", "<max mass count>", -1, "maximum mass count per spec", false, true);
    //registerDoubleOption_("minCDScore", "<...>", .0, "minimum charge distribution score threshold", false, true);
    registerDoubleOption_("maxM", "<max mass>", 200000.0, "maximum mass (Da)", false, false);
    registerDoubleOption_("tol", "<tolerance>", 10.0, "ppm tolerance", false, false);
    registerDoubleOption_("minInt", "<min intensity>", 0.0, "intensity threshold", false, true);
    registerDoubleOption_("RTwindow", "<seconds>", 120.0,
                          "RT window for feature extraction", false, true);
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
    param.intensityThreshold = getDoubleOption_("minInt");
    param.minContinuousChargePeakCount = getIntOption_("minCC");
    //param.maxIsotopeCount = getIntOption_("maxIC");
    param.maxMassCount = getIntOption_("maxMC");
    //param.chargeDistributionScoreThreshold = getDoubleOption_("minCDScore");
    param.RTwindow = getDoubleOption_("RTwindow");
    param.threads = getIntOption_("threads");
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
    fstream fs, fsf, fsm, fsp;

    if (!isOutPathDir)
    {
      fs.open(outfilePath + ".tsv", fstream::out);
      if (param.RTwindow > 0)
      {
        fsf.open(outfilePath + "feature.tsv", fstream::out);
      }

      writeHeader(fs, fsf, param.RTwindow > 0);
      fsm.open(outfilePath + ".m", fstream::out);
      fsm << "m=[";

      fsp.open(outfilePath + "peak.m", fstream::out);
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
      cout << "Processing : " << infile.toStdString() << endl;

      mzml.setLogType(log_type_);
      mzml.load(infile, map);

      param.fileName = QFileInfo(infile).fileName().toStdString();

      int ms1Cntr = 0;
      for (auto it = map.begin(); it != map.end(); ++it)
      {
        if (it->getMSLevel() != 1)
        {
          continue;
        }
        ms1Cntr++;
      }

      double rtDuration = (map[map.size() - 1].getRT() - map[0].getRT()) / ms1Cntr;
      param.numOverlappedScans = max(param.minNumOverLappedScans, (int) (.5 + param.RTwindow / rtDuration));
      cout << "# Overlapped MS1 scans:" << param.numOverlappedScans << endl;
      if (isOutPathDir)
      {
        std::string outfileName(param.fileName);
        std::size_t found = outfileName.find_last_of(".");
        outfileName = outfileName.substr(0, found);
        fs.open(outfilePath + outfileName + ".tsv", fstream::out);
        if (param.RTwindow > 0)
        {
          fsf.open(outfilePath + "m" + outfileName + "feature.tsv", fstream::out);
        }
        writeHeader(fs, fsf, param.RTwindow > 0);

        outfileName.erase(std::remove(outfileName.begin(), outfileName.end(), '_'), outfileName.end());
        outfileName.erase(std::remove(outfileName.begin(), outfileName.end(), '-'), outfileName.end());
        fsm.open(outfilePath + "m" + outfileName + ".m", fstream::out);
        fsm << "m=[";
        fsp.open(outfilePath + "m" + outfileName + "peak.m", fstream::out);
      }

      cout << "Running FLASHDeconv ... " << endl;
      auto begin = clock();
      auto t_start = chrono::high_resolution_clock::now();
      //continue;
      auto peakGroups = Deconvolution(map, param, averagines, specCntr, qspecCntr, massCntr);
      auto t_end = chrono::high_resolution_clock::now();
      auto end = clock();
      elapsed_cpu_secs = double(end - begin) / CLOCKS_PER_SEC;
      elapsed_wall_secs = chrono::duration<double>(t_end - t_start).count();

      cout << endl << "writing results ...";
      cout.flush();


      for (auto &pg : peakGroups)
      {
        //  writePeakGroup(pg, param, fs, fsm, fsp);
      }
      cout << "done" << endl;

      if (param.RTwindow > 0 && !peakGroups.empty() && specCntr > 0 && map.size() > 1)
      {
        findFeatures(peakGroups, map, featureCntr, fsf, param);
      }

      if (isOutPathDir)
      {
        cout << "In this run, FLASHDeconv found " << massCntr << " masses in " << qspecCntr
             << " MS1 spectra out of "
             << specCntr << endl;
        if (featureCntr > 0)
        {
          cout << "Mass tracer found " << featureCntr << " features" << endl;
        }

        fsm << "];";
        fsm.close();
        fsp.close();
        fs.close();
        if (param.RTwindow > 0)
        {
          fsf.close();
        }

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
      total_elapsed_cpu_secs += elapsed_cpu_secs;
      total_elapsed_wall_secs += elapsed_wall_secs;
    }

    cout << "-- done [took " << total_elapsed_cpu_secs << " s (CPU), " << total_elapsed_wall_secs
         << " s (Wall)] --"
         << endl;
    cout << "-- per spectrum [took " << 1000.0 * total_elapsed_cpu_secs / total_specCntr
         << " ms (CPU), " << 1000.0 * total_elapsed_wall_secs / total_specCntr << " ms (Wall)] --" << endl;

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
      fsm << "];";
      fsm.close();
      fsp.close();
      fs.close();
      if (param.RTwindow > 0)
      {
        fsf.close();
      }
    }
    return EXECUTION_OK;
  }

  void
  findFeatures(vector<PeakGroup> &peakGroups, MSExperiment &map, int &featureCntr, fstream &fsf, Parameter &param)
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
      auto &spec = pg.spec;

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
    mtd_param.setValue("mass_error_ppm", param.tolerance * 1e6, "");
    mtd_param.setValue("trace_termination_criterion", "outlier", "");

    mtd_param.setValue("reestimate_mt_sd", "true", "");
    mtd_param.setValue("quant_method", "area", "");
    mtd_param.setValue("noise_threshold_int", .0, "");

    //double rtDuration = (map[map.size() - 1].getRT() - map[0].getRT()) / ms1Cntr;
    mtd_param.setValue("min_sample_rate", 0.01, "");
    mtd_param.setValue("trace_termination_outliers", param.numOverlappedScans, "");
    mtd_param.setValue("min_trace_length", .001, "");
    //mtd_param.setValue("max_trace_length", 1000.0, "");
    mtdet.setParameters(mtd_param);

    vector<MassTrace> m_traces;
    mtdet.run(map, m_traces);  // m_traces : output of this function

    for (auto &mt : m_traces)
    {
      int minCharge = param.chargeRange + param.minCharge + 1;
      int maxCharge = 0;
      boost::dynamic_bitset<> charges(param.chargeRange + param.minCharge + 1);

      for (auto &p2 : mt)
      {
        int specIndex = rtSpecMap[(float) p2.getRT()];
        auto &pgMap = peakGroupMap[specIndex];
        auto &pg = pgMap[(float) p2.getMZ()];
        minCharge = min(minCharge, pg.minCharge);
        maxCharge = max(maxCharge, pg.maxCharge);
        for (auto &p : pg.peaks)
        {
          charges[p.charge] = true;
        }
      }


      auto mass = mt.getCentroidMZ();
      fsf << ++featureCntr << "\t" << param.fileName << "\t" << to_string(mass) << "\t"
          //fsf << ++featureCntr << "\t" << param.fileName << "\t" << mass << "\t"
          << getNominalMass(mass) << "\t"
          << mt.begin()->getRT() << "\t"
          << mt.rbegin()->getRT() << "\t"
          << mt.getTraceLength() << "\t"
          << mt[mt.findMaxByIntPeak()].getRT() << "\t"
          << mt.getMaxIntensity(false) << "\t"
          << mt.computePeakArea() << "\t"
          << minCharge << "\t"
          << maxCharge << "\t"
          << charges.count() << "\n";
    }
    delete[] peakGroupMap;
  }

  static void writeHeader(fstream &fs, fstream &fsf, bool featureOut = false)
  {
    fs
        << "MassIndex\tSpecIndex\tFileName\tSpecID\tMassCountInSpec\tExactMass\tAvgMass\tNominalMass(round(ExactMass*0.999497))\t"
           "PeakChargeRange\tPeakMinCharge\tPeakMaxCharge\t"
           "AggregatedIntensity\tRetentionTime\tPeakCount\tPeakMZs\tPeakCharges\tPeakMasses\tPeakIsotopeIndices\t"
           "PeakIntensities\tIsotopeCosineScore\n";
    if (!featureOut)
    {
      return;
    }
    fsf << "ID\tFileName\tExactMass\tNominalMass\tStartRetentionTime"
           "\tEndRetentionTime\tRetentionTimeDuration\tApexRetentionTime"
           "\tMaxIntensity\tQuantity\tMinCharge\tMaxCharge\tChargeCount\n";
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
    long **hBinOffsets = new long *[param.chargeRange];

    for (int i = 0; i < param.chargeRange; i++)
    {
      filter[i] = log(1.0 / (i + param.minCharge)); // should be descending, and negative!
      hBinOffsets[i] = new long[param.hCharges.size()];
      for (Size k = 0; k < param.hCharges.size(); k++)
      {
        auto &hc = param.hCharges[k];
        float n = (float) (hc > 0 ? (hc / 2) : (-hc / 2) + 1);
        auto harmonicFilter = log(1.0 / (i + n / hc + param.minCharge));
        hBinOffsets[i][k] = (long) round((filter[i] - harmonicFilter) * param.binWidth);
      }
    }

    float prevProgress = .0;
    vector<PeakGroup> allPeakGroups;

    //to overlap previous mass bins.
    vector<vector<Size>> prevMassBinVector;
    vector<double> prevMinBinLogMassVector;
    //vector <MSSpectrum> prevSpectra;

    for (auto it = map.begin(); it != map.end(); ++it)
    {
      if (it->getMSLevel() != 1)
      {
        continue;
      }

      float progress = (float) (it - map.begin()) / map.size();
      if (progress > prevProgress + .01)
      {
        printProgress(progress); //
        prevProgress = progress;
      }

      specCntr++;

      auto logMzPeaks = getLogMzPeaks(*it, param);
      auto peakGroups = getPeakGroupsFromSpectrum(logMzPeaks, filter, hBinOffsets,
                                                  prevMassBinVector, prevMinBinLogMassVector,
                                                  averagines,
                                                  param);


      if (peakGroups.empty())
      {
        continue;
      }
      //vector<PeakGroup>().swap(peakGroups);
      //prevPgs = peakGroups;
      qspecCntr++;
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
    delete[] hBinOffsets;
    printProgress(1); //
    allPeakGroups.shrink_to_fit();
    return allPeakGroups; //
  }

  static void writePeakGroup(PeakGroup &pg, Parameter &param, fstream &fs, fstream &fsm, fstream &fsp)
  {
    //return;//
    if (pg.peaks.empty())
    {
      return;
    }
    double &m = pg.monoisotopicMass;
    double &am = pg.avgMass;
    double &intensity = pg.intensity;
    int nm = getNominalMass(m);
    sort(pg.peaks.begin(), pg.peaks.end());
    int minCharge = 100;
    int maxCharge = 0;
    for (auto &p : pg.peaks)
    {
      minCharge = minCharge < p.charge ? minCharge : p.charge;
      maxCharge = maxCharge > p.charge ? maxCharge : p.charge;
    }
    //cout<<1<<endl;
    fs << pg.massIndex << "\t" << pg.specIndex << "\t" << param.fileName << "\t" << pg.spec->getNativeID() << "\t"
       << pg.massCntr << "\t"
       << fixed << setprecision(3) << m << "\t" << am << "\t" << nm << "\t" <<
       (maxCharge - minCharge + 1) << "\t" << minCharge << "\t" << maxCharge << "\t"
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
    fs << "\t" << pg.isotopeCosineScore
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
    fsm << m << "," << nm << "," << intensity << "," << pg.spec->getRT() << "\n";
    //cout<<4<<endl;


  }

  static void printProgress(float progress)
  {
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
    for (auto &peak: spec)
    {
      if (peak.getIntensity() <= param.intensityThreshold)
      {
        continue;
      }
      LogMzPeak logMzPeak(peak);
      logMzPeaks.push_back(logMzPeak);
    }

    logMzPeaks.shrink_to_fit();

    return logMzPeaks;
  }


  static vector<PeakGroup>
  getPeakGroupsFromSpectrum(vector<LogMzPeak> &logMzPeaks, double *filter, long **hBinOffsets,
                            vector<vector<Size>> &prevMassBinVector,
                            vector<double> &prevMinBinLogMassVector,
                            PrecalcularedAveragine &averagines,
                            const Parameter &param)
  {
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
      binOffsets[i] = (long) round((mzBinMinValue - filter[i] - massBinMinValue) * param.binWidth);
    }

    Size mzBinNumber = getBinNumber(mzBinMaxValue, mzBinMinValue, param.binWidth) + 1;
    float *logIntensities = new float[mzBinNumber];
    fill_n(logIntensities, mzBinNumber, 0);

    auto mzBins = getMzBins(logMzPeaks, mzBinMinValue, mzBinNumber, param.binWidth, logIntensities);
    boost::dynamic_bitset<> massBins(massBinNumber);

    auto unionMassBins = getUnionMassBin(massBins, massBinMinValue, prevMassBinVector, prevMinBinLogMassVector,
                                         param);
    auto perMassChargeRanges = getMassBins(massBins, mzBins, massBinMinValue,
                                           binOffsets,
                                           hBinOffsets,
                                           unionMassBins,
                                           logIntensities,
                                           param);

    auto peakGroups = getPeakGroupsWithMassBins(unionMassBins, logMzPeaks, mzBinMinValue,
                                                binOffsets, perMassChargeRanges,
                                                param);

    auto filteredPeakGroups = scoreAndFilterPeakGroups(peakGroups, averagines, param);
    if (prevMassBinVector.size() > 0 && prevMassBinVector.size() >= (Size) param.numOverlappedScans)
    {
      prevMassBinVector.erase(prevMassBinVector.begin());
      prevMinBinLogMassVector.erase(prevMinBinLogMassVector.begin());
    }

    vector<Size> mb;
    mb.reserve(filteredPeakGroups.size());
    for (auto &pg : filteredPeakGroups)
    {
      if (massBins[pg.massBinIndex])
      {
        mb.push_back(pg.massBinIndex);
      }
    }

    prevMassBinVector.push_back(mb);
    prevMinBinLogMassVector.push_back(massBinMinValue);

    prevMassBinVector.shrink_to_fit();
    prevMinBinLogMassVector.shrink_to_fit();

    delete[] binOffsets;
    for (int i = 0; i < 2; i++)
    {
      delete[] perMassChargeRanges[i]; // delete array within matrix
    }// delete actual matrix
    delete[] perMassChargeRanges;
    delete[] logIntensities;

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
    int *currentPeakIndex = new int[param.chargeRange];
    fill_n(currentPeakIndex, param.chargeRange, 0); //

    vector<PeakGroup> peakGroups;
    auto &minChargeRanges = chargeRanges[0];
    auto &maxChargeRanges = chargeRanges[1];
    vector<Size> toRemove;

    auto massBinIndex = unionedMassBins.find_first();
    // Size lastSetMassBinIndex = unionedMassBins.size();
    Size* peakBinNumbers = new Size[logMzPeakSize];
    for(int i=0;i<logMzPeakSize;i++){
      peakBinNumbers[i] =  getBinNumber(logMzPeaks[i].logMz, mzBinMinValue, binWidth);
    }

    while (massBinIndex != unionedMassBins.npos)
    {
      int isoOff = 0;
      PeakGroup pg;
      pg.reserve(chargeRange * 2);

      for (auto j = minChargeRanges[massBinIndex]; j <= maxChargeRanges[massBinIndex]; j++)
      {
        int charge = j + minCharge;
        long &binOffset = binOffsets[j];
        auto &cpi = currentPeakIndex[j];
        double maxIntensity = 0.0;
        int maxIntensityPeakIndex = -1;

        while (cpi < logMzPeakSize - 1)
        {
          auto bi = peakBinNumbers[cpi] + binOffset;
          if (bi == massBinIndex)
          {
            auto intensity = logMzPeaks[cpi].orgPeak->getIntensity();
            if (intensity > maxIntensity)
            {
              maxIntensity = intensity;
              maxIntensityPeakIndex = cpi;
            }
          }
          else if (bi > massBinIndex)
          {
            break;
          }
          cpi++;
        }

        if (maxIntensityPeakIndex >= 0)
        {
          double mz = logMzPeaks[maxIntensityPeakIndex].orgPeak->getMZ() - Constants::PROTON_MASS_U;
          double &logMz = logMzPeaks[maxIntensityPeakIndex].logMz;
          double isof = Constants::C13C12_MASSDIFF_U / charge / mz;
          int maxI = 0;
          for (int d = -1; d <= 1; d += 2)
          { // negative then positive direction.
            int peakIndex = maxIntensityPeakIndex + (d < 0 ? d : 0);
            int lastPeakIndex = -100;
            for (int i = 0; i < maxIsotopeCount && (peakIndex >= 0 && peakIndex < logMzPeakSize); i++)
            {
              maxI = max(maxI, i);
              double centerLogMz = logMz + isof * i * d;
              double centerLogMzMin = centerLogMz - tol;
              double centerLogMzMax = centerLogMz + tol;
              bool isotopePeakPresent = false;
              if (lastPeakIndex >= 0)
              {
                peakIndex = lastPeakIndex;
              }//maxIntensityPeakIndex + (d < 0 ? d : 0);
              for (; peakIndex >= 0 && peakIndex < logMzPeakSize; peakIndex += d)
              {
                double &observedLogMz = logMzPeaks[peakIndex].logMz;
                if (observedLogMz < centerLogMzMin)
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
                if (observedLogMz > centerLogMzMax)
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
                  LogMzPeak p(*logMzPeaks[peakIndex].orgPeak, charge, i * d);
                  auto bin = peakBinNumbers[peakIndex] + binOffset;

                  //if(massBinIndex == bin || unionedMassBins[bin]) {
                  pg.push_back(p);
                  lastPeakIndex = peakIndex;
                  //}
                  if (massBinIndex != bin && unionedMassBins[bin])
                  {
                    toRemove.push_back(bin);
                  }
                }
              }
              if (!isotopePeakPresent)
              {
                break;
              }
              //if (d < 0) {
              //    isoOff = -i > isoOff ? isoOff : -i;
              //}
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
                double centerLogMz = logMz + isof * (p.isotopeIndex + i * d);
                double delta = abs(centerLogMz - p.logMz);

                if (delta > minMzDelta)
                {
                  break;
                }
                maxId = i * d;
                minMzDelta = delta;
              }
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
        for (auto &p : pg.peaks)
        {
          p.isotopeIndex -= isoOff;
        }

        pg.updateMassesAndIntensity();
        pg.massBinIndex = massBinIndex;
        /*if(lastSetMassBinIndex == massBinIndex-1){// remove duplicate
            auto prevIntensity = peakGroups[peakGroups.size()-1].intensity;
            auto currentIntensity = pg.intensity;
            if(prevIntensity < currentIntensity){
                peakGroups[peakGroups.size()-1] = pg;
            }
        }else*/
        peakGroups.push_back(pg);
        //lastSetMassBinIndex = massBinIndex;
      }

      if (massBinIndex < unionedMassBins.size() - 1 && !unionedMassBins[massBinIndex + 1])
      {
        for (auto &b : toRemove)
        {
          unionedMassBins[b] = false;
        }
        toRemove.clear();
      }

      massBinIndex = unionedMassBins.find_next(massBinIndex);
    }
    delete[] currentPeakIndex;
    delete[] peakBinNumbers;
    return peakGroups;
  }

  static boost::dynamic_bitset<>
  getMzBins(vector<LogMzPeak> &logMzPeaks, double &mzBinMinValue, Size &binNumber, double binWidth,
            float *logIntensities
  )
  {
    boost::dynamic_bitset<> mzBins(binNumber);
    double *intensities = new double[binNumber];
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

      if (delta * binWidth > 0)
      {
        //add bin + 1
        if (bi < binNumber - 1)
        {
          mzBins.set(bi + 1);
          intensities[bi + 1] += p.orgPeak->getIntensity();
        }
      }
      else if (delta * binWidth < 0)
      {
        if (bi > 0)
        {
          mzBins.set(bi - 1);
          intensities[bi - 1] += p.orgPeak->getIntensity();
        }
      }
    }

    for (Size i = 0; i < binNumber; i++)
    {
      if (intensities[i] <= 0)
      {
        continue;
      }
      logIntensities[i] = //(Byte) round
          (float) (log(intensities[i])/log(4.0));
    }

    delete[] intensities;
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
                            long *binOffsets,
                            long **hBinOffsets,
                            boost::dynamic_bitset<> &unionMassBins,
                            float *logIntensities,
                            const Parameter &param)
  {

    long binThresholdMinMass = (long) getBinNumber(log(param.minMass), massBinMinValue, param.binWidth);

    boost::dynamic_bitset<> isQualified(massBins.size());
    Byte *continuousChargePeakPairCount = new Byte[massBins.size()];
    fill_n(continuousChargePeakPairCount, massBins.size(), 0);

    float *sumLogIntensities = new float[massBins.size()];
    fill_n(sumLogIntensities, massBins.size(), 0);

    getInitialMassBins(massBins, mzBins, isQualified,
                       sumLogIntensities,
                       continuousChargePeakPairCount,
                       hBinOffsets, binOffsets,
                       logIntensities,
                       param,
                       binThresholdMinMass);

    //printMasses(isQualified, massBinMinValue, continuousChargePeakPairCount, param);
    //auto toSkip = (isQualified | isQualified).flip();
    auto perMassChargeRanges = getFinalMassBins(massBins, mzBins, isQualified, unionMassBins,
                                                sumLogIntensities,
        // continuousChargePeakPairCount,
                                                binOffsets,
                                                param,
                                                binThresholdMinMass);

    //printMasses(massBins, massBinMinValue, continuousChargePeakPairCount, param);


    delete[] continuousChargePeakPairCount;
    delete[] sumLogIntensities;
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
                                 float *sumLogIntensities,
                                 Byte *continuousChargePeakPairCount,
                                 long **hBinOffsets,
                                 long *binOffsets,
                                 float *logIntensities,
                                 const Parameter &param,
                                 long &binStart)
  {

    int chargeRange = param.chargeRange;
    int hChargeSize = (int) param.hCharges.size();
    int minContinuousChargePeakCount = param.minContinuousChargePeakCount;
    long mzBinSize = (long) mzBins.size();
    long binEnd = (long) massBins.size();

    Byte *prevCharges = new Byte[massBins.size()];
    fill_n(prevCharges, massBins.size(), (Byte) (chargeRange + 2));

    float *prevIntensities = new float[massBins.size()];
    fill_n(prevIntensities, massBins.size(), 0);

    //float idThreshold = 1.0;
    auto mzBinIndex = mzBins.find_first();
    while (mzBinIndex != mzBins.npos)
    {
      auto &logIntensity = logIntensities[mzBinIndex];
      for (Byte j = 0; j < chargeRange; j++)
      {
        long massBinIndex = mzBinIndex + binOffsets[j];
        if (massBinIndex < binStart)
        {
          continue;
        }
        if (massBinIndex >= binEnd)
        {
          break;
        }

        //int idt = j < ct ? 2 : 1;
        auto cd = prevCharges[massBinIndex] - j;
        prevCharges[massBinIndex] = j;
        auto &prevIntensity = prevIntensities[massBinIndex];
        auto id = prevIntensity - logIntensity;
        //prevIntensities[massBinIndex] = logIntensity;
        if (cd != 1 || id > 1 || id < -1)
        {
          prevIntensity = logIntensity;
          continue;
        }

        //if (!isQualified[massBinIndex])

        bool h = false;
        auto &hbOffsets = hBinOffsets[j];
        for (int k = 0; k < hChargeSize; k++)
        {
          long hbi = mzBinIndex - hbOffsets[k];// + rand() % 10000 - 5000 ;
          for (int i = -2; i <= 2; i++)
          {
            auto bin = hbi + i;
            if (bin < 0 || bin > mzBinSize)
            {
              continue;
            }
            if (mzBins[bin]
                && (abs(logIntensity - logIntensities[bin]) <= 1))// ||
                   // abs(prevIntensity - logIntensities[bin]) <= 1))
            {
              h = true;
              break;
            }
          }
          if (h)
          {
            break;
          }
        }
        prevIntensity = logIntensity;
        if (h)
        {
          continue;
        }
        sumLogIntensities[massBinIndex] += logIntensity;
        if (isQualified[massBinIndex])
        {
          continue;
        }
        isQualified[massBinIndex] =
            ++continuousChargePeakPairCount[massBinIndex] >= minContinuousChargePeakCount;

      }
      mzBinIndex = mzBins.find_next(mzBinIndex);
    }
    delete[] prevCharges;
    delete[] prevIntensities;
    //delete[] continuousChargePeakPairCount;
  }

  static Byte **getFinalMassBins(boost::dynamic_bitset<> &massBins, boost::dynamic_bitset<> &mzBins,
                                 boost::dynamic_bitset<> &isQualified,
                                 boost::dynamic_bitset<> &unionMassBins,
                                 float *sumLogIntensities,
      // Byte *continuousChargePeakPairCount,
                                 long *binOffsets,
                                 const Parameter &param,
                                 long binStart)
  {
    int chargeRange = param.chargeRange;
    Byte *maxChargeRanges = new Byte[massBins.size()];
    fill_n(maxChargeRanges, massBins.size(), 0);

    Byte *minChargeRanges = new Byte[massBins.size()];
    fill_n(minChargeRanges, massBins.size(), chargeRange + 1);

    //Byte *selected = new Byte[massBins.size()];
    //fill_n(selected, massBins.size(), 0);

    auto mzBinIndex = mzBins.find_first();
    long binEnd = (long) massBins.size();
    //int minContinuousChargePeakCount = param.minContinuousChargePeakCount;

    auto toSkip = (isQualified | unionMassBins).flip();
    unionMassBins.reset();

    while (mzBinIndex != mzBins.npos)
    {
      long maxIndex = -1;
      int maxCount = 0;
      Byte charge = 0;

      //vector<Byte> setCharges;
      //setCharges.reserve(chargeRange);
      for (Byte j = 0; j < chargeRange; j++)
      {
        long massBinIndex = mzBinIndex + binOffsets[j];
        if (massBinIndex < binStart)
        {
          continue;
        }
        if (massBinIndex >= binEnd)
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

        if (maxCount < t)
        {
          maxCount = t;
          maxIndex = massBinIndex;
          charge = j;
        }
      }

      if (maxIndex >= 0)
      {
        maxChargeRanges[maxIndex] = max(maxChargeRanges[maxIndex], charge);
        minChargeRanges[maxIndex] = min(minChargeRanges[maxIndex], charge);
        //++selected[maxIndex];
        //int t = max(continuousChargePeakPairCount[maxIndex],minContinuousChargePeakCount);
        //if (selected[maxIndex] >= 2)
          //&& selected[maxIndex] >= continuousChargePeakPairCount[maxIndex])
        {
          massBins[maxIndex] = isQualified[maxIndex];
          unionMassBins[maxIndex] = true;
        }
      }
      mzBinIndex = mzBins.find_next(mzBinIndex);
    }

    Byte **chargeRanges = new Byte *[2];
    chargeRanges[0] = minChargeRanges;
    chargeRanges[1] = maxChargeRanges;
    //delete[] selected;
    return chargeRanges;
  }

  static double
  getIsotopeCosineThreshold(double &mass, double minCosine, double maxCosine, double minMass, double maxMass)
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

    return isotopeCosineThreshold;
  }

  static vector<PeakGroup> scoreAndFilterPeakGroups(vector<PeakGroup> &peakGroups,
                                                    PrecalcularedAveragine &averagines,
                                                    const Parameter &param)
  {

    //removeOverlappingPeakGroups(peakGroups, param.tolerance);
    vector<double> intensities;
    vector<PeakGroup> filteredPeakGroups;
    filteredPeakGroups.reserve(peakGroups.size());
    intensities.reserve(peakGroups.size());
    auto perChargeIsotopeIntensity = new double *[param.chargeRange];//check
    auto perIsotopeMaxCharge = new int[param.maxIsotopeCount];
    auto perIsotopeMinCharge = new int[param.maxIsotopeCount];
    auto perChargeMaxIsotope = new int[param.chargeRange];
    auto perChargeMinIsotope = new int[param.chargeRange];// polish later

    auto perIsotopeIntensity = new double[param.maxIsotopeCount];

    for (int i = 0; i < param.chargeRange; i++)
    {
      perChargeIsotopeIntensity[i] = new double[param.maxIsotopeCount];
    }

    for (auto &pg : peakGroups)
    {
      updatePerChargeIsotopeIntensity(perChargeIsotopeIntensity,
                                      perIsotopeMinCharge, perIsotopeMaxCharge,
                                      perChargeMinIsotope, perChargeMaxIsotope,
                                      perIsotopeIntensity, pg, param);

      bool isIsotopeSpanGood = checkSpanDistribution(perChargeMinIsotope, perChargeMaxIsotope,
                                                     param.chargeRange, param.minContinuousChargePeakCount);

      if (!isIsotopeSpanGood)
      {
        continue;
      }

      bool isChargeSpanGood = checkSpanDistribution(perIsotopeMinCharge, perIsotopeMaxCharge,
                                                    param.maxIsotopeCount, 1);

      if (!isChargeSpanGood)
      {
        continue;
      }

      bool isChargeWellDistributed = checkChargeDistribution(perChargeIsotopeIntensity,
                                                             param.chargeRange, param.maxIsotopeCount,
                                                             param.minContinuousChargePeakCount);

      if (!isChargeWellDistributed)
      {
        continue;
      }

      pg.isotopeCosineScore = getIsotopeCosineAndDetermineExactMass(pg,
                                                                    perIsotopeIntensity,
                                                                    param.maxIsotopeCount,
                                                                    averagines);

      //double isotopeCosineThreshold = getIsotopeCosineThreshold(pg.monoisotopicMass, .8, .8, 10000, 100000);

      if (pg.peaks.empty() || pg.isotopeCosineScore <= .6)
      {
        continue;
      }

      filteredPeakGroups.push_back(pg);
      intensities.push_back(pg.intensity);
    }

    removeOverlappingPeakGroups(filteredPeakGroups, param.tolerance);

    for (int i = 0; i < param.chargeRange; i++)
    {
      delete[] perChargeIsotopeIntensity[i];
    }
    delete[] perChargeIsotopeIntensity;
    delete[] perIsotopeIntensity;
    delete[] perIsotopeMinCharge;
    delete[] perIsotopeMaxCharge;
    delete[] perChargeMinIsotope;
    delete[] perChargeMaxIsotope;

    if (filteredPeakGroups.empty())
    {
      return filteredPeakGroups;
    }
    filterPeakGroupsByIntensity(peakGroups, intensities, param);
    vector<double>().swap(intensities);
    filteredPeakGroups.shrink_to_fit();
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
      for (Size j = i + 1; j < pgs.size(); j++)
      {
        auto &pgo = pgs[j];
        if (pgo.monoisotopicMass - pg.monoisotopicMass > pg.monoisotopicMass * tol * 2)
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
        if (pg.monoisotopicMass - pgo.monoisotopicMass > pg.monoisotopicMass * tol * 2)
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
  }

  static void updatePerChargeIsotopeIntensity(double **perChargeIsotopeIntensity,
      //int *perChargeMinIsotope, int *perChargeMaxIsotope,
                                              int *perIsotopeMinCharge, int *perIsotopeMaxCharge,
                                              int *perChargeMinIsotope, int *perChargeMaxIsotope,
                                              double *perIsotopeIntensity, PeakGroup &pg,
                                              const Parameter &param)
  {
    for (int i = 0; i < param.chargeRange; i++)
    {
      //perChargeIsotopeIntensity[i] = new double[param.maxIsotopeCount];
      fill_n(perChargeIsotopeIntensity[i], param.maxIsotopeCount, 0);

    }

    fill_n(perIsotopeMinCharge, param.maxIsotopeCount, param.chargeRange);
    fill_n(perIsotopeMaxCharge, param.maxIsotopeCount, -1);

    fill_n(perChargeMinIsotope, param.chargeRange, param.maxIsotopeCount);
    fill_n(perChargeMaxIsotope, param.chargeRange, -1);

    fill_n(perIsotopeIntensity, param.maxIsotopeCount, 0);

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
      perChargeIsotopeIntensity[index][p.isotopeIndex] += p.orgPeak->getIntensity();
      perChargeMinIsotope[index] = min(perChargeMinIsotope[index], p.isotopeIndex);
      perChargeMaxIsotope[index] = max(perChargeMaxIsotope[index], p.isotopeIndex);
      perIsotopeMinCharge[p.isotopeIndex] = min(perIsotopeMinCharge[p.isotopeIndex], index);
      perIsotopeMaxCharge[p.isotopeIndex] = max(perIsotopeMaxCharge[p.isotopeIndex], index);
      //  tmp[p.isotopeIndex + param.maxIsotopeCount * index]=true;
    }
    pg.maxCharge = maxCharge;
    pg.minCharge = minCharge;

  }


  /*
  bool isIsotopeIntensityQualified(double *perIsotopeIntensity, int *perChargeMaxIsotope, const Parameter &param) {
      int maxSetIntensityCounter = 0;
      int setIntensityCounter = 0;


      for (int i = 0; i < param.chargeRange; i++) {
          if (perChargeMaxIsotope[i] < 0) {
              setIntensityCounter = 0;
              continue;
          }
          setIntensityCounter++;
          maxSetIntensityCounter =
                  maxSetIntensityCounter > setIntensityCounter ? maxSetIntensityCounter : setIntensityCounter;
      }

      if(maxSetIntensityCounter < param.minContinuousChargePeakCount) return false;

      maxSetIntensityCounter = 0;
      setIntensityCounter = 0;

      for (int i = 0; i < param.maxIsotopeCount; i++) {
          if (perIsotopeIntensity[i] <= 0) {
              setIntensityCounter = 0;
              continue;
          }
          setIntensityCounter++;
          maxSetIntensityCounter =
                  maxSetIntensityCounter > setIntensityCounter ? maxSetIntensityCounter : setIntensityCounter;
      }

      return maxSetIntensityCounter >= param.minContinuousIsotopeCount;
  }*/

  static double getIsotopeCosineAndDetermineExactMass(PeakGroup &pg,
                                                      double *perIsotopeIntensities,
                                                      int perIsotopeIntensitiesSize,
                                                      PrecalcularedAveragine &averagines)
  {
    auto iso = averagines.get(pg.peaks[0].getMass());
    auto isoNorm = averagines.getNorm(pg.peaks[0].getMass());
    int isoSize = (int) iso.size();

    //double isoDiff = Constants::C13C12_MASSDIFF_U;// OpenMS::Math::mean(diffs.begin(), diffs.end());

    int offset = 0;
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

      if (maxCosine < cos)
      {
        maxCosine = cos;
        offset = f;
      }
    }

    vector<LogMzPeak> tmpPeaks;
    tmpPeaks.swap(pg.peaks);
    pg.reserve(tmpPeaks.size());

    for (auto &p : tmpPeaks)
    {
      p.isotopeIndex -= offset;
      if (p.isotopeIndex < 0 || p.isotopeIndex >= isoSize)
      {
        continue;
      }
      pg.push_back(p);
    }

    int mostAbundantIndex = 0;
    double mostAbundantInt = 0;
    //auto mostAbundant = iso.getMostAbundant();
    for (int i = 0; i < isoSize; i++)
    {
      if (mostAbundantInt >= iso[i].getIntensity())
      {
        continue;
      }
      mostAbundantInt = iso[i].getIntensity();
      mostAbundantIndex = i;

    }

    pg.updateMassesAndIntensity(mostAbundantIndex);
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
    int n_h = 0;

    double spanThreshold = maxSpan / 4.0;//

    for (int k = nonZeroStart + 1; k <= nonZeroEnd; k++)
    {
      if (maxs[k] < 0)
      {
        continue;
      }

      int intersectSpan = min(maxs[prevCharge], maxs[k])
                          - max(mins[prevCharge], mins[k]);

      if (k - prevCharge == 1)
      {
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

          int intersectSpan = min(maxs[prevCharge], maxs[k])
                              - max(mins[prevCharge], mins[k]);

          if (k - prevCharge == i)
          {
            if (spanThreshold <= intersectSpan)
            {
              t++;
            }
          }
          prevCharge = k;
        }
        n_h = max(n_h, t);
      }
    }

    return n_r > n_h;
  }

  static bool checkChargeDistribution(double **perChargeIsotopeIntensity,
                                      int range, int range2,
      //double cost,
                                      int threshold)
  {
    auto perChargeIntensity = new double[range];
    double maxPerChargeIntensity = .0;
    int nonZeroStart = -1, nonZeroEnd = 0;
    for (int i = 0; i < range; i++)
    {
      double intensity = 0;
      for (int j = 0; j < range2; j++)
      {
        intensity += perChargeIsotopeIntensity[i][j];
      }
      perChargeIntensity[i] = intensity;
      if (intensity > 0)
      {
        // intensities.push_back(intensity);
        maxPerChargeIntensity = max(maxPerChargeIntensity, intensity);
        if (nonZeroStart < 0)
        {
          nonZeroStart = i;
        }
        nonZeroEnd = i;
      }
    }

    int prevCharge = nonZeroStart;

    int n_r = .0;
    int n_h = 0;

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
      delete[] perChargeIntensity;
      return false;
    }

    for (int i = 2; i < min(12, range); i++)
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
        n_h = max(n_h, t);
      }
    }

    delete[] perChargeIntensity;
    return n_r > n_h;
  }

  static double
  getCosine(double *a, int &aStart, int &aEnd, IsotopeDistribution &b, int &bSize, double &bNorm, int offset = 0)
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
      n += a[j] * b[i].getIntensity();
    }

    double d = (d1 * bNorm);
    if (d <= 0)
    {
      return 0;
    }
    return n / sqrt(d);
  }

  static double
  getCosine(double *a, double *b, int &size)
  {
    double n = .0, d1 = .0, d2 = .0;

    //int overlapCntr = 0;
    for (int j = 0; j < size; j++)
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
