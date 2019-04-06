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
//#include <Eigen/Dense>
//#include <cmath>



using namespace OpenMS;
using namespace std;
//using namespace Eigen;

class FLASHDeconv :
        public TOPPBase {

public:
    FLASHDeconv() :
            TOPPBase("FLASHDeconv",
                     "Ultra-fast high-quality deconvolution enables online processing of top-down MS data",
                     false) {}

    struct Parameter {
        int minCharge;
        double minMass;
        double maxMass;
        double tolerance;
        String fileName;// up to here: ordinary user accessible parameters

        double intensityThreshold;// advanced parameters
        int chargeRange;
        int minContinuousChargePeakCount;

        //int minContinuousIsotopeCount;
        int maxIsotopeCount;
        int maxMassCount;
        //double isotopeCosineThreshold;
        double chargeDistributionScoreThreshold;
        double minRTspan;

        vector<int> hCharges{2, 3, 5, 7}; // automated or fixed parameters
        double binWidth;
        int numOverlappedScans = 20;
        int threads = 1;
    };

    struct PrecalcularedAveragine {
        vector<IsotopeDistribution> isotopes;
        vector<double> norms;

        double massInterval;
        double minMass;

        PrecalcularedAveragine(double m, double M, double delta, CoarseIsotopePatternGenerator *generator)
                : massInterval(delta), minMass(m) {
            int i = 0;
            while (true) {
                double a = i * massInterval;
                i++;
                if (a < m) continue;
                if (a > M) break;
                auto iso = generator->estimateFromPeptideWeight(a);
                iso.trimRight(0.01 * iso.getMostAbundant().getIntensity());
                isotopes.push_back(iso);
                double norm = .0;
                for (Size k = 0; k < iso.size(); k++) {
                    norm += iso[k].getIntensity() * iso[k].getIntensity();
                }
                norms.push_back(norm);
            }
        }

        IsotopeDistribution get(double mass) {
            Size i = (Size) (.5 + (mass - minMass) / massInterval);
            i = i >= isotopes.size() ? isotopes.size() - 1 : i;
            return isotopes[i];
        }

        double getNorm(double mass) {
            Size i = (Size) (.5 + (mass - minMass) / massInterval);
            i = i >= isotopes.size() ? isotopes.size() - 1 : i;
            return norms[i];
        }
    };

    struct LogMzPeak {
        Peak1D *orgPeak;
        double logMz;
        double mass = .0;
        int charge;
        int isotopeIndex = -1;
        //double score;

        LogMzPeak() : orgPeak(nullptr), logMz(-1000), charge(0), isotopeIndex(0) {}

        explicit LogMzPeak(Peak1D &peak) : orgPeak(&peak), logMz(getLogMz(peak.getMZ())), charge(0), isotopeIndex(0) {}

        LogMzPeak(Peak1D &peak, int c, int i) : orgPeak(&peak), logMz(getLogMz(peak.getMZ())), charge(c),
                                                isotopeIndex(i) {}

        ~LogMzPeak() {

        }

        double getMass() {
            if (mass <= 0) mass = exp(logMz) * charge;
            return mass;
        }

        bool operator<(const LogMzPeak &a) const {
            return this->logMz < a.logMz;
        }

        bool operator>(const LogMzPeak &a) const {
            return this->logMz >= a.logMz;
        }

    };

    struct PeakGroup {
        vector<LogMzPeak> peaks;
        double monoisotopicMass = .0;
        double avgMass = .0;
        double intensity = .0;
        double chargeDistributionScore = .0;

        double isotopeCosineScore = .0;
        int massIndex, specIndex, massCntr;
        MSSpectrum *spec;

        ~PeakGroup() {
            //vector<LogMzPeak>().swap(peaks);
        }

        void push_back(LogMzPeak &p) {
            peaks.push_back(p);
        }

        void reserve(Size n) {
            peaks.reserve(n);
        }

        bool operator<(const PeakGroup &a) const {
            return this->spec->getRT() < a.spec->getRT();
        }

        bool operator>(const PeakGroup &a) const {
            return this->spec->getRT() >= a.spec->getRT();
        }

    };

protected:

    static double getLogMz(double mz) {
        return log(mz - Constants::PROTON_MASS_U);
    }

    static int getNominalMass(double &m) {
        return (int) (m * 0.999497 + .5);
    }

    static double getBinValue(Size bin, double minV, double binWidth) {
        return minV + bin / binWidth;
    }

    static Size getBinNumber(double v, double minV, double binWidth) {
        if (v < minV) return 0;
        return (Size) (((v - minV) * binWidth) + .5);

    }


    void registerOptionsAndFlags_() override {
        registerInputFile_("in", "<input file>", "", "Input file");
        //setValidFormats_("in", ListUtils::create<String>("mzML"));
        registerOutputFile_("out", "<output file prefix/output dir>", "",
                            "Output file prefix or output dir (if prefix, [file prefix].tsv , [file prefix]feature.tsv, and [file prefix].m will be generated. "
                            "if dir, [dir]/[inputfile].tsv, [dir]/[inputfile]feature.tsv, and [dir]/[inputfile].m are generated per [inputfile])");

        registerIntOption_("minC", "<min charge>", 2, "minimum charge state", false, false);
        registerIntOption_("maxC", "<max charge>", 100, "maximum charge state", false, false);
        registerDoubleOption_("minM", "<min mass>", 4000.0, "minimum mass (Da)", false, false);
        registerIntOption_("minCC", "<min continuous charge peak count>", 3,
                           "minimum number of peaks of continuous charges per mass", false, true);
        //registerIntOption_("minCC", "<min charge count>", 4,
        //                  "minimum number of peaks of distinct charges per mass (recommended - ~25% of (maxC - minC))",
        //                 false, true);
        //registerIntOption_("minIC", "<min isotope count>", 3, "minimum continuous isotope count", false, true);
        registerIntOption_("maxIC", "<max isotope count>", 300, "maximum isotope count", false, true);
        registerIntOption_("maxMC", "<max mass count>", -1, "maximum mass count per spec", false, true);
        registerDoubleOption_("minCDScore", "<...>", .3, "minimum charge distribution score threshold",
                              false, true);

        registerDoubleOption_("maxM", "<max mass>", 150000.0, "maximum mass (Da)", false, false);
        registerDoubleOption_("tol", "<tolerance>", 5.0, "ppm tolerance", false, false);
        registerDoubleOption_("minInt", "<min intensity>", 0.0, "intensity threshold", false, true);
        //registerDoubleOption_("minIsoScore", "<score 0-1>", .6, "minimum isotope cosine score threshold (0-1)", false,
        //                     true);
        registerDoubleOption_("minRTspan", "<min RT span in seconds>", 20.0,
                              "minimum RT span for feature", false, true);

        //registerDoubleOption_("RTwindow", "<max RT delta>", 10.0, "max retention time duration with no peak in a feature (seconds); if negative, no feature finding performed", false, true);
    }

    Parameter setParameter() {
        Parameter param;
        param.minCharge = getIntOption_("minC");
        param.chargeRange = getIntOption_("maxC") - param.minCharge + 1;
        param.maxMass = getDoubleOption_("maxM");
        param.minMass = getDoubleOption_("minM");
        param.tolerance = getDoubleOption_("tol") * 1e-6;
        param.binWidth = .5 / param.tolerance;
        param.intensityThreshold = getDoubleOption_("minInt");
        param.minContinuousChargePeakCount = getIntOption_("minCC");
        //param.minChargeCount = getIntOption_("minCC");
        //param.minContinuousIsotopeCount = getIntOption_("minIC");
        param.maxIsotopeCount = getIntOption_("maxIC");
        param.maxMassCount = getIntOption_("maxMC");
        //param.isotopeCosineThreshold = getDoubleOption_("minIsoScore");
        param.chargeDistributionScoreThreshold = getDoubleOption_("minCDScore");
        param.minRTspan = getDoubleOption_("minRTspan");
        param.threads = getIntOption_("threads");
        //param.minChargeCoverage = .5;
        return param;
    }

    ExitCodes main_(int, const char **) override {
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

        if (check_file.isDir()) {
            QDirIterator it(path, QStringList() << "*.mzml", QDir::Files, QDirIterator::Subdirectories);
            while (it.hasNext()) {
                infileArray.push_back(it.next());
            }
        } else {
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

        if (!isOutPathDir) {
            fs.open(outfilePath + ".tsv", fstream::out);
            if (param.minRTspan > 0) fsf.open(outfilePath + "feature.tsv", fstream::out);

            writeHeader(fs, fsf, param.minRTspan > 0);
            fsm.open(outfilePath + ".m", fstream::out);
            fsm << "m=[";

            fsp.open(outfilePath + "peak.m", fstream::out);
        }

        for (auto &infile : infileArray) {
            if (isOutPathDir) {
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
            for (auto it = map.begin(); it != map.end(); ++it) {
                if (it->getMSLevel() != 1) continue;
                ms1Cntr++;
            }

            double rtDuration = (map[map.size() - 1].getRT() - map[0].getRT()) / ms1Cntr;
            param.numOverlappedScans = max(20, (int) (.5 + param.minRTspan * 2 / rtDuration)); // TODO if 0 ,segmentation fault!
            //cout<<param.numOverlappedScans<<endl;
            if (isOutPathDir) {
                std::string outfileName(param.fileName);
                std::size_t found = outfileName.find_last_of(".");
                outfileName = outfileName.substr(0, found);
                fs.open(outfilePath + outfileName + ".tsv", fstream::out);
                if (param.minRTspan > 0) fsf.open(outfilePath + "m" + outfileName + "feature.tsv", fstream::out);
                writeHeader(fs, fsf, param.minRTspan > 0);

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
                writePeakGroup(pg, param, fs, fsm, fsp);
            cout << "done" << endl;

            if (param.minRTspan > 0 && !peakGroups.empty() && specCntr > 0 && map.size() > 1) {
                findFeatures(peakGroups, map, featureCntr, fsf, param);
            }

            if (isOutPathDir) {
                cout << "In this run, FLASHDeconv found " << massCntr << " masses in " << qspecCntr
                     << " MS1 spectra out of "
                     << specCntr << endl;
                if (featureCntr > 0) cout << "Mass tracer found " << featureCntr << " features" << endl;

                fsm << "];";
                fsm.close();
                fsp.close();
                fs.close();
                if (param.minRTspan > 0) fsf.close();

                total_specCntr += specCntr;
                total_qspecCntr += qspecCntr;
                total_massCntr += massCntr;
                total_featureCntr += featureCntr;
            } else {
                cout << "So far, FLASHDeconv found " << massCntr << " masses in " << qspecCntr
                     << " MS1 spectra out of "
                     << specCntr << endl;
                if (featureCntr > 0) cout << "Mass tracer found " << featureCntr << " features" << endl;

                total_specCntr = specCntr;
                total_qspecCntr = qspecCntr;
                total_massCntr = massCntr;
                total_featureCntr = featureCntr;
            }
            total_elapsed_cpu_secs += elapsed_cpu_secs;
            total_elapsed_wall_secs += elapsed_wall_secs;

            //vector<PeakGroup>().swap(peakGroups);
        }

        cout << "-- done [took " << total_elapsed_cpu_secs << " s (CPU), " << total_elapsed_wall_secs
             << " s (Wall)] --"
             << endl;
        cout << "-- per spectrum [took " << 1000.0 * total_elapsed_cpu_secs / total_specCntr
             << " ms (CPU), " << 1000.0 * total_elapsed_wall_secs / total_specCntr << " ms (Wall)] --" << endl;

        if (massCntr < total_massCntr) {
            cout << "In total, FLASHDeconv found " << total_massCntr << " masses in " << total_qspecCntr
                 << " MS1 spectra out of "
                 << total_specCntr << endl;
            if (featureCntr > 0) cout << "Mass tracer found " << total_featureCntr << " features" << endl;
        }

        if (!isOutPathDir) {
            fsm << "];";
            fsm.close();
            fsp.close();
            fs.close();
            if (param.minRTspan > 0) fsf.close();
        }
        return EXECUTION_OK;
    }

    void
    findFeatures(vector<PeakGroup> &peakGroups, MSExperiment &map, int &featureCntr, fstream &fsf, Parameter &param) {

        for (auto it = map.begin(); it != map.end(); ++it) {
            it->clear(false);
        }

        for (auto &pg : peakGroups) {
            auto &spec = pg.spec;
            Peak1D tp(pg.monoisotopicMass, (float) pg.intensity);//
            spec->push_back(tp);
        }
        for (auto it = map.begin(); it != map.end(); ++it) {
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

        //cout<<(int) (param.RTwindow / rtDuration)<<endl;
        mtd_param.setValue("min_sample_rate", 0.01, "");
        mtd_param.setValue("trace_termination_outliers", param.numOverlappedScans, "");
        mtd_param.setValue("min_trace_length", param.minRTspan, "");
        //mtd_param.setValue("max_trace_length", 1000.0, "");
        mtdet.setParameters(mtd_param);

        vector<MassTrace> m_traces;
        mtdet.run(map, m_traces);  // m_traces : output of this function

        for (auto &mt : m_traces) {
            auto mass = mt.getCentroidMZ();

            fsf << ++featureCntr << "\t" << param.fileName << "\t" << mass << "\t"
                << getNominalMass(mass) << "\t"
                << mt.begin()->getRT() << "\t"
                << mt.rbegin()->getRT() << "\t"
                << mt.getTraceLength() << "\t"
                << mt[mt.findMaxByIntPeak()].getRT() << "\t"
                << mt.getMaxIntensity(false) << "\t"
                << mt.computePeakArea() << "\n";
        }
        //vector<MassTrace>().swap(m_traces);
        //fsf.flush();
    }

    void writeHeader(fstream &fs, fstream &fsf, bool featureOut = false) {
        fs
                << "MassIndex\tSpecIndex\tFileName\tSpecID\tMassCountInSpec\tExactMass\tAvgMass\tNominalMass(round(ExactMass*0.999497))\t"
                   "PeakChargeRange\tPeakMinCharge\tPeakMaxCharge\t"
                   "AggregatedIntensity\tRetentionTime\tPeakCount\tPeakMZs\tPeakCharges\tPeakMasses\tPeakIsotopeIndices\t"
                   "PeakIntensities\tChargeDistScore\tIsotopeCosineScore\n";
        if (!featureOut) return;
        fsf << "ID\tFileName\tExactMass\tNominalMass\tStartRetentionTime"
               "\tEndRetentionTime\tRetentionTimeDuration\tApexRetentionTime"
               "\tMaxIntensity\tQuantity"
            << endl;

        return;
    }

    PrecalcularedAveragine getPrecalculatedAveragines(Parameter &param) {
        auto generator = new CoarseIsotopePatternGenerator();
        auto maxIso = generator->estimateFromPeptideWeight(param.maxMass);
        maxIso.trimRight(0.01 * maxIso.getMostAbundant().getIntensity());
        param.maxIsotopeCount = min(param.maxIsotopeCount, (int) maxIso.size() - 1);
        generator->setMaxIsotope((Size) param.maxIsotopeCount);
        return PrecalcularedAveragine(100, param.maxMass, 50, generator);
    }

    vector<PeakGroup> Deconvolution(MSExperiment &map, Parameter &param,
                                    PrecalcularedAveragine &averagines, int &specCntr, int &qspecCntr, int &massCntr) {

        double *filter = new double[param.chargeRange];
        long **hBinOffsets = new long *[param.chargeRange];

        for (int i = 0; i < param.chargeRange; i++) {
            filter[i] = log(1.0 / (i + param.minCharge)); // should be descending, and negative!
            hBinOffsets[i] = new long[param.hCharges.size()];
            for (Size k = 0; k < param.hCharges.size(); k++) {
                auto &hc = param.hCharges[k];
                float n = (float) (hc>0?   (hc / 2) : (-hc/2) + 1);
                auto harmonicFilter = log(1.0 / (i + n / hc + param.minCharge));
                hBinOffsets[i][k] = (long) round((filter[i] - harmonicFilter) * param.binWidth);
                //harmonicFilter = log(1.0 / (i - n / hc + param.minCharge));
                //hBinOffsets[i][k+1] = (long) round((filter[i] - harmonicFilter) * param.binWidth);
                //harmonicFilter = log(1.0 / (i - n / hc + param.minCharge));
                //hBinOffsets[i][k+1] = (long) floor((filter[i] - harmonicFilter) * param.binWidth);
            }
        }

        float prevProgress = .0;
        vector<PeakGroup> allPeakGroups;

        //to overlap previous mass bins.
        vector<vector<Size>> prevMassBinVector;
        vector<double> prevMinBinLogMassVector;
        //vector <MSSpectrum> prevSpectra;

        for (auto it = map.begin(); it != map.end(); ++it) {
            if (it->getMSLevel() != 1) continue;

            float progress = (float) (it - map.begin()) / map.size();
            if (progress > prevProgress + .01) {
                printProgress(progress); //
                prevProgress = progress;
            }

            specCntr++;

            auto logMzPeaks = getLogMzPeaks(*it, param);
            auto peakGroups = getPeakGroupsFromSpectrum(logMzPeaks, filter, hBinOffsets,
                                                        prevMassBinVector, prevMinBinLogMassVector,
                                                        param);

            if (peakGroups.empty()) continue;

            auto filteredPeakGroups = scoreAndFilterPeakGroups(peakGroups, averagines, param); //

            if (filteredPeakGroups.empty()) continue;
            //vector<PeakGroup>().swap(peakGroups);
            //prevPgs = peakGroups;
            qspecCntr++;
            for (auto &pg : filteredPeakGroups) {
                massCntr++;
                pg.spec = &(*it);
                pg.massIndex = massCntr;
                pg.specIndex = qspecCntr;
                pg.massCntr = (int) filteredPeakGroups.size();
                allPeakGroups.push_back(pg);
                //if (param.minRTspan <= 0) continue;

                // Peak1D tp(pg.monoisotopicMass, (float) pg.intensity);//
                //it->push_back(tp);
            }


            //vector<PeakGroup>().swap(filteredPeakGroups);
            //if (param.minRTspan > 0) it->sortByPosition();
        }

        delete[] filter;
        delete[] hBinOffsets;
        printProgress(1); //
        allPeakGroups.shrink_to_fit();
        return allPeakGroups; //
    }

    void writePeakGroup(PeakGroup &pg, Parameter &param, fstream &fs, fstream &fsm, fstream &fsp) {
        if (pg.peaks.empty()) return;
        double &m = pg.monoisotopicMass;
        double &am = pg.avgMass;
        double &intensity = pg.intensity;
        int nm = getNominalMass(m);
        sort(pg.peaks.begin(), pg.peaks.end());
        int minCharge = 100;
        int maxCharge = 0;
        for (auto &p : pg.peaks) {
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
        for (auto &p : pg.peaks) {
            fs << p.orgPeak->getMZ() << ";";
        }
        fs << "\t";
        for (auto &p : pg.peaks) {
            fs << p.charge << ";";
        }
        fs << "\t";
        for (auto &p : pg.peaks) {
            fs << p.getMass() << ";";
        }
        fs << "\t";
        for (auto &p : pg.peaks) {
            fs << p.isotopeIndex << ";";
        }
        fs << "\t";
        fs << fixed << setprecision(1);
        for (auto &p : pg.peaks) {
            fs << p.orgPeak->getIntensity() << ";";
        }
        fs << "\t" << pg.chargeDistributionScore << "\t" << pg.isotopeCosineScore
           << "\n";

        //cout<<1<<endl;

        fsp << "pg" << (int) (pg.monoisotopicMass * 10)  <<  "rt"  <<   (int)( pg.spec->getRT())
            << "=[";

        for (auto &p : pg.peaks) {
            fsp << p.charge << "," << p.isotopeIndex << "," << p.orgPeak->getIntensity() << ";";
        }

        fsp << "];\n";
        //cout<<3<<endl;

        fsm << m << "," << nm << "," << intensity << "," << pg.spec->getRT() << "\n";
        //cout<<4<<endl;


    }

    void printProgress(float progress) {
        int barWidth = 70;
        cout << "[";
        int pos = (int) (barWidth * progress);
        for (int i = 0; i < barWidth; ++i) {
            if (i < pos) std::cout << "=";
            else if (i == pos) std::cout << ">";
            else std::cout << " ";
        }
        cout << "] " << int(progress * 100.0) << " %\r";
        cout.flush();
    }

    vector<LogMzPeak> getLogMzPeaks(MSSpectrum &spec, const Parameter &param) {
        vector<LogMzPeak> logMzPeaks;
        //logMzPeaks.reserve(spec.size());

        /*if (prevSpectra.size() >= (Size) param.numOverlappedScans) {
            prevSpectra.erase(prevSpectra.begin());
        }
        prevSpectra.push_back(spec);
*/
        //for(auto &pspec: prevSpectra){
        for (auto &peak: spec) {
            if (peak.getIntensity() <= param.intensityThreshold) continue;
            LogMzPeak logMzPeak(peak);
            logMzPeaks.push_back(logMzPeak);
        }
        // }

        //sort(logMzPeaks.begin(), logMzPeaks.end());

        logMzPeaks.shrink_to_fit();

        return logMzPeaks;
    }


    vector<PeakGroup>
    getPeakGroupsFromSpectrum(vector<LogMzPeak> &logMzPeaks, double *filter, long **hBinOffsets,
                              vector<vector<Size>> &prevMassBinVector,
                              vector<double> &prevMinBinLogMassVector,
                              const Parameter &param) {
        double massBinMaxValue = min(
                logMzPeaks[logMzPeaks.size() - 1].logMz -
                filter[param.chargeRange - param.minContinuousChargePeakCount - 1],
                log(param.maxMass));

        double massBinMinValue = logMzPeaks[0].logMz - filter[param.minContinuousChargePeakCount];
        double mzBinMinValue = logMzPeaks[0].logMz;
        double mzBinMaxValue = logMzPeaks[logMzPeaks.size() - 1].logMz;
        Size massBinNumber = getBinNumber(massBinMaxValue, massBinMinValue, param.binWidth) + 1;

        long *binOffsets = new long[param.chargeRange];

        for (int i = 0; i < param.chargeRange; i++) {
            binOffsets[i] = (long) round((mzBinMinValue - filter[i] - massBinMinValue) * param.binWidth);
        }

        Size mzBinNumber = getBinNumber(mzBinMaxValue, mzBinMinValue, param.binWidth) + 1;
        Byte *logIntensities = new Byte[mzBinNumber];
        fill_n(logIntensities, mzBinNumber, 0);

        auto mzBins = getMzBins(logMzPeaks, mzBinMinValue, mzBinNumber, param.binWidth, logIntensities);
        boost::dynamic_bitset<> massBins(massBinNumber);

        auto unionPrevMassBins = getUnionMassBin(massBins, massBinMinValue, prevMassBinVector, prevMinBinLogMassVector,
                                                 param);
        auto perMassChargeRanges = getMassBins(massBins, mzBins, massBinMinValue,
                                               binOffsets,
                                               hBinOffsets,
                                               unionPrevMassBins,
                                               logIntensities,
                                               param);

        auto unionMassBins = unionPrevMassBins | massBins;
        auto peakGroups = getPeakGroupsWithMassBins(unionMassBins, massBins, logMzPeaks, mzBinMinValue,
                                                    binOffsets, perMassChargeRanges,
                                                    param);

        if (prevMassBinVector.size() >0 && prevMassBinVector.size() >= (Size) param.numOverlappedScans) {
//            auto &p = prevMassBinVector[0];
            //vector<Size>().swap(p);
            prevMassBinVector.erase(prevMassBinVector.begin());
            prevMinBinLogMassVector.erase(prevMinBinLogMassVector.begin());
        }

        auto index = massBins.find_first();
        vector<Size> mb;
        mb.reserve(massBins.count());
        while (index != massBins.npos) {
            mb.push_back(index);
            index = massBins.find_next(index);
        }
        prevMassBinVector.push_back(mb);
        prevMinBinLogMassVector.push_back(massBinMinValue);

        //clear memory
        prevMassBinVector.shrink_to_fit();
        prevMinBinLogMassVector.shrink_to_fit();
        //boost::dynamic_bitset<>().swap(mzBins);
        //boost::dynamic_bitset<>().swap(unionPrevMassBins);
        //boost::dynamic_bitset<>().swap(unionMassBins);
        //boost::dynamic_bitset<>().swap(massBins);
        delete[] binOffsets;
        for (int i = 0; i < 2; i++) {
            delete[] perMassChargeRanges[i]; // delete array within matrix
        }// delete actual matrix
        delete[] perMassChargeRanges;
        delete[] logIntensities;
        return peakGroups;
    }

    boost::dynamic_bitset<> getUnionMassBin(boost::dynamic_bitset<> &massBins, double &massBinMinValue,
                                            vector<vector<Size>> &prevMassBinVector,
                                            vector<double> &prevMassBinMinValue,
                                            const Parameter &param) {
        boost::dynamic_bitset<> u(massBins.size());
        if (u.size() == 0) return u;
        for (Size i = 0; i < prevMassBinVector.size(); i++) {
            auto &pmb = prevMassBinVector[i];
            if (pmb.empty()) continue;
            long shift = (long) (round((massBinMinValue - prevMassBinMinValue[i]) * param.binWidth));

            for (auto &index : pmb) {
                long j = index - shift;
                if (j < 0) continue;
                if ((Size) j >= u.size()) break;
                u[j] = true;
            }
        }
        return u;
    }

    vector<PeakGroup> getPeakGroupsWithMassBins(boost::dynamic_bitset<> &unionedMassBins,//double &binScoreThreshold,
                                                boost::dynamic_bitset<> &massBins,
                                                vector<LogMzPeak> &logMzPeaks,
                                                double &mzBinMinValue,
                                                long *binOffsets,
                                                Byte **chargeRanges,
                                                const Parameter &param) {
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

        auto massBinIndex = unionedMassBins.find_first();
        while (massBinIndex != unionedMassBins.npos) {
            int isoOff = 0;
            PeakGroup pg;
            pg.reserve(chargeRange * 2);

            for (auto j = minChargeRanges[massBinIndex]; j <= maxChargeRanges[massBinIndex]; j++) {
                int charge = j + minCharge;
                long &binOffset = binOffsets[j];
                auto &cpi = currentPeakIndex[j];
                double maxIntensity = 0.0;
                int maxIntensityPeakIndex = -1;

                while (cpi < logMzPeakSize - 1) {
                    auto bi = getBinNumber(logMzPeaks[cpi].logMz, mzBinMinValue, binWidth) + binOffset;
                    if (bi == massBinIndex) {
                        auto intensity = logMzPeaks[cpi].orgPeak->getIntensity();
                        if (intensity > maxIntensity) {
                            maxIntensity = intensity;
                            maxIntensityPeakIndex = cpi;
                        }
                    } else if (bi > massBinIndex) {
                        break;
                    }
                    cpi++;
                }

                if (maxIntensityPeakIndex >= 0) {
                    //PeakGroup pgc;
                    //pgc.reserve(100);
                    double mz = logMzPeaks[maxIntensityPeakIndex].orgPeak->getMZ() - Constants::PROTON_MASS_U;
                    double &logMz = logMzPeaks[maxIntensityPeakIndex].logMz;
                    double isof = Constants::C13C12_MASSDIFF_U / charge / mz;
                    int maxI = 0;
                    for (int d = -1; d <= 1; d += 2) { // negative then positive direction.
                        int peakIndex = maxIntensityPeakIndex + (d < 0 ? d : 0);
                        int lastPeakIndex = -100;
                        for (int i = 0; i < maxIsotopeCount && (peakIndex >= 0 && peakIndex < logMzPeakSize); i++) {
                            maxI = max(maxI, i);
                            double centerLogMz = logMz + isof * i * d;
                            double centerLogMzMin = centerLogMz - tol;
                            double centerLogMzMax = centerLogMz + tol;
                            bool isotopePeakPresent = false;
                            if (lastPeakIndex >= 0) peakIndex = lastPeakIndex;//maxIntensityPeakIndex + (d < 0 ? d : 0);
                            for (; peakIndex >= 0 && peakIndex < logMzPeakSize; peakIndex += d) {
                                double &observedLogMz = logMzPeaks[peakIndex].logMz;
                                if (observedLogMz < centerLogMzMin) { if (d < 0) { break; } else { continue; }}
                                if (observedLogMz > centerLogMzMax) { if (d < 0) { continue; } else { break; }}
                                isotopePeakPresent = true;
                                if (peakIndex != lastPeakIndex) {
                                    LogMzPeak p(*logMzPeaks[peakIndex].orgPeak, charge, i * d);
                                    auto bin = getBinNumber(p.logMz, mzBinMinValue, binWidth) + binOffset;

                                    //if(massBinIndex == bin || unionedMassBins[bin]) {
                                    pg.push_back(p);
                                    lastPeakIndex = peakIndex;
                                    //}
                                    if (massBinIndex != bin) {
                                        unionedMassBins[bin] = false; //
                                         massBins[bin] = false; //
                                    }
                                }
                            }
                            if (!isotopePeakPresent) break;
                            if (d < 0) {
                                isoOff = -i > isoOff ? isoOff : -i;
                            }
                        }
                    }

                    for (auto &p : pg.peaks) {// assign the nearest isotope index..
                        if (p.charge != charge) continue;
                        double minMzDelta = 100000.0;
                        for (int d = -1; d <= 1; d += 2) { // negative then positive direction.
                            for (int i = 0; i <= maxI; i++) {
                                double centerLogMz = logMz + isof * i * d;
                                double delta = abs(centerLogMz - p.logMz);
                                if (delta > minMzDelta) continue;
                                minMzDelta = delta;
                                p.isotopeIndex = i * d;
                            }
                        }
                        //pg.push_back(p);
                    }
                }
            }
            //cout<<pg.peaks.size()<<endl;
            pg.peaks.shrink_to_fit();
            if (!pg.peaks.empty()) {
                for (LogMzPeak &p : pg.peaks) {
                    //cout<<p.orgPeak->getMZ()<<" " << p.charge<< " " << p.isotopeIndex << endl;
                    p.isotopeIndex -= isoOff;
                }
                peakGroups.push_back(pg);
            }
            massBinIndex = unionedMassBins.find_next(massBinIndex);
        }
        delete[] currentPeakIndex;
        return peakGroups;
    }

    boost::dynamic_bitset<>
    getMzBins(vector<LogMzPeak> &logMzPeaks, double &mzBinMinValue, Size &binNumber, double binWidth,
              Byte *logIntensities
    ) {
        boost::dynamic_bitset<> mzBins(binNumber);
        double* intensities = new double[binNumber];
        fill_n(intensities, binNumber, .0);

        for (auto &p : logMzPeaks) {
            Size bi = getBinNumber(p.logMz, mzBinMinValue, binWidth);
            if (bi >= binNumber) continue;
            mzBins.set(bi);
            intensities[bi] += p.orgPeak->getIntensity();

            auto delta = (p.logMz - getBinValue(bi, mzBinMinValue, binWidth));

            //cout<<delta<<endl;
            if(delta * binWidth > .5){
                //add bin + 1
                if(bi<binNumber-1){
                    mzBins.set(bi + 1);
                    intensities[bi + 1] += p.orgPeak->getIntensity();
                }
            }else if(delta * binWidth< -.5){
                if(bi>0){
                    mzBins.set(bi - 1);
                    intensities[bi - 1] += p.orgPeak->getIntensity();
                }
            }


            //if (bi > 0) {
            //    mzBins.set(bi - 1);
            //}
            //if(p.orgPeak->getMZ() < 1310 && p.orgPeak->getMZ() > 1305){
            //    cout << fixed<<setprecision(4)<< 52*(p.orgPeak->getMZ() - Constants::PROTON_MASS_U) << " " << exp(getBinValue(bi, mzBinMinValue, binWidth))<<endl;
            //}
            //
        }

        for(Size i=0; i<binNumber;i++){
            if(intensities[i]<=0) continue;
            logIntensities[i] = (Byte)round(log2(intensities[i]));
        }

        delete[] intensities;
        //cout<<endl;
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


    Byte **getMassBins(boost::dynamic_bitset<> &massBins, boost::dynamic_bitset<> &mzBins,
                       double &massBinMinValue,
            //double &mzBinMinValue,
                       long *binOffsets,
            //double *filter,
                       long **hBinOffsets,
                       boost::dynamic_bitset<> &unionPrevMassBins,
                       Byte *logIntensities,
                       const Parameter &param) {

        long binThresholdMinMass = (long) getBinNumber(log(param.minMass), massBinMinValue, param.binWidth);
        //long binNumber = (long) massBins.size();


        boost::dynamic_bitset<> isQualified(massBins.size());
        Byte *continuousChargePeakPairCount = new Byte[massBins.size()];
        fill_n(continuousChargePeakPairCount, massBins.size(), 0);
        //Byte *noneContinuousChargePeakPairCount = new Byte[massBins.size()];
        //fill_n(noneContinuousChargePeakPairCount, massBins.size(), 0);

        getInitialMassBins(massBins, mzBins, isQualified, unionPrevMassBins, continuousChargePeakPairCount,
                           hBinOffsets, binOffsets,
                           logIntensities,
                           param,
                //massBinMinValue,
                           binThresholdMinMass);
        //cout<<"*"<<endl;
        // printMasses(isQualified, massBinMinValue, continuousChargePeakPairCount, param);

        auto perMassChargeRanges = getFinalMassBins(massBins, mzBins, isQualified, unionPrevMassBins,
                                                    continuousChargePeakPairCount,
                                                    binOffsets,
                                                    param,
                                                    binThresholdMinMass);
        // printMasses(massBins, massBinMinValue, continuousChargePeakPairCount, param);

        //delete[] noneContinuousChargePeakPairCount;
        delete[] continuousChargePeakPairCount;
        //boost::dynamic_bitset<>().swap(isQualified);
        return perMassChargeRanges;
    }

    void printMasses(boost::dynamic_bitset<> &massBins, double &massBinMinValue, Byte *continuousChargePeakPairCount,
                     const Parameter &param) {
        auto index = massBins.find_first();
        while (index != massBins.npos) {
            auto m = exp(getBinValue(index, massBinMinValue, param.binWidth));
            if (m < 3497 && m > 3496.5) {

                cout << m <<
                     " " << (int) continuousChargePeakPairCount[index] <<
                     //" " << (int) noneContinuousChargePeakPairCount[index] <<
                     endl;
            }
            index = massBins.find_next(index);
        }
        //cout << endl;
    }


    void getInitialMassBins(boost::dynamic_bitset<> &massBins,
                            boost::dynamic_bitset<> &mzBins,
                            boost::dynamic_bitset<> &isQualified,
                            boost::dynamic_bitset<> &unionPrevMassBins,
                            Byte *continuousChargePeakPairCount,
                            long **hBinOffsets,
                            long *binOffsets,
                            Byte *logIntensities,
                            const Parameter &param,
                            long &binStart) {

        int chargeRange = param.chargeRange;
        int hChargeSize = (int) param.hCharges.size();
        int minContinuousChargePeakCount = param.minContinuousChargePeakCount;
        long mzBinSize = (long) mzBins.size();
        long binEnd = (long) massBins.size();

        Byte *prevCharges = new Byte[massBins.size()];
        fill_n(prevCharges, massBins.size(), (Byte) (chargeRange + 2));

        Byte *prevIntensities = new Byte[massBins.size()];
        fill_n(prevIntensities, massBins.size(), 0);

        auto mzBinIndex = mzBins.find_first();
        while (mzBinIndex != mzBins.npos) {
            auto &logIntensity = logIntensities[mzBinIndex];
            for (Byte j = 0; j < chargeRange; j++) {
                long massBinIndex = mzBinIndex + binOffsets[j];
                if (massBinIndex < binStart) continue;
                if (massBinIndex >= binEnd) break;

                auto cd = prevCharges[massBinIndex] - j;
                prevCharges[massBinIndex] = j;
                auto id = prevIntensities[massBinIndex] - logIntensity;
                prevIntensities[massBinIndex] = logIntensity;
                if (cd != 1 || abs(id) > 1) {
                    continue;
                }

                bool h = false;
                auto &hbOffsets = hBinOffsets[j];
                for (int k = 0; k < hChargeSize; k++) {
                    long hbi = mzBinIndex - hbOffsets[k];// + rand() % 10000 - 5000 ;
                    for (int i = -1; i <= 1; i++) { //
                        auto bin = hbi + i;
                        if (bin < 0 || bin > mzBinSize) continue;
                        if (mzBins[bin] && abs(logIntensity - logIntensities[bin]) <= 1) { //
                            h = true;
                            unionPrevMassBins[massBinIndex] = false;
                            break;
                        }
                    }
                    if (h) break;
                }
                if (h) continue;

                isQualified[massBinIndex] =
                        ++continuousChargePeakPairCount[massBinIndex] >= minContinuousChargePeakCount; //
            }
            mzBinIndex = mzBins.find_next(mzBinIndex);
        }
        delete[] prevCharges;
        delete[] prevIntensities;
    }

    Byte **getFinalMassBins(boost::dynamic_bitset<> &massBins, boost::dynamic_bitset<> &mzBins,
                            boost::dynamic_bitset<> &isQualified,
                            boost::dynamic_bitset<> &unionPrevMassBins,
                            Byte *continuousChargePeakPairCount,
                            long *binOffsets,
                            const Parameter &param,
                            long binStart) {

        int chargeRange = param.chargeRange;
        Byte *maxChargeRanges = new Byte[massBins.size()];
        fill_n(maxChargeRanges, massBins.size(), 0);

        Byte *minChargeRanges = new Byte[massBins.size()];
        fill_n(minChargeRanges, massBins.size(), 200);

        auto mzBinIndex = mzBins.find_first();
        long binEnd = (long) massBins.size();

        auto toSkip = (isQualified | unionPrevMassBins).flip();

        while (mzBinIndex != mzBins.npos) {
            long maxIndex = -1;
            Byte max = 0;
            Byte charge = 0;

            //vector<Byte> setCharges;
            //setCharges.reserve(chargeRange);
            for (Byte j = 0; j < chargeRange; j++) {
                long massBinIndex = mzBinIndex + binOffsets[j];
                if (massBinIndex < binStart) continue;
                if (massBinIndex >= binEnd) break;
                if (toSkip[massBinIndex]) continue;

                auto &t = continuousChargePeakPairCount[massBinIndex];// + noneContinuousChargePeakPairCount[massBinIndex];//

                //if (isQualified[massBinIndex]) massBins[massBinIndex] = true;
                //maxChargeRanges[massBinIndex] = std::max(maxChargeRanges[massBinIndex], j);
                //minChargeRanges[massBinIndex] = std::min(minChargeRanges[massBinIndex], j);//

                if (max < t) {
                    //if (max < t) {
                    //    setCharges.clear();
                    max = t;
                    //}
                    maxIndex = massBinIndex;
                    //setCharges.push_back(j);
                    charge = j;
                }
            }
            //if (!setCharges.empty()) {
            //  for(auto & j : setCharges){
            if (maxIndex >= 0) {
                if (isQualified[maxIndex]) massBins[maxIndex] = true;
                maxChargeRanges[maxIndex] = std::max(maxChargeRanges[maxIndex], charge);
                minChargeRanges[maxIndex] = std::min(minChargeRanges[maxIndex], charge);
            }
            //     }
            //}

            mzBinIndex = mzBins.find_next(mzBinIndex);
        }

        Byte **chargeRanges = new Byte *[2];
        chargeRanges[0] = minChargeRanges;
        chargeRanges[1] = maxChargeRanges;
        return chargeRanges;
    }


    vector<PeakGroup> scoreAndFilterPeakGroups(vector<PeakGroup> &peakGroups,
                                               PrecalcularedAveragine &averagines,
                                               const Parameter &param) {
        vector<double> intensities;
        vector<PeakGroup> filteredPeakGroups;
        filteredPeakGroups.reserve(peakGroups.size());
        intensities.reserve(peakGroups.size());
        auto perChargeIsotopeIntensity = new double*[param.chargeRange];//check
        auto perIsotopeMaxCharge = new int[param.maxIsotopeCount];
        auto perIsotopeMinCharge = new int[param.maxIsotopeCount];
        auto perIsotopeIntensity = new double[param.maxIsotopeCount];

        for(int i=0;i<param.chargeRange;i++){
            perChargeIsotopeIntensity[i] = new double[param.maxIsotopeCount];
            //fill_n(perChargeIsotopeIntensity[i], param.maxIsotopeCount, 0);

        }

        for (auto &pg : peakGroups) {
            updatePerChargeIsotopeIntensity(perChargeIsotopeIntensity,
                                            perIsotopeMinCharge,  perIsotopeMaxCharge,
                                            perIsotopeIntensity, pg, param);

            double monoIsotopeMass, avgMass;
            pg.isotopeCosineScore = getIsotopeCosineAndDetermineExactMass(monoIsotopeMass, avgMass, pg,
                                                                          perIsotopeIntensity,
                                                                          param.maxIsotopeCount,
                                                                          averagines);
            double isotopeCosineThreshold;
            if (monoIsotopeMass < 10000) isotopeCosineThreshold = .8;
            else if (monoIsotopeMass > 150000) isotopeCosineThreshold = .0;
            else {
                isotopeCosineThreshold = .8 / (150000 - 10000) * (150000 - monoIsotopeMass);
            }

            //isotopeCosineThreshold = 0;

            if (pg.peaks.empty() || pg.isotopeCosineScore <= isotopeCosineThreshold) {
                continue;
            }
            //isotopeCosineThreshold = .6;
            updatePerChargeIsotopeIntensity(perChargeIsotopeIntensity,
                    //perChargeMinIsotope, perChargeMaxIsotope,
                                            perIsotopeMinCharge,  perIsotopeMaxCharge,
                                            perIsotopeIntensity, pg, param);

            pg.chargeDistributionScore = getDistributionScore2(perChargeIsotopeIntensity,
                    param.chargeRange, param.maxIsotopeCount, param.minContinuousChargePeakCount);

            if (pg.chargeDistributionScore < param.chargeDistributionScoreThreshold) {
                continue;
            }

            pg.chargeDistributionScore = getDistributionScore(perIsotopeMinCharge, perIsotopeMaxCharge, param.maxIsotopeCount, 0);

            if (pg.chargeDistributionScore < 0) {
                continue;
            }

            pg.monoisotopicMass = monoIsotopeMass;
            pg.avgMass = avgMass;
            pg.intensity = accumulate(perIsotopeIntensity, perIsotopeIntensity + param.maxIsotopeCount, .0);
            filteredPeakGroups.push_back(pg);
            intensities.push_back(pg.intensity);
            //delete[] perChargeMaxIsotope;
            //delete[] perChargeMinIsotope;

        }

        for(int i=0;i<param.chargeRange;i++) {
            delete[] perChargeIsotopeIntensity[i];
        }
        delete[] perChargeIsotopeIntensity;
        delete[] perIsotopeIntensity;
        delete[] perIsotopeMinCharge;
        delete[] perIsotopeMaxCharge;

        if (filteredPeakGroups.empty()) return filteredPeakGroups;
        filterPeakGroupsByIntensity(peakGroups, intensities, param);
        vector<double>().swap(intensities);
        filteredPeakGroups.shrink_to_fit();
        return filteredPeakGroups;
    }

    void updatePerChargeIsotopeIntensity(double **perChargeIsotopeIntensity,
            //int *perChargeMinIsotope, int *perChargeMaxIsotope,
                                         int *perIsotopeMinCharge, int *perIsotopeMaxCharge,
                                         double *perIsotopeIntensity, PeakGroup &pg,
                                         const Parameter &param) {
        //fill_n(perChargeMinIsotope, param.chargeRange, param.maxIsotopeCount);
        //fill_n(perChargeMaxIsotope, param.chargeRange, -1);
        for(int i=0;i<param.chargeRange;i++){
            //perChargeIsotopeIntensity[i] = new double[param.maxIsotopeCount];
            fill_n(perChargeIsotopeIntensity[i], param.maxIsotopeCount, 0);

        }

        fill_n(perIsotopeMinCharge, param.maxIsotopeCount, param.chargeRange);
        fill_n(perIsotopeMaxCharge, param.maxIsotopeCount, -1);

        fill_n(perIsotopeIntensity, param.maxIsotopeCount, 0);

// bool *tmp = new bool[param.chargeRange * param.maxIsotopeCount];
        for (auto &p : pg.peaks) {
            if (p.isotopeIndex < 0 || p.isotopeIndex >= param.maxIsotopeCount) continue;
            int index = p.charge - param.minCharge;
            perIsotopeIntensity[p.isotopeIndex] += p.orgPeak->getIntensity();
            perChargeIsotopeIntensity[index][p.isotopeIndex] += p.orgPeak->getIntensity();
            //perChargeMinIsotope[index] = min(perChargeMinIsotope[index], p.isotopeIndex);
            //perChargeMaxIsotope[index] = max(perChargeMaxIsotope[index], p.isotopeIndex);
            perIsotopeMinCharge[p.isotopeIndex] = min(perIsotopeMinCharge[p.isotopeIndex], index);
            perIsotopeMaxCharge[p.isotopeIndex] = max(perIsotopeMaxCharge[p.isotopeIndex], index);
            //  tmp[p.isotopeIndex + param.maxIsotopeCount * index]=true;
        }

        /*boost::dynamic_bitset<> discardedCharges(param.chargeRange);
        for(int i=0;i<param.chargeRange;i++){
            if(perChargeMaxIsotope[i] < 0)continue;
            int isoSpan = perChargeMaxIsotope[i] - perChargeMinIsotope[i];
            if(isoSpan < param.minContinuousIsotopeCount){ // discard corresponding peaks
                discardedCharges.set(i);
            }
        }

        vector<LogMzPeak> tmpPeaks;
        tmpPeaks.swap(pg.peaks);
        pg.reserve(tmpPeaks.size());

        for (auto &p : tmpPeaks) {
            if (discardedCharges[p.charge - param.minCharge]) continue;
            pg.push_back(p);
        }*/



        /*for(int j=0;j<param.chargeRange;j++){
            for(int i=0;i<param.maxIsotopeCount;i++){
                   // if(tmp[i + j * param.maxIsotopeCount])
                   //     perChargePeakCount[j] ++;
            }
        }*/

        // delete[] tmp;
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

    double getIsotopeCosineAndDetermineExactMass(double &monoIsotopeMass, double &avgMass, PeakGroup &pg,
                                                 double *perIsotopeIntensities,
                                                 int perIsotopeIntensitiesSize,
                                                 PrecalcularedAveragine &averagines) {
        auto iso = averagines.get(pg.peaks[0].getMass());
        auto isoNorm = averagines.getNorm(pg.peaks[0].getMass());
        int isoSize = (int) iso.size();

        double isoDiff = Constants::C13C12_MASSDIFF_U;// OpenMS::Math::mean(diffs.begin(), diffs.end());

        int offset = 0;
        double maxCosine = -1;
        int maxIsotopeIndex = 0, minIsotopeIndex = -1;

        for (int i = 0; i < perIsotopeIntensitiesSize; i++) {
            if (perIsotopeIntensities[i] <= 0) continue;
            maxIsotopeIndex = i;
            if (minIsotopeIndex < 0) minIsotopeIndex = i;
            //maxIsotopeIndex = p.isotopeIndex < maxIsotopeIndex ? maxIsotopeIndex : p.isotopeIndex;
            //minIsotopeIndex = p.isotopeIndex < minIsotopeIndex ? p.isotopeIndex : minIsotopeIndex;
        }

        for (int f = -isoSize + minIsotopeIndex; f <= maxIsotopeIndex; f++) {
            auto cos = getCosine(perIsotopeIntensities, minIsotopeIndex, maxIsotopeIndex, iso, isoSize, isoNorm, f);

            if (maxCosine < cos) {
                maxCosine = cos;
                offset = f;
            }
        }

        vector<LogMzPeak> tmpPeaks;
        tmpPeaks.swap(pg.peaks);
        pg.reserve(tmpPeaks.size());

        for (auto &p : tmpPeaks) {
            p.isotopeIndex -= offset;
            if (p.isotopeIndex < 0 || p.isotopeIndex >= isoSize) continue;
            pg.push_back(p);
        }

        int mostAbundantIndex = 0;
        auto mostAbundant = iso.getMostAbundant();
        for (int i = 0; i < isoSize; i++) {
            if (iso[i] == mostAbundant) {
                mostAbundantIndex = i;
                break;
            }
        }

        double maxIntensityForMonoIsotopeMass = -1;
        for (auto &p : pg.peaks) {
            if (p.isotopeIndex > maxIsotopeIndex - offset) continue;
            double intensity = p.orgPeak->getIntensity();
            if (maxIntensityForMonoIsotopeMass > intensity) continue;
            maxIntensityForMonoIsotopeMass = intensity;
            monoIsotopeMass = p.getMass() - p.isotopeIndex * isoDiff;
            avgMass = p.getMass() + (mostAbundantIndex - p.isotopeIndex) * isoDiff;
        }

        //monoIsotopeMass += offset * isoDiff;
        return maxCosine;
    }

    /*

    double getChargeDistributionScore2(int *perChargePeakCount, const Parameter &param) {

        int sum = 0;
        int nonZeroStart = -1, nonZeroEnd = 0;

        for (int i = 0; i < param.chargeRange; i++) {
            if(perChargePeakCount[i]<=0) continue;
            sum += perChargePeakCount[i];

            if (nonZeroStart < 0) nonZeroStart = i;
            nonZeroEnd = i;
        }

        int prevCharge = nonZeroStart;
        //double n2 = .0;
        int n2=0;
        boost::dynamic_bitset<> setCharge(param.chargeRange);
        for (int k = nonZeroStart + 1; k <= nonZeroEnd; k++) {
            if(perChargePeakCount[k]  <= 0) continue;

            //int &peakCount = perChargePeakCount[k];
            //if (peakCount <= 0) continue;

            if (k - prevCharge == 1) {
                setCharge[k]=true;
                setCharge[prevCharge]=true;
                n2++;
            }
            prevCharge = k;
        }
        if(n2<param.minContinuousChargePeakCount) return -100;

        double n1 = .0;

        auto i = setCharge.find_first();
        while(i!=setCharge.npos){

            n1 += perChargePeakCount[i];

            i = setCharge.find_next(i);
        }

        return n1 / sum;
    }
*/


    double getDistributionScore(int *mins, int* maxs, int range, int threshold){
        int nonZeroStart = -1, nonZeroEnd = 0;

        for (int i = 0; i < range; i++) {
            if (maxs[i] >= 0) {
                if (nonZeroStart < 0) nonZeroStart = i;
                nonZeroEnd = i;
            }
        }
        int maxSpan = nonZeroEnd - nonZeroStart;
        if (maxSpan <= 0) return -100.0;

        int prevCharge = nonZeroStart;
        double n1 = .0;
        double n2 = .0;
        int *n3 = new int[range];
        fill_n(n3, range, 0);
        double spanThreshold = maxSpan / 4.0;//

        for (int k = nonZeroStart + 1; k <= nonZeroEnd; k++) {
            if (maxs[k] < 0) {
                continue;
            }

            int intersectSpan = min(maxs[prevCharge], maxs[k])
                                - max(mins[prevCharge], mins[k]);

            if (k - prevCharge == 1) {
                if (spanThreshold <= intersectSpan) { //
                    n1++;
                } else {
                    n2++;
                }
            } else{// if(spanThreshold <= intersectSpan){
                n3[k - prevCharge]++;

            }
            prevCharge = k;
        }

        int n4 = 0;
        for (int i = 0; i < range; i++) {
            n4 = max(n4, n3[i]);
        }
        delete[] n3;
        //n2 = max(n2, n3);
        if (n4>2 && n1 < n4) return -100.0;
        if (n1 < threshold) return -100.0;
        return n1 / (n1 + n2);
    }

    double getDistributionScore2(double **perChargeIsotopeIntensity,
             int range, int range2, int threshold){

        auto perChargeIntensity = new double[range];

        int nonZeroStart = -1, nonZeroEnd = 0;
        //int maxIndex = 0;
//        double maxIntensity = .0;


        for (int i = 0; i < range; i++) {
            double intensity = 0;
            for(int j=0;j< range2 ; j++){
                intensity += perChargeIsotopeIntensity[i][j];
            }
            perChargeIntensity[i] = intensity;
            if (intensity > 0) {
                if (nonZeroStart < 0) nonZeroStart = i;
                nonZeroEnd = i;
            }

        }

        int prevCharge = nonZeroStart;
        double n1 = .0;
        double n2 = .0;
        int *n3 = new int[range];
        fill_n(n3, range, 0);

        /*auto cosines = new double[range];
        fill_n(cosines, range, 0);
        //double maxCosine = -1.0;
        for (int k = nonZeroStart; k <= nonZeroEnd; k++) {
            if (perChargeIntensity[k] <= 0) {
                continue;
            }

            cosines[k] = getCosine(perChargeIsotopeIntensity[k], perChargeIsotopeIntensity[maxIndex], range2);

        }*/
        //if(maxCosine <=0) return -100.0;
       // int* maxRange = getRange(perChargeIsotopeIntensity[maxIndex]);

        for (int k = nonZeroStart + 1; k <= nonZeroEnd; k++) {
            if (perChargeIntensity[k] <= 0) {
                continue;
            }

            if (k - prevCharge == 1) {
                double cos = getCosine(perChargeIsotopeIntensity[k], perChargeIsotopeIntensity[prevCharge], range2);

                if (cos >  .7) { //
                    n1++;
                } else {
                    n2++;
                }
            } else{//} if (cos>.8){// if(spanThreshold <= intersectSpan){
                n3[k - prevCharge]++;

            }
            prevCharge = k;
        }
        //delete[] maxRange;
        delete[] perChargeIntensity;

        int n4 = 0;
        for (int i = 0; i < range; i++) {
            n4 = max(n4, n3[i]);
        }
        delete[] n3;
       // delete[] cosines;
        //n2 = max(n2, n3);
        if (n4>1 && n1 <= n4) return -100.0;
        if (n1 < threshold) return -100.0;
        return n1 / (n1 + n2);
    }


/*
    double getDistributionScore3(double **perChargeIsotopeIntensity,  PrecalcularedAveragine &averagines,
            double& mass, int range, int range2, double isotopeCosineThreshold, int threshold){
        auto iso = averagines.get(mass);
        auto isoNorm = averagines.getNorm(mass);
        int isoSize = (int) iso.size();

        auto perChargeIntensity = new double[range];

        int nonZeroStart = -1, nonZeroEnd = 0;

        for (int i = 0; i < range; i++) {
            double intensity = 0;
            for(int j=0;j< range2 ; j++){
                intensity += perChargeIsotopeIntensity[i][j];
            }
            perChargeIntensity[i] = intensity;
            if (intensity > 0) {
                if (nonZeroStart < 0) nonZeroStart = i;
                nonZeroEnd = i;
            }

        }

        int prevCharge = nonZeroStart;
        double n1 = .0;
        double n2 = .0;
        int *n3 = new int[range];
        fill_n(n3, range, 0);


        for (int k = nonZeroStart + 1; k <= nonZeroEnd; k++) {
            if (perChargeIntensity[k] <= 0) {
                continue;
            }

            if (k - prevCharge == 1) {
                double a = .0, b = .0;

                for(int i=0;i<range2;i++){
                    a += min(perChargeIsotopeIntensity[k][i], perChargeIsotopeIntensity[prevCharge][i]);
                    b += max(perChargeIsotopeIntensity[k][i], perChargeIsotopeIntensity[prevCharge][i]);
                }

                double overlab = a/b;

                if (overlab > isotopeCosineThreshold) { //
                    n1++;
                } else {
                    n2++;
                }
            } else{//} if (cos>.8){// if(spanThreshold <= intersectSpan){
                n3[k - prevCharge]++;

            }
            prevCharge = k;
        }

        delete[] perChargeIntensity;

        int n4 = 0;
        for (int i = 0; i < range; i++) {
            n4 = max(n4, n3[i]);
        }
        delete[] n3;

        if (n4>1 && n1 <= n4) return -100.0;
        if (n1 < threshold) return -100.0;
        return n1 / (n1 + n2);
    }
*/

    double
    getCosine(double *a, int &aStart, int &aEnd, IsotopeDistribution &b, int &bSize, double &bNorm, int offset = 0) {
        double n = .0, d1 = .0;

        for (int j = aStart; j < aEnd; j++) {
            d1 += a[j] * a[j];
            int i = j - offset;
            if (i < 0 || i >= bSize) continue;
            n += a[j] * b[i].getIntensity();
        }

        double d = (d1 * bNorm);
        if (d <= 0) return 0;
        return n / sqrt(d);
    }

    double
    getCosine(double *a, double *b, int &size) {
        double n = .0, d1 = .0, d2=.0;

        for (int j = 0; j < size; j++) {
            d1 += a[j] * a[j];
            d2 += b[j] * b[j];
            n += a[j] * b[j];
        }

        double d = (d1 * d2);
        if (d <= 0) return 0;
        return n / sqrt(d);
    }

    void
    filterPeakGroupsByIntensity(vector<PeakGroup> &peakGroups, vector<double> &intensities, const Parameter &param) {
        if (param.maxMassCount < 0 || intensities.size() <= (Size) param.maxMassCount) return;
        Size mc = (Size) param.maxMassCount;
        sort(intensities.begin(), intensities.end());
        auto threshold = intensities[intensities.size() - mc];
        for (auto pg = peakGroups.begin(); pg != peakGroups.end();) {
            if (peakGroups.size() <= mc) break;
            if (pg->intensity < threshold) {
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

int main(int argc, const char **argv) {
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
