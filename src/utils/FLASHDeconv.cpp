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
        double tolerance; // up to here: ordinary user accessible parameters
        double binWidth;
        String fileName;

        double intensityThreshold;// advanced parameters
        int chargeRange;
        int minContinuousChargePeakPairCount;

        int minContinuousIsotopeCount;
        int maxIsotopeCount;
        int maxMassCount;
        double isotopeCosineThreshold;
        int chargeDistributionScoreThreshold;
        double maxRTDelta;
        vector<Byte> hCharges{2, 3, 5, 7};
        int numOverlappedScans = 10;
        int threads = 1;
    };

    struct PrecalcularedAveragine {
        vector<IsotopeDistribution> isotopes;
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
                iso.trimRight(.01 * iso.getMostAbundant().getIntensity());
                isotopes.push_back(iso);
            }
        }

        IsotopeDistribution get(double mass) {
            Size i = (Size) ((mass - minMass) / massInterval);
            i = i >= isotopes.size() ? isotopes.size() - 1 : i;
            return isotopes[i];
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

        double getMass() {
            if (mass <= 0) mass = exp(logMz) * charge;
            return mass;
        }

        //double getMonoIsotopeMass() {
        //    return getMass() - isotopeIndex * Constants::C13C12_MASSDIFF_U;
        //}

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
        double intensity = .0;
        int chargeDistributionScore = 0;
        double isotopeCosineScore = .0;
        int massIndex, specIndex, massCntr;
        MSSpectrum *spec;

        void push_back(LogMzPeak &p) {
            peaks.push_back(p);
        }

        void reserve(Size n) {
            peaks.reserve(n);
        }

        bool operator<(const PeakGroup &a) const {
            // if(this->spec->getRT() == a.spec->getRT()) return this->monoisotopicMass < a.monoisotopicMass;
            return this->spec->getRT() < a.spec->getRT();
        }

        bool operator>(const PeakGroup &a) const {
            // if(this->spec->getRT() == a.spec->getRT()) return this->monoisotopicMass > a.monoisotopicMass;
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
        registerIntOption_("maxC", "<max charge>", 60, "maximum charge state", false, false);
        registerDoubleOption_("minM", "<min mass>", 1000.0, "minimum mass (Da)", false, false);
        registerIntOption_("minCC", "<min continuous charge peak pair count>", 3,
                           "minimum number of peaks of continuous charges per mass", false, true);
        //registerIntOption_("minCC", "<min charge count>", 4,
        //                  "minimum number of peaks of distinct charges per mass (recommended - ~25% of (maxC - minC))",
        //                 false, true);
        registerIntOption_("minIC", "<min isotope count>", 3, "minimum continuous isotope count", false, true);
        registerIntOption_("maxIC", "<max isotope count>", 100, "maximum isotope count", false, true);
        registerIntOption_("maxMC", "<max mass count>", -1, "maximum mass count per spec", false, true);
        registerIntOption_("minCDScore", "<score 0,1,2,...>", 1, "minimum charge distribution score threshold (>= 0)",
                           false, true);

        registerDoubleOption_("maxM", "<max mass>", 100000.0, "maximum mass (Da)", false, false);
        registerDoubleOption_("tol", "<tolerance>", 5.0, "ppm tolerance", false, false);
        registerDoubleOption_("minInt", "<min intensity>", 0.0, "intensity threshold", false, true);
        registerDoubleOption_("minIsoScore", "<score 0-1>", .5, "minimum isotope cosine score threshold (0-1)", false,
                              true);
        registerDoubleOption_("maxRTDelta", "<maximum RT between masses for feature finding>", 20.0,
                              "maximum RT between masses for feature finding", false, true);

        //registerDoubleOption_("maxRTDelta", "<max RT delta>", 10.0, "max retention time duration with no peak in a feature (seconds); if negative, no feature finding performed", false, true);
    }

    Parameter setParameter() {
        Parameter param;
        param.minCharge = getIntOption_("minC");
        param.chargeRange = getIntOption_("maxC") - param.minCharge + 1;
        param.maxMass = getDoubleOption_("maxM");
        param.minMass = getDoubleOption_("minM");
        param.tolerance = getDoubleOption_("tol") * 1e-6;
        param.binWidth = 2.0 / param.tolerance;
        param.intensityThreshold = getDoubleOption_("minInt");
        param.minContinuousChargePeakPairCount = getIntOption_("minCC");
        //param.minChargeCount = getIntOption_("minCC");
        param.minContinuousIsotopeCount = getIntOption_("minIC");
        param.maxIsotopeCount = getIntOption_("maxIC");
        param.maxMassCount = getIntOption_("maxMC");
        param.isotopeCosineThreshold = getDoubleOption_("minIsoScore");
        param.chargeDistributionScoreThreshold = getIntOption_("minCDScore");
        param.maxRTDelta = getDoubleOption_("maxRTDelta");
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
            QDirIterator it(path, QStringList() << "*.mzml", QDir::Files);
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
        fstream fs, fsm, fsf;
        MSExperiment map;
        MzMLFile mzml;

        if (!isOutPathDir) {
            fs.open(outfilePath + ".tsv", fstream::out);
            if (param.maxRTDelta > 0) fsf.open(outfilePath + "feature.tsv", fstream::out);

            writeHeader(fs, fsf, param.maxRTDelta > 0);
            fsm.open(outfilePath + ".m", fstream::out);
            fsm << "m=[";
        }

        for (auto &infile : infileArray) {
            if (isOutPathDir) {
                specCntr = qspecCntr = massCntr = featureCntr = 0;
            }

            double elapsed_cpu_secs = 0, elapsed_wall_secs = 0;
            cout << "Processing : " << infile.toStdString() << endl;

            mzml.setLogType(log_type_);
            mzml.load(infile, map);

            param.fileName = QFileInfo(infile).fileName().toStdString();

            if (isOutPathDir) {
                std::string outfileName(param.fileName);
                std::size_t found = outfileName.find_last_of(".");
                outfileName = outfileName.substr(0, found);
                fs.open(outfilePath + outfileName + ".tsv", fstream::out);
                if (param.maxRTDelta > 0) fsf.open(outfilePath + "m" + outfileName + "feature.tsv", fstream::out);
                writeHeader(fs, fsf, param.maxRTDelta > 0);

                outfileName.erase(std::remove(outfileName.begin(), outfileName.end(), '_'), outfileName.end());
                outfileName.erase(std::remove(outfileName.begin(), outfileName.end(), '-'), outfileName.end());
                fsm.open(outfilePath + "m" + outfileName + ".m", fstream::out);
                fsm << "m=[";
            }

            //param.numOverlappedScans = param.numOverlappedScans > param.maxNumOverlappedScans ? param.maxNumOverlappedScans : param.numOverlappedScans;
            //cout<<param.numOverlappedScans<<endl;
            cout << "Running FLASHDeconv ... " << endl;
            auto begin = clock();
            auto t_start = chrono::high_resolution_clock::now();
            auto peakGroups = Deconvolution(map, param, averagines, specCntr, qspecCntr, massCntr);
            auto t_end = chrono::high_resolution_clock::now();
            auto end = clock();
            elapsed_cpu_secs = double(end - begin) / CLOCKS_PER_SEC;
            elapsed_wall_secs = chrono::duration<double>(t_end - t_start).count();


            // if (!peakGroups.empty()) {
            cout << endl << "writing results ...";
            cout.flush();


            for (auto &pg : peakGroups)
                writePeakGroup(pg, param, fs, fsm);
            cout << "done\n";

            if (param.maxRTDelta > 0 && !peakGroups.empty() && specCntr > 0 && map.size() > 1) {
                findFeatures(map, featureCntr, fsf, specCntr, param);
            }


            cout << "In this run, FLASHDeconv found " << massCntr << " masses in " << qspecCntr
                 << " MS1 spectra out of "
                 << specCntr << endl;
            if (featureCntr > 0) cout << "Mass tracer found " << featureCntr << " features" << endl;

            if (isOutPathDir) {
                fsm << "];";
                fsm.close();

                fs.close();
                if (param.maxRTDelta > 0) fsf.close();

                total_specCntr += specCntr;
                total_qspecCntr += qspecCntr;
                total_massCntr += massCntr;
                total_featureCntr += featureCntr;
            } else {
                total_specCntr = specCntr;
                total_qspecCntr = qspecCntr;
                total_massCntr = massCntr;
                total_featureCntr = featureCntr;
            }
            total_elapsed_cpu_secs += elapsed_cpu_secs;
            total_elapsed_wall_secs += elapsed_wall_secs;
            peakGroups.clear();
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
            fs.close();
            if (param.maxRTDelta > 0) fsf.close();
        }
        return EXECUTION_OK;
    }

    void findFeatures(MSExperiment &map, int &featureCntr, fstream &fsf, int &ms1Cntr, Parameter &param) {
        Param common_param = getParam_().copy("algorithm:common:", true);
        writeDebug_("Common parameters passed to sub-algorithms (mtd and ffm)", common_param, 3);

        Param mtd_param = getParam_().copy("algorithm:mtd:", true);
        writeDebug_("Parameters passed to MassTraceDetection", mtd_param, 3);

        MassTraceDetection mtdet;
        mtd_param.insert("", common_param);
        mtd_param.remove("chrom_fwhm");

        //mtd_param.setValue("mass_error_da", .2,// * (param.chargeRange+ param.minCharge),
        //                   "Allowed mass deviation (in da).");
        mtd_param.setValue("mass_error_ppm", param.tolerance * 1e6 * 10, "");
        mtd_param.setValue("trace_termination_criterion", "outlier", "");
        //mtd_param.setValue("trace_termination_criterion", "sample_rate", "");

        mtd_param.setValue("reestimate_mt_sd", "true", "");
        mtd_param.setValue("quant_method", "area", "");
        mtd_param.setValue("noise_threshold_int", .0, "");

        double rtDuration = (map[map.size() - 1].getRT() - map[0].getRT()) / ms1Cntr;

        mtd_param.setValue("min_sample_rate", 0.01, "");
        mtd_param.setValue("trace_termination_outliers", (int) (param.maxRTDelta / rtDuration), "");
        mtd_param.setValue("min_trace_length", param.maxRTDelta / 2, "");
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
        fsf.flush();
    }

    void writeHeader(fstream &fs, fstream &fsf, bool featureOut = false) {
        fs
                << "MassIndex\tSpecIndex\tFileName\tSpecID\tMassCountInSpec\tExactMass\tNominalMass(round(ExactMass*0.999497))\t"
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
        maxIso.trimRight(.01 * maxIso.getMostAbundant().getIntensity());
        param.maxIsotopeCount = min(param.maxIsotopeCount, (int) maxIso.size() - 1);
        generator->setMaxIsotope((Size) param.maxIsotopeCount);
        return PrecalcularedAveragine(param.minMass, param.maxMass, max(10.0, (param.maxMass - param.minMass) / 100.0),
                                      generator);
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
                float n = (float) (hc / 2);
                auto harmonicFilter = log(1.0 / (i - n / hc + param.minCharge));
                hBinOffsets[i][k] = (long) floor((filter[i] - harmonicFilter) * param.binWidth);
                //harmonicFilter = log(1.0 / (i - n / hc + param.minCharge));
                //hBinOffsets[i][k+1] = (long) floor((filter[i] - harmonicFilter) * param.binWidth);
            }
        }

        float prevProgress = .0;
        vector<PeakGroup> allPeakGroups;

        //to overlap previous mass bins.
        vector<vector<Size>> prevMassBinVector;
        vector<double> prevMinBinLogMassVector;

        for (auto it = map.begin(); it != map.end(); ++it) {
            if (it->getMSLevel() != 1) continue;

//            if(it->getRT() < 3000.0 || it->getRT() > 4000.0)
            //              continue;

            float progress = (float) (it - map.begin()) / map.size();
            if (progress > prevProgress + .05) {
                printProgress(progress);
                prevProgress = progress;
            }

            specCntr++;

            auto logMzPeaks = getLogMzPeaks(*it, param.intensityThreshold);
            auto peakGroups = getPeakGroupsFromSpectrum(logMzPeaks, filter, hBinOffsets,
                                                        prevMassBinVector, prevMinBinLogMassVector,
                    //fsmm,
                    //(it->getRT() > 360 && it->getRT() < 380),
                                                        param);
            if (param.maxRTDelta > 0) it->clear(false);

            if (peakGroups.empty()) continue;

            auto filteredPeakGroups = scoreAndFilterPeakGroups(peakGroups, averagines, param);

            if (filteredPeakGroups.empty()) continue;

            //prevPgs = peakGroups;
            qspecCntr++;
            for (auto &pg : filteredPeakGroups) {
                massCntr++;
                pg.spec = &(*it);
                pg.massIndex = massCntr;
                pg.specIndex = qspecCntr;
                pg.massCntr = (int) filteredPeakGroups.size();
                allPeakGroups.push_back(pg);
                if (param.maxRTDelta <= 0) continue;

                /*for(auto &p : pg.peaks){
                    Peak1D tp(p.getMonoIsotopeMass(), (float) p.orgPeak->getIntensity());//
                    it->push_back(tp);
                }*/


                Peak1D tp(pg.monoisotopicMass, (float) pg.intensity);//
                it->push_back(tp);
            }
            peakGroups.clear();
            filteredPeakGroups.clear();
            if (param.maxRTDelta > 0) it->sortByPosition();
        }

        printProgress(1);
        return allPeakGroups;
    }

    void writePeakGroup(PeakGroup &pg, Parameter &param, fstream &fs, fstream &fsm) {
        double m = pg.monoisotopicMass;
        double intensity = pg.intensity;
        int nm = getNominalMass(m);
        sort(pg.peaks.begin(), pg.peaks.end());
        int minCharge = 100;
        int maxCharge = 0;
        for (auto &p : pg.peaks) {
            minCharge = minCharge < p.charge ? minCharge : p.charge;
            maxCharge = maxCharge > p.charge ? maxCharge : p.charge;
        }


        fs << pg.massIndex << "\t" << pg.specIndex << "\t" << param.fileName << "\t" << pg.spec->getNativeID() << "\t"
           << pg.massCntr << "\t"
           << fixed << setprecision(3) << m << "\t" << nm << "\t" <<
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
        fs << "\t" << pg.chargeDistributionScore << "\t" << pg.isotopeCosineScore << "\n";

        fsm << m << "," << nm << "," << intensity << "," << pg.spec->getRT() << "\n";
        fs.flush();
        fsm.flush();
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

    vector<LogMzPeak> getLogMzPeaks(MSSpectrum &spec, double &threshold) {
        vector<LogMzPeak> logMzPeaks;
        logMzPeaks.reserve(spec.size());
/*
        for (int cpi = 1; cpi<spec.size()-1;cpi++ ) {
            auto &peak = spec[cpi];
            if (peak.getMZ() - spec[cpi - 1].getMZ() < peak.getMZ() * 1e-5) {
                if (peak.getIntensity() < spec[cpi - 1].getIntensity()) {
                    continue;
                }
            }

            if (spec[cpi + 1].getMZ() - peak.getMZ() < peak.getMZ() * 1e-5) {
                if (peak.getIntensity() < spec[cpi + 1].getIntensity()) {
                    continue;
                }
            }
            */
        for (auto &peak: spec) {
            if (peak.getIntensity() <= threshold) continue;
            LogMzPeak logMzPeak(peak);
            logMzPeaks.push_back(logMzPeak);
        }

        return logMzPeaks;
    }


    vector<PeakGroup>
    getPeakGroupsFromSpectrum(vector<LogMzPeak> &logMzPeaks, double *filter, long **hBinOffsets,
                              vector<vector<Size>> &prevMassBinVector,
                              vector<double> &prevMinBinLogMassVector,
            //fstream &fsmm,
                              const Parameter &param) {
        double massBinMaxValue = min(
                logMzPeaks[logMzPeaks.size() - 1].logMz -
                filter[param.chargeRange - param.minContinuousChargePeakPairCount - 1],
                log(param.maxMass));

        double massBinMinValue = logMzPeaks[0].logMz - filter[param.minContinuousChargePeakPairCount];
        double mzBinMinValue = logMzPeaks[0].logMz;
        double mzBinMaxValue = logMzPeaks[logMzPeaks.size() - 1].logMz;
        Size massBinNumber = getBinNumber(massBinMaxValue, massBinMinValue, param.binWidth) + 1;


        long *binOffsets = new long[param.chargeRange];

        for (int i = 0; i < param.chargeRange; i++) {
            binOffsets[i] = (long) round((mzBinMinValue - filter[i] - massBinMinValue) * param.binWidth);
        }

        auto mzBins = getMzBins(logMzPeaks, mzBinMinValue, mzBinMaxValue, param.binWidth);
        boost::dynamic_bitset<> massBins(massBinNumber);

        auto unionPrevMassBins = getUnionMassBin(massBins, massBinMinValue, prevMassBinVector, prevMinBinLogMassVector, param);
        auto perMassChargeRanges = getMassBins(massBins, mzBins, massBinMinValue,
                                               binOffsets,
                                               hBinOffsets,
                                               unionPrevMassBins,
                                               param);
        auto unionMassBins = unionPrevMassBins | massBins;
        auto peakGroups = getPeakGroupsWithMassBins(unionMassBins, massBins, logMzPeaks, mzBinMinValue,
                                                    binOffsets, perMassChargeRanges,
                                                    param);

        // prev mass update here.
        if (prevMassBinVector.size() >= (Size) param.numOverlappedScans) {
            prevMassBinVector.erase(prevMassBinVector.begin());
            prevMinBinLogMassVector.erase(prevMinBinLogMassVector.begin());
        }
        auto index = massBins.find_first();
        vector<Size> mb;
        while (index != massBins.npos) {
            mb.push_back(index);
            index = massBins.find_next(index);
        }
        prevMassBinVector.push_back(mb);
        prevMinBinLogMassVector.push_back(massBinMinValue);
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
            //double &massBinMinValue,
                                                long *binOffsets,
            //double *filter,
                                                Byte **chargeRanges,
                                                const Parameter &param) {
        double binWidth = param.binWidth;
        double tol = param.tolerance;
        int minCharge = param.minCharge;

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
            int maxChargeRange = maxChargeRanges[massBinIndex];
            int minChargeRange = minChargeRanges[massBinIndex];

            /*if (maxChargeRange == 0) {
                maxChargeRange = param.chargeRange - 1;
                minChargeRange = 0;
            }*/
            //vector<Size> isoMassBins;
            for (int j = minChargeRange; j <= maxChargeRange; j++) {
                int charge = j + minCharge;
                long &binOffset = binOffsets[j];
                auto &cpi = currentPeakIndex[j];
                double maxIntensity = 0.0;
                int maxIntensityPeakIndex = 0;

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

                if (maxIntensityPeakIndex > 0) {
                    double mz = logMzPeaks[maxIntensityPeakIndex].orgPeak->getMZ() - Constants::PROTON_MASS_U;
                    double &logMz = logMzPeaks[maxIntensityPeakIndex].logMz;
                    double isof = Constants::C13C12_MASSDIFF_U / charge / mz;

                    for (int d = -1; d <= 1; d += 2) { // negative then positive direction.
                        int peakIndex = maxIntensityPeakIndex + (d < 0 ? d : 0);
                        for (int i = 0; peakIndex >= 0 && peakIndex < logMzPeakSize; i++) {
                            double centerLogMz = logMz + isof * i * d;
                            double centerLogMzMin = centerLogMz - tol;
                            double centerLogMzMax = centerLogMz + tol;
                            bool isotopePeakPresent = false;
                            for (; peakIndex >= 0 && peakIndex < logMzPeakSize; peakIndex += d) {
                                double &observedLogMz = logMzPeaks[peakIndex].logMz;
                                if (observedLogMz < centerLogMzMin) { if (d < 0) { break; } else { continue; }}
                                if (observedLogMz > centerLogMzMax) { if (d < 0) { continue; } else { break; }}
                                isotopePeakPresent = true;
                                LogMzPeak p(*logMzPeaks[peakIndex].orgPeak, charge, i * d);
                                pg.push_back(p);
                                auto bin = getBinNumber(p.logMz, mzBinMinValue, binWidth) + binOffset;
                                if (massBinIndex != bin) {
                                    unionedMassBins[bin] = false;
                                    massBins[bin] = false;
                                }
                            }
                            if (!isotopePeakPresent) break;
                            if (d < 0) {
                                isoOff = -i > isoOff ? isoOff : -i;
                            }
                        }
                    }
                }
            }
            //cout<<pg.peaks.size()<<endl;
            if (!pg.peaks.empty()) {
                for (LogMzPeak &p : pg.peaks) {
                    //cout<<p.orgPeak->getMZ()<<" " << p.charge<< " " << p.isotopeIndex << endl;
                    p.isotopeIndex -= isoOff;
                }
                peakGroups.push_back(pg);
            }
            massBinIndex = unionedMassBins.find_next(massBinIndex);
        }
        return peakGroups;
    }

    /*
    void addPeaksAroundPeakIndexToPeakGroup(PeakGroup &peakGroup, vector<Size> &peakMassBins,
                                            int &currentPeakIndex, double &mzToMassOffset,
                                            vector<LogMzPeak> &logMzPeaks,
                                            int &charge, double &tol, double &binWidth) {
        double &logMz = logMzPeaks[currentPeakIndex].logMz;
        int logMzPeaksSize = (int) logMzPeaks.size();
        for (int i = -1; i < 2; i++) {
            int peakIndex = currentPeakIndex + i;
            if (peakIndex < 0 || peakIndex >= logMzPeaksSize) continue;
            if (abs(logMz - logMzPeaks[peakIndex].logMz) > tol) continue;

            LogMzPeak p(*logMzPeaks[peakIndex].orgPeak, charge, 0);
            peakGroup.push_back(p);

            if (i > 0) peakMassBins.push_back((Size) ((p.logMz + mzToMassOffset) * binWidth));
        }
    }

    void addIsotopicPeaksToPeakGroup(PeakGroup &peakGroup, vector<Size> &peakMassBins,
                                     int &currentPeakIndex, double &mzToMassOffset, vector<LogMzPeak> &logMzPeaks,
                                     int &charge, double &tol, double &binWidth, int &isoOff) {
        double mz = logMzPeaks[currentPeakIndex].orgPeak->getMZ() - Constants::PROTON_MASS_U;
        double &logMz = logMzPeaks[currentPeakIndex].logMz;

        int logMzPeaksSize = (int) logMzPeaks.size();
        double isof = Constants::C13C12_MASSDIFF_U / charge / mz;

        for (int d = -1; d <= 1; d += 2) { // negative then positive direction.
            int peakIndex = currentPeakIndex + d;
            for (int i = 1; peakIndex >= 0 && peakIndex < logMzPeaksSize; i++) {
                double isotopeLogMz = logMz + isof * i * d;

                double isotopeLogMzMin = isotopeLogMz - tol;
                double isotopeLogMzMax = isotopeLogMz + tol;
                bool isotopePeakPresent = false;
                while (peakIndex >= 0 && peakIndex < logMzPeaksSize) {
                    double &logMzForIsotopePeak = logMzPeaks[peakIndex].logMz;
                    peakIndex += d;
                    if (logMzForIsotopePeak < isotopeLogMzMin) { if (d < 0) { break; } else { continue; }}
                    if (logMzForIsotopePeak > isotopeLogMzMax) { if (d < 0) { continue; } else { break; }}
                    isotopePeakPresent = true;
                    LogMzPeak p(*logMzPeaks[peakIndex - d].orgPeak, charge, i * d);
                    peakGroup.push_back(p);

                    if (d > 0) peakMassBins.push_back((Size) ((p.logMz + mzToMassOffset) * binWidth));
                }
                if (!isotopePeakPresent) break;
                if (d < 0) {
                    isoOff = -i > isoOff ? isoOff : -i;
                }
            }
        }
    }
*/

    boost::dynamic_bitset<>
    getMzBins(vector<LogMzPeak> &logMzPeaks, double &mzBinMinValue, double &mzBinMaxValue, double binWidth) {
        Size binNumber =
                getBinNumber(mzBinMaxValue, mzBinMinValue, binWidth) +
                1;//(Size)((mzBinMaxValue - mzBinMinValue) / tol);
        boost::dynamic_bitset<> mzBins(binNumber);
        for (auto &p : logMzPeaks) {
            Size bi = getBinNumber(p.logMz, mzBinMinValue, binWidth);
            mzBins.set(bi);
            //if (bi > 0) mzBins.set(bi - 1);
        }
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
        int binScoreThreshold = param.minContinuousChargePeakPairCount;
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
        // int* minContinuousChargePeakPairCount = new int[param.chargeRange];
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
                       const Parameter &param) {

        long binThresholdMinMass = (long) getBinNumber(log(param.minMass), massBinMinValue, param.binWidth);
        //long binNumber = (long) massBins.size();


        boost::dynamic_bitset<> isQualified(massBins.size());
        Byte *continuousChargePeakPairCount = new Byte[massBins.size()];
        fill_n(continuousChargePeakPairCount, massBins.size(), 0);
        Byte *noneContinuousChargePeakPairCount = new Byte[massBins.size()];
        fill_n(noneContinuousChargePeakPairCount, massBins.size(), 0);

        getInitialMassBins(massBins, mzBins, isQualified, continuousChargePeakPairCount,
                                                         noneContinuousChargePeakPairCount, hBinOffsets, binOffsets,
                                                         param,
                //massBinMinValue,
                                                         binThresholdMinMass);

        auto perMassChargeRanges = getFinalMassBins(massBins, mzBins, isQualified, unionPrevMassBins,
                                                    continuousChargePeakPairCount, noneContinuousChargePeakPairCount,
                                                    binOffsets,
                                                    param,
                                                    binThresholdMinMass);
        delete[] noneContinuousChargePeakPairCount;
        delete[] continuousChargePeakPairCount;
        // isQualified.~dynamic_bitset();
        return perMassChargeRanges;
    }


     void getInitialMassBins(boost::dynamic_bitset<> &massBins,
                             boost::dynamic_bitset<> &mzBins,
                             boost::dynamic_bitset<> &isQualified,
                             Byte *continuousChargePeakPairCount,
                             Byte *noneContinuousChargePeakPairCount,
                             long **hBinOffsets,
                             long *binOffsets,
                             const Parameter &param,
            //double &massBinMinValue,
                             long &binStart) {

        boost::dynamic_bitset<> hasHarmony(massBins.size());
        int chargeRange = param.chargeRange;
        int hChargeSize = (int) param.hCharges.size();
        //double binWidth = param.binWidth;
        int minContinuousChargePeakPairCount = param.minContinuousChargePeakPairCount;
        long mzBinSize = (long) mzBins.size();
        long binEnd = (long) massBins.size();

        Byte *prevCharges = new Byte[massBins.size()];
        fill_n(prevCharges, massBins.size(), (Byte) (chargeRange + 2));
        //Byte *_continuousChargePeakPairCount = new Byte[massBins.size()];
        //fill_n(_continuousChargePeakPairCount, massBins.size(), 0);

        auto mzBinIndex = mzBins.find_first();
        while (mzBinIndex != mzBins.npos) {
            for (Byte j = 0; j < chargeRange; j++) {
                long massBinIndex = mzBinIndex + binOffsets[j];
                if (massBinIndex < binStart) continue;
                if (massBinIndex >= binEnd) break;
                if (hasHarmony[massBinIndex]) {
                    continue;
                }

                auto cd = prevCharges[massBinIndex] - j;
                prevCharges[massBinIndex] = j;

                if (cd == 1) {
                    auto &hbOffsets = hBinOffsets[j];
                    for (int k = 0; k < hChargeSize; k++) {
                        long hbi = mzBinIndex - hbOffsets[k];// + rand() % 100000 - 50000 ;

                        for (int i = -2; i <= 2; i++) {
                            auto bin = hbi + i;
                            if (bin < 0 || bin > mzBinSize) continue;
                            if (mzBins[bin]) {
                                hasHarmony[massBinIndex] = true;
                                isQualified[massBinIndex] = false;
                                break;
                            }
                        }
                        if (hasHarmony[massBinIndex]) break;
                    }
                    if (hasHarmony[massBinIndex]) continue;

                    if (++continuousChargePeakPairCount[massBinIndex] >= minContinuousChargePeakPairCount) {
                        isQualified[massBinIndex] = true;
                    }
                } else {
                    ++noneContinuousChargePeakPairCount[massBinIndex];
                    //if(maxChargeRanges[mzBinIndex] == 0) maxChargeRanges[mzBinIndex] = j;// std::max(maxChargeRanges[mzBinIndex], j);
                }

                //maxChargeRanges[mzBinIndex]= std::max(maxChargeRanges[mzBinIndex], j);
            }
            mzBinIndex = mzBins.find_next(mzBinIndex);
        }

        //hasHarmony.~dynamic_bitset();
        delete[] prevCharges;
        //return maxChargeRanges;
    }

    Byte **getFinalMassBins(boost::dynamic_bitset<> &massBins, boost::dynamic_bitset<> &mzBins,
                            boost::dynamic_bitset<> &isQualified,
                            boost::dynamic_bitset<> &unionPrevMassBins,
                            Byte *continuousChargePeakPairCount,
                            Byte *noneContinuousChargePeakPairCount,
                            long *binOffsets,
                            const Parameter &param,
                            long binStart) {

        int chargeRange = param.chargeRange;
        Byte *maxChargeRanges = new Byte[massBins.size()];
        fill_n(maxChargeRanges, massBins.size(), 0);

        Byte *minChargeRanges = new Byte[massBins.size()];
        fill_n(minChargeRanges, massBins.size(), 200);

        auto mzBinIndex = mzBins.find_first();
        const int minChargeScore = -1;
        long binEnd = (long) massBins.size();

        while (mzBinIndex != mzBins.npos) {
            long maxIndex = -1;
            int max = minChargeScore;
            Byte maxCharge = 0;

            for (Byte j = 0; j < chargeRange; j++) {
                long massBinIndex = mzBinIndex + binOffsets[j];
                if (massBinIndex < binStart) continue;
                if (massBinIndex >= binEnd) break;
                if (!isQualified[massBinIndex]) {
                    if(!unionPrevMassBins[massBinIndex])continue;
                }

                int t = continuousChargePeakPairCount[massBinIndex] - noneContinuousChargePeakPairCount[massBinIndex];//
                if (max <= t) {
                    max = t;
                    maxIndex = massBinIndex;
                    maxCharge = j;
                }
            }
            if (maxIndex > 0) {
                if(isQualified[maxIndex]) massBins[maxIndex] = true;
                maxChargeRanges[maxIndex] = std::max(maxChargeRanges[maxIndex], maxCharge);
                minChargeRanges[maxIndex] = std::min(minChargeRanges[maxIndex], maxCharge);
            }
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
        for (auto &pg : peakGroups) {
            auto perChargeIntensities = new double[param.chargeRange];
            auto perIsotopeIntensities = new double[param.maxIsotopeCount];
            updatePerChargeIsotopeIntensities(perChargeIntensities, perIsotopeIntensities, pg, param);
            if (!isIntensitiesQualified(perChargeIntensities, perIsotopeIntensities, param)) {
                continue;
            }

            pg.chargeDistributionScore = getChargeDistributionScore(perChargeIntensities, param);
            if (pg.chargeDistributionScore < param.chargeDistributionScoreThreshold) {
                continue;
            }

            double monoIsotopeMass;
            pg.isotopeCosineScore = getIsotopeCosineAndDetermineExactMass(monoIsotopeMass, pg, perIsotopeIntensities,
                                                                          averagines);

            if (pg.isotopeCosineScore <= param.isotopeCosineThreshold) {
                continue;
            }
            pg.monoisotopicMass = monoIsotopeMass;
            pg.intensity = accumulate(perChargeIntensities, perChargeIntensities + param.chargeRange, .0);
            filteredPeakGroups.push_back(pg);
            intensities.push_back(pg.intensity);
        }
        if (filteredPeakGroups.empty()) return filteredPeakGroups;
        filterPeakGroupsByIntensity(peakGroups, intensities, param);
        return filteredPeakGroups;
    }

    void updatePerChargeIsotopeIntensities(double *perChargeIntensities, double *perIsotopeIntensities, PeakGroup &pg,
                                           const Parameter &param) {
        fill_n(perChargeIntensities, param.chargeRange, 0);
        fill_n(perIsotopeIntensities, param.maxIsotopeCount, 0);

        for (auto &p : pg.peaks) {
            if (p.isotopeIndex >= param.maxIsotopeCount) continue;
            int index = p.charge - param.minCharge;
            double intensity = p.orgPeak->getIntensity();
            perChargeIntensities[index] += intensity;
            perIsotopeIntensities[p.isotopeIndex] += intensity;
        }
    }

    bool isIntensitiesQualified(double *perChargeIntensities, double *perIsotopeIntensities, const Parameter &param) {
        int maxSetIntensityCounter = 0;
        int setIntensityCounter = 0;
        for (int i = 1; i < param.chargeRange; i++) {

            //if (perChargeIntensities[i] > 0 && perChargeIntensities[i-1] > 0) {
            //    setIntensityCounter ++;
            //    continue;
            //}
            if (perChargeIntensities[i] <= 0) {
                setIntensityCounter = 0;
                continue;
            }
            setIntensityCounter++;
            maxSetIntensityCounter =
                    maxSetIntensityCounter > setIntensityCounter ? maxSetIntensityCounter : setIntensityCounter;
        }
        if (maxSetIntensityCounter < param.minContinuousChargePeakPairCount) return false;

        maxSetIntensityCounter = 0;
        setIntensityCounter = 0;
        for (int i = 0; i < param.maxIsotopeCount; i++) {
            if (perIsotopeIntensities[i] <= 0) {
                setIntensityCounter = 0;
                continue;
            }
            setIntensityCounter++;
            maxSetIntensityCounter =
                    maxSetIntensityCounter > setIntensityCounter ? maxSetIntensityCounter : setIntensityCounter;
        }
        return maxSetIntensityCounter >= param.minContinuousIsotopeCount;
    }

    double getIsotopeCosineAndDetermineExactMass(double &monoIsotopeMass, PeakGroup &pg, double *perIsotopeIntensities,
                                                 PrecalcularedAveragine &averagines) {
        auto iso = averagines.get(pg.peaks[0].getMass());

        int isoSize = (int) iso.size();

        /*
        map<int, LogMzPeak> mostIntensePeakPerIsotopeMap;
        for (auto &p : pg.peaks) {
            int i = p.isotopeIndex;
            auto it = mostIntensePeakPerIsotopeMap.find(i);
            if(it != mostIntensePeakPerIsotopeMap.end()){
                auto op = it->second;
                if(op.orgPeak->getIntensity() < p.orgPeak->getIntensity()){
                    mostIntensePeakPerIsotopeMap[i] = p;
                }
            }else mostIntensePeakPerIsotopeMap[i] = p;
        }


        LogMzPeak pp;
        vector<double> diffs;
        for (auto it=mostIntensePeakPerIsotopeMap.begin(); it!=mostIntensePeakPerIsotopeMap.end(); ++it) {
            auto p = it->second;
            if(pp.isotopeIndex>=0 && p.isotopeIndex - pp.isotopeIndex == 1){
                diffs.push_back((p.getMass() - pp.getMass()));
            }
            pp = p;
        }*/
        double isoDiff = Constants::C13C12_MASSDIFF_U;// OpenMS::Math::mean(diffs.begin(), diffs.end());

        int mostAbundantIndex = 0;
        auto mostAbundant = iso.getMostAbundant();
        for (int i = 0; i < isoSize; i++) {
            if (iso[i] == mostAbundant) {
                mostAbundantIndex = i;
                break;
            }
        }

        int offset = 0;
        double maxCosine = -1;
        int maxIsotopeIndex = 0, minIsotopeIndex = isoSize;

        for (auto &p : pg.peaks) {
            maxIsotopeIndex = p.isotopeIndex < maxIsotopeIndex ? maxIsotopeIndex : p.isotopeIndex;
            minIsotopeIndex = p.isotopeIndex < minIsotopeIndex ? p.isotopeIndex : minIsotopeIndex;
        }
        for (int f = -mostAbundantIndex + 1; f <= 3; f++) {
            if (minIsotopeIndex < f) continue;
            if (maxIsotopeIndex - f > isoSize)continue;

            auto cos = getCosine(perIsotopeIntensities, iso, isoSize, f);

            if (maxCosine < cos) {
                maxCosine = cos;
                offset = f;
            }
        }
        for (auto &p : pg.peaks) {
            p.isotopeIndex -= offset;
        }

        double maxIntensityForMonoIsotopeMass = -1;
        for (auto &p : pg.peaks) {
            if (p.isotopeIndex > maxIsotopeIndex - offset) continue;
            double intensity = p.orgPeak->getIntensity();
            if (maxIntensityForMonoIsotopeMass > intensity) continue;
            maxIntensityForMonoIsotopeMass = intensity;
            monoIsotopeMass = p.getMass() - p.isotopeIndex * isoDiff;
        }

        /*double minError = 100.0;\

        for (auto &p : pg.peaks) {
            //    if (p.isotopeIndex != -offset) continue;
            //double intensity = p.orgPeak->getIntensity();
            double m = (p.getMass() - p.isotopeIndex * isoDiff);
            double error = m - getNominalMass(m);
            error = error<0?-error:error;
            if(error>minError) continue;
            minError = error;
            monoIsotopeMass = m;
        }*/

        //monoIsotopeMass += offset * isoDiff;
        return maxCosine;
    }


    int getChargeDistributionScore(double *perChargeIntensities, const Parameter &param) {
        int maxIntensityIndex = 0;
        double maxIntensity = -1;
        for (int i = 0; i < param.chargeRange; i++) {
            if (maxIntensity < perChargeIntensities[i]) {
                maxIntensity = perChargeIntensities[i];
                maxIntensityIndex = i;
            }
        }

        int score = 0;
        for (int k = 1; k < param.chargeRange; k++) {
            int d1 = k <= maxIntensityIndex ? 0 : -1;
            int d2 = k <= maxIntensityIndex ? -1 : 0;
            double &int1 = perChargeIntensities[k + d1];
            double &int2 = perChargeIntensities[k + d2];
            //if (int2 <= 0) continue;
            if (int1 == int2) continue;
            if (int1 == 0) score -= 2;
            else score += int1 > int2 ? 1 : -1;
        }
        return score;
    }


    double getCosine(double *a, IsotopeDistribution &b, int &isoSize, int offset = 0) {
        double n = 0, d1 = 0, d2 = 0;
        for (int i = 0; i < isoSize; i++) {
            int j = i + offset;
            double bInt = b[i].getIntensity();
            d2 += bInt * bInt;
            if (j < 0 || j >= isoSize) continue;
            n += a[j] * bInt;
            d1 += a[j] * a[j];

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

        double delta = param.maxRTDelta;

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
