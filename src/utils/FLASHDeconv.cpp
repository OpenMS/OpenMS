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
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <QDirIterator>
#include <QFileInfo>
#include <chrono>
#include <OpenMS/FILTERING/DATAREDUCTION/MassTraceDetection.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

using namespace OpenMS;
using namespace std;

class FLASHDeconv:
        public TOPPBase {

public:
    FLASHDeconv() :
            TOPPBase("FLASHDeconv", "Ultra-fast high-quality deconvolution enables online processing of top-down MS data",
                     false)
 {}

    struct Parameter{
        int minCharge;
        double minMass;
        double maxMass;
        double tolerance; // up to here: ordinary user accessible parameters
        String fileName;

        double intensityThreshold;// advanced parameters
        int chargeRange;
        int minContinuousChargeCount;
        int minChargeCount;
        int minContinuousIsotopeCount;
        int maxIsotopeCount;
        int maxMassCount;
        double isotopeCosineThreshold;
        int chargeDistributionScoreThreshold;
        double maxRTDelta;
    };

    struct PrecalcularedAveragine{
        vector<IsotopeDistribution> isotopes;
        double massInterval;
        double minMass;

        PrecalcularedAveragine(double m, double M, double delta, CoarseIsotopePatternGenerator *generator): massInterval(delta), minMass(m)
        {
            int i = 0;
            while(true){
                double a = i * massInterval;
                i++;
                if(a < m) continue;
                if(a > M) break;
                auto iso = generator->estimateFromPeptideWeight(a);
                iso.trimRight(.01*iso.getMostAbundant().getIntensity());
                isotopes.push_back(iso);
            }
        }

        IsotopeDistribution get(double mass){
            Size i = (Size)((mass - minMass)/massInterval);
            i = i >= isotopes.size() ? isotopes.size()-1 : i;
            return isotopes[i];
        }
    } ;

    struct LogMzPeak {
        Peak1D *orgPeak;
        double logMz;
        double mass = .0;
        int charge;
        int isotopeIndex = -1;
        double score;

        LogMzPeak() : orgPeak(nullptr), logMz(-10000), charge(0), isotopeIndex(0), score(0) {}
        explicit LogMzPeak(Peak1D &peak) : orgPeak(&peak), logMz(getLogMz(peak.getMZ())), charge(0), isotopeIndex(0),
                                           score(0) {}
        LogMzPeak(Peak1D &peak, int c, int i) : orgPeak(&peak), logMz(getLogMz(peak.getMZ())), charge(c),
                                                isotopeIndex(i), score(0) {}
        double getMass() {
            if(mass <= 0) mass = exp(logMz) * charge;
            return mass;
        }
        /*double getMonoIsotopeMass(){
            return getMass() - isotopeIndex * Constants::C13C12_MASSDIFF_U;
        }*/

        bool operator<(const LogMzPeak &a) const {
            return this->logMz < a.logMz;
        }
        bool operator>(const LogMzPeak &a) const {
            return this->logMz >= a.logMz;
        }

    };

    struct PeakGroup{
        vector<LogMzPeak> peaks;
        double monoisotopicMass = .0;
        double intensity = .0;
        int chargeDistributionScore = 0;
        double isotopeCosineScore = .0;
        int massIndex, specIndex;
        MSSpectrum *spec;

        void push_back(LogMzPeak & p){
            peaks.push_back(p);
        }

        void reserve(Size n){
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

    static int getNominalMass(double &m){
        return (int)(m*0.999497+.5);
    }

    void registerOptionsAndFlags_() override {
        registerInputFile_("in", "<inputfile>", "", "Input file");
        //setValidFormats_("in", ListUtils::create<String>("mzML"));
        registerOutputFile_("out", "<output file prefix/output dir>", "", "Output file prefix or output dir (if prefix, [file prefix].tsv , [file prefix]feature.tsv, and [file prefix].m will be generated. "
                                                        "if dir, [dir]/[inputfile].tsv, [dir]/[inputfile]feature.tsv, and [dir]/[inputfile].m are generated per [inputfile])");

        registerIntOption_("minC", "<max charge>", 2, "minimum charge state", false, false);
        registerIntOption_("maxC", "<min charge>", 35, "maximum charge state", false, false);
        registerDoubleOption_("minM", "<min mass>", 500.0, "minimum mass (Da)", false, false);
        registerDoubleOption_("maxM", "<max mass>", 100000.0, "maximum mass (Da)", false, false);
        registerDoubleOption_("tol", "<tolerance>", 5.0, "ppm tolerance", false, false);

        registerDoubleOption_("minInt", "<min intensity>", 100.0, "intensity threshold", false, true);
        registerIntOption_("minCCC", "<min continuous charge count>", 3, "minimum number of peaks of continuous charges per mass", false, true);
        registerIntOption_("minCC", "<min charge count>", 8, "minimum number of peaks of distinct charges per mass (recommended - ~25% of (maxC - minC))", false, true);
        registerIntOption_("minIC", "<min isotope count>", 3, "minimum continuous isotope count", false, true);
        registerIntOption_("maxIC", "<max isotope count>", 100, "maximum isotope count", false, true);
        registerIntOption_("maxMC", "<max mass count>", -1, "maximum mass count per spec", false, true);
        registerDoubleOption_("minIsoScore", "<score 0-1>", .75, "minimum isotope cosine score threshold (0-1)", false, true);
        registerIntOption_("minCDScore", "<score 0,1,2,...>", 1, "minimum charge distribution score threshold (>= 0)", false, true);
        registerDoubleOption_("maxRTDelta", "<max RT delta>", -1, "max retention time duration with no peak in a feature (seconds); if negative, no feature finding performed", false, true);
    }

    Parameter setParameter(){
        Parameter param;
        param.minCharge = getIntOption_("minC");
        param.chargeRange = getIntOption_("maxC") - param.minCharge;
        param.maxMass = getDoubleOption_("maxM");
        param.minMass = getDoubleOption_("minM");
        param.tolerance = getDoubleOption_("tol")*1e-6;

        param.intensityThreshold = getDoubleOption_("minInt");
        param.minContinuousChargeCount = getIntOption_("minCCC");
        param.minChargeCount = getIntOption_("minCC");
        param.minContinuousIsotopeCount = getIntOption_("minIC");
        param.maxIsotopeCount = getIntOption_("maxIC");
        param.maxMassCount = getIntOption_("maxMC");
        param.isotopeCosineThreshold = getDoubleOption_("minIsoScore");
        param.chargeDistributionScoreThreshold = getIntOption_("minCDScore");
        param.maxRTDelta = getDoubleOption_("maxRTDelta");
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
        QString path = QString::fromUtf8( infilePath.data(), (int)infilePath.size());
        QFileInfo check_file(path);

        if (check_file.isDir()){
            QDirIterator it(path, QStringList() << "*.mzml", QDir::Files, QDirIterator::Subdirectories);
            while (it.hasNext()) {
                infileArray.push_back( it.next());
            }
        }else{
            infileArray.push_back(path);
        }
        sort(infileArray.begin(), infileArray.end());

        bool isOutPathDir = (QFileInfo(QString::fromUtf8(outfilePath.data(), (int)outfilePath.size())).isDir());

        cout << "Initializing ... " <<endl;
        auto param = setParameter();
        auto averagines = getPrecalculatedAveragines(param);
        int specCntr = 0, qspecCntr = 0, massCntr = 0, featureCntr = 0;
        int total_specCntr = 0, total_qspecCntr = 0, total_massCntr = 0;//, total_featureCntr = 0;
        double total_elapsed_cpu_secs = 0, total_elapsed_wall_secs = 0;
        fstream fs, fsm, fsf;

        if(!isOutPathDir) {
            fs.open(outfilePath + ".tsv", fstream::out);
            if(param.maxRTDelta>=0) fsf.open(outfilePath + "feature.tsv", fstream::out);

            writeHeader(fs, fsf, param.maxRTDelta>=0);
            fsm.open(outfilePath + ".m", fstream::out);
            fsm << "m=[";
        }

        for (auto& infile : infileArray){
            param.fileName = QFileInfo(infile).fileName().toStdString();

            if(isOutPathDir){
                std::string outfileName(param.fileName);
                std::size_t found = outfileName.find_last_of(".");
                outfileName = outfileName.substr(0,found);
                fs.open(outfilePath + outfileName + ".tsv", fstream::out);
                if(param.maxRTDelta>=0) fsf.open(outfilePath + "m" + outfileName + "feature.tsv", fstream::out);
                writeHeader(fs, fsf, param.maxRTDelta>=0);

                outfileName.erase(std::remove(outfileName.begin(), outfileName.end(), '_'), outfileName.end());
                outfileName.erase(std::remove(outfileName.begin(), outfileName.end(), '-'), outfileName.end());
                fsm.open(outfilePath + "m" + outfileName + ".m", fstream::out);
                fsm << "m=[";
            }

            double elapsed_cpu_secs = 0, elapsed_wall_secs = 0;
            cout << "Processing : " << infile.toStdString() << endl;

            MSExperiment map;
            MzMLFile mzml;
            mzml.setLogType(log_type_);
            mzml.load(infile, map);

            cout << "Running FLASHDeconv ... "<<endl;
            auto begin = clock();
            auto t_start = chrono::high_resolution_clock::now();
            auto peakGroups = Deconvolution(map, param, averagines, specCntr, qspecCntr, massCntr);
            auto t_end = chrono::high_resolution_clock::now();
            auto end = clock();
            elapsed_cpu_secs = double(end - begin) / CLOCKS_PER_SEC;
            elapsed_wall_secs = chrono::duration<double>(t_end-t_start).count();

            /*cout << endl << "-- done [took " << elapsed_cpu_secs << " s (CPU), " <<elapsed_wall_secs << " s (Wall)] --" << endl;
            cout << "-- per spectrum [took " << 1000.0*elapsed_cpu_secs/(specCntr - prevSpecCntr)
                << " ms (CPU), " << 1000.0*elapsed_wall_secs /(specCntr - prevSpecCntr) << " ms (Wall)] --"  << endl;
            cout << "Found " << massCntr - prevMassCntr << " masses in "<< qspecCntr - prevQspecCntr << " MS1 spectra out of "
            << specCntr - prevSpecCntr << endl;
            */
            if(peakGroups.empty()) continue;
            cout<< "writing results ...";
            cout.flush();
            for(auto &pg : peakGroups)
                writePeakGroup(pg, peakGroups.size(), param, fs, fsm);
            cout<<"done\n";
/*
            Param common_param = getParam_().copy("algorithm:common:", true);
            writeDebug_("Common parameters passed to sub-algorithms (mtd and ffm)", common_param, 3);

            Param mtd_param = getParam_().copy("algorithm:mtd:", true);
            writeDebug_("Parameters passed to MassTraceDetection", mtd_param, 3);

            MassTraceDetection mtdet;
            mtd_param.insert("", common_param);
            mtd_param.remove("chrom_fwhm");

            // trace_termination_outliers : MassTraceDetection.cpp // * (param.minCharge + param.chargeRange)
            mtd_param.setValue("mass_error_ppm", param.tolerance * (param.chargeRange + param.minCharge), "Allowed mass deviation (in ppm).");
            mtd_param.setValue("trace_termination_criterion", "sample_rate", "");
            mtd_param.setValue("reestimate_mt_sd", "false","");
            mtd_param.setValue("quant_method", "area", "");
            mtd_param.setValue("noise_threshold_int", .0, "");

            double rtd = (map[map.size()-1].getRT() - map[0].getRT() ) / (specCntr - prevSpecCntr);
            double delta = param.maxRTDelta < 0? (map[map.size()-1].getRT() - map[0].getRT() )/10.0 : param.maxRTDelta;
            mtd_param.setValue("min_sample_rate", rtd/delta, "");
            mtd_param.setValue("min_trace_length", .001, "");

            mtdet.setParameters(mtd_param);

            vector<MassTrace> m_traces;
            mtdet.run(map, m_traces);  // m_traces : output of this function
            //cout<<map.size();
            for(auto &mt : m_traces){
                auto mass = mt.getCentroidMZ();
                fsf<< fixed << setprecision(4)<<++featureCntr<<"\t"<<param.fileName<<"\t"<<mass<<"\t"<<mt.begin()->getRT()
                <<"\t"<<mt.rbegin()->getRT()<<"\t"<<mt.getTraceLength()<<"\t"<< mt[mt.findMaxByIntPeak()].getRT()
                <<"\t"<<mt.getMaxIntensity(false)<<"\t"<<mt.computePeakArea()<<"\n";
            }
*/
            if(param.maxRTDelta>=0) findNominalMassFeatures(peakGroups, featureCntr, fsf, param);

            if(isOutPathDir){
                fsm << "];";
                fsm.close();
                fs.close();
                if(param.maxRTDelta>=0) fsf.close();
                total_specCntr += specCntr;
                total_qspecCntr += qspecCntr;
                total_massCntr += massCntr;
                specCntr = qspecCntr = massCntr = featureCntr = 0;
            }else{
                total_specCntr = specCntr;
                total_qspecCntr = qspecCntr;
                total_massCntr = massCntr;
            }

            //prevFeatureCntr = featureCntr;
            //prevSpecCntr = specCntr; prevQspecCntr = qspecCntr; prevMassCntr = massCntr;
            total_elapsed_cpu_secs += elapsed_cpu_secs; total_elapsed_wall_secs += elapsed_wall_secs;
        }

        //if(infileArray.size() > 1){
            cout << "-- done [took " << total_elapsed_cpu_secs << " s (CPU), " <<total_elapsed_wall_secs << " s (Wall)] --" << endl;
            cout << "-- per spectrum [took " << 1000.0*total_elapsed_cpu_secs/total_specCntr
                 << " ms (CPU), " << 1000.0*total_elapsed_wall_secs /total_specCntr << " ms (Wall)] --" << endl;
            cout << "In total, found " << total_massCntr << " masses in "<< total_qspecCntr << " MS1 spectra out of "
                 << total_specCntr << endl;
       // }

        if(!isOutPathDir) {
            fsm << "];";
            fsm.close();
            fs.close();
            if(param.maxRTDelta>=0) fsf.close();
        }
        return EXECUTION_OK;
    }

    void writeHeader(fstream &fs, fstream &fsf, bool featureOut = false){
        fs << "MassIndex\tSpecIndex\tFileName\tSpecID\tMassCountInSpec\tExactMass\tNominalMass(round(ExactMass*0.999497))\t"
              "AggregatedIntensity\tRetentionTime\tPeakCount\tPeakMZs\tPeakCharges\tPeakMasses\tPeakIsotopeIndices\t"
              "PeakIntensities\tChargeDistScore\tIsotopeCosineScore\n";
        if(!featureOut) return;
        fsf << "FeatureID\tFileName\tNominalMass\tStartRetentionTime"
                   "\tEndRetentionTime\tRetentionTimeDuration\tApexRetentionTime"
                   "\tApexIntensity\tAbundance"
                << endl;

        return;
    }

    PrecalcularedAveragine getPrecalculatedAveragines(Parameter &param){
        auto generator = new CoarseIsotopePatternGenerator();
        auto maxIso = generator->estimateFromPeptideWeight(param.maxMass);
        maxIso.trimRight(.01*maxIso.getMostAbundant().getIntensity());
        param.maxIsotopeCount = min(param.maxIsotopeCount, (int)maxIso.size() - 1);
        generator->setMaxIsotope((Size)param.maxIsotopeCount);
        return PrecalcularedAveragine(param.minMass, param.maxMass, max(10.0, (param.maxMass - param.minMass)/100.0), generator);
    }

    vector<PeakGroup> Deconvolution(MSExperiment &map, Parameter &param,
            PrecalcularedAveragine& averagines, int &specCntr, int &qspecCntr, int &massCntr) {

        double* filter = new double[param.chargeRange];
        double* harmonicFilter = new double[param.chargeRange];
        double* harmonicFilter2 = new double[param.chargeRange];

        for (int i = 0; i < param.chargeRange; i++) {
            filter[i] = log(1.0 / (i + param.minCharge)); // should be descending, and negative!
            harmonicFilter[i] = log(1.0 / (i + .5 + param.minCharge));
            harmonicFilter2[i] = log(1.0 / (i + 1.0/3.0 + param.minCharge));
        }

        float prevProgress = .0;
        vector<PeakGroup> allPeakGroups;
        //vector<PeakGroup> prevPgs;

        for (auto it = map.begin(); it != map.end(); ++it) {
            if (it->getMSLevel() != 1) continue;

            float progress = (float)(it - map.begin())/map.size();
            if(progress>prevProgress+.05){
                printProgress(progress);
                prevProgress = progress;
            }

            specCntr++;

            auto logMzPeaks = getLogMzPeaks(*it, param.intensityThreshold);
            auto peakGroups = getPeakGroupsFromSpectrum(logMzPeaks, filter, harmonicFilter, harmonicFilter2, param);
            scoreAndFilterPeakGroups(peakGroups, averagines, param);

            //it->clear(false);
            if(peakGroups.empty()){
                continue;
            }
            //prevPgs = peakGroups;
            qspecCntr++;
            for (auto &pg : peakGroups) {
                massCntr++;
                pg.spec = &(*it);
                pg.massIndex = massCntr;
                pg.specIndex = qspecCntr;
                allPeakGroups.push_back(pg);

                /*
                for (auto &p : pg.peaks) {
                    auto mass = p.getMass() - p.isotopeIndex * Constants::C13C12_MASSDIFF_U;
                    //mass = getNominalMass(mass);
                    Peak1D tp(mass, (float) p.orgPeak->getIntensity());//
                    it->push_back(tp);
                    //Peak1D tp2(mass+.0000001, (float) p.orgPeak->getIntensity());
                    //it->push_back(tp2);
                }
                auto mass = pg.monoisotopicMass;
                mass = getNominalMass(mass);
                Peak1D tp(mass, (float) pg.intensity);//
                it->push_back(tp);*/
            }
            //it->sortByPosition();
        }
        printProgress(1);
        return allPeakGroups;
    }

    void writePeakGroup(PeakGroup &pg, Size peakGroupSize, Parameter& param, fstream &fs, fstream &fsm){
        double m = pg.monoisotopicMass;
        double intensity = pg.intensity;
        int nm = getNominalMass(m);

        fs<<fixed<<setprecision(4);

        fs <<pg.massIndex<<"\t"<<pg.specIndex<<"\t"<<param.fileName<<"\t"<<pg.spec->getNativeID()<<"\t"<<peakGroupSize<<"\t"
            << m << "\t" << nm<<"\t"<< intensity<<"\t"<<pg.spec->getRT()<<"\t"<<pg.peaks.size()<<"\t";
        sort(pg.peaks.begin(), pg.peaks.end());


        for(auto &p : pg.peaks){
            fs<<p.orgPeak->getMZ()<<";";
        }
        fs<<"\t";
        for(auto &p : pg.peaks){
            fs<<p.charge<<";";
        }
        fs<<"\t";
        for(auto &p : pg.peaks){
            fs<<p.getMass()<<";";
        }
        fs<<"\t";
        for(auto &p : pg.peaks){
            fs<<p.isotopeIndex<<";";
        }
        fs<<"\t";
        for(auto &p : pg.peaks){
            fs<<p.orgPeak->getIntensity()<<";";
        }
        fs<<"\t"<<pg.chargeDistributionScore<<"\t"<<pg.isotopeCosineScore<<"\n";

        fsm<<fixed << setprecision(4) << m<<","<<nm<<","<< intensity<<","<<pg.spec->getRT()<<"\n";
    }

    void printProgress(float progress){
        int barWidth = 70;
        cout << "[";
        int pos = (int)(barWidth * progress);
        for (int i = 0; i < barWidth; ++i) {
            if (i < pos) std::cout << "=";
            else if (i == pos) std::cout << ">";
            else std::cout << " ";
        }
        cout << "] " << int(progress * 100.0) << " %\r";
        cout.flush();
    }

    vector<LogMzPeak> getLogMzPeaks(MSSpectrum &spec, double &threshold){
        vector<LogMzPeak> logMzPeaks;
        logMzPeaks.reserve(spec.size());
        for (auto &peak : spec) {
            if (peak.getIntensity() <= threshold) continue;
            LogMzPeak logMzPeak(peak);
            logMzPeaks.push_back(logMzPeak);
        }

        //vector<LogMzPeak> tlogMzPeaks;
        //tlogMzPeaks.reserve(spec.size() + prev.size());
        //for(auto &peak : logMzPeaks) tlogMzPeaks.push_back(peak);
        //for(auto &peak : prev) tlogMzPeaks.push_back(peak);
        //sort(tlogMzPeaks.begin(), tlogMzPeaks.end());

        //prev = logMzPeaks;
        return logMzPeaks;
    }


    vector<PeakGroup> getPeakGroupsFromSpectrum(vector<LogMzPeak> &logMzPeaks, double *filter, double *hFilter,double *hFilter2,
                                                const Parameter &param){
        double maxBinLogMass = std::min(logMzPeaks[logMzPeaks.size() - 1].logMz - filter[param.chargeRange - param.minChargeCount], log(param.maxMass));
        double minBinLogMass = logMzPeaks[0].logMz - filter[param.minContinuousChargeCount];

        Size massBinNumber = (Size)((maxBinLogMass - minBinLogMass) / param.tolerance + 1);
        boost::dynamic_bitset<> mzBins = getMzBins(logMzPeaks, param.tolerance);
        boost::dynamic_bitset<> massBins(massBinNumber);
        //Byte *massBinScores = new Byte[massBinNumber];
        int chargeRange = getMassBins(massBins, mzBins, minBinLogMass, logMzPeaks[0].logMz, filter, hFilter, hFilter2, massBinNumber, param);

        //double binScoreThreshold = calculateBinScoreThreshold(massBins, massBinScores, param);
        //auto c = param.chargeRange;
        return getPeakGroupsWithMassBins(massBins, logMzPeaks, minBinLogMass, filter, chargeRange, param);
    }

    vector<PeakGroup> getPeakGroupsWithMassBins(boost::dynamic_bitset<> &massBins,//double &binScoreThreshold,
                                                vector<LogMzPeak> &logMzPeaks,
                                                double &minBinLogMass,
                                                double *filter, int& chargeRange, const Parameter &param) {
        int *currentPeakIndex = new int[chargeRange];
        double tol = param.tolerance;

        fill_n(currentPeakIndex, chargeRange, 0);

        vector<PeakGroup> peakGroups;
        vector<Size> peakMassBins;
        peakMassBins.reserve(logMzPeaks.size());

        Size binThreshold = (Size)max(.0, (log(param.minMass) - minBinLogMass)/param.tolerance);
        //peakGroups.reserve((Size)(param.maxMassCount * 10));
        auto setBinIndex = massBins.find_first();
        while (setBinIndex != massBins.npos) {
            if(setBinIndex < binThreshold){
               setBinIndex = massBins.find_next(setBinIndex);
               continue;
           }
           // if (massBinScores[setBinIndex] >= binScoreThreshold) {
            int isoOff = 0, maxi = 0;
            PeakGroup pg;
            peakMassBins.clear();
            pg.reserve((Size)(param.chargeRange * 10));

            for (int j = 0; j < chargeRange; j++) {
                int charge = j + param.minCharge;
                double mzToMassOffset = - filter[j] - minBinLogMass;

                while (currentPeakIndex[j] < (int)logMzPeaks.size()) {
                    Size bi = (Size)((logMzPeaks[currentPeakIndex[j]].logMz - filter[j] - minBinLogMass) / tol);
                    if (bi > setBinIndex) break;
                    if (setBinIndex == bi) {
                        addPeaksAroundPeakIndexToPeakGroup(pg, peakMassBins,currentPeakIndex[j],
                                mzToMassOffset, logMzPeaks, charge, tol);
                        addIsotopicPeaksToPeakGroup(pg, peakMassBins, currentPeakIndex[j],
                                mzToMassOffset, logMzPeaks, charge, tol, isoOff, maxi);
                    }
                    currentPeakIndex[j]++;
                }
            }

            if(maxi + isoOff >= param.minContinuousIsotopeCount - 1){
                for(auto& i :peakMassBins){
                    massBins[i] = false;
                }
                for (LogMzPeak &p : pg.peaks) {
                    p.isotopeIndex -= isoOff;
                }
                peakGroups.push_back(pg);
            }
            setBinIndex = massBins.find_next(setBinIndex);
        }
        return peakGroups;
    }

    void addPeaksAroundPeakIndexToPeakGroup(PeakGroup &peakGroup, vector<Size> &peakMassBins,
            int &currentPeakIndex, double &mzToMassOffset, vector<LogMzPeak> &logMzPeaks,
            int &charge, double &tol){
        double& logMz = logMzPeaks[currentPeakIndex].logMz;
        int logMzPeaksSize = (int)logMzPeaks.size();
        for(int i=-1; i < 2; i++){
            int peakIndex = currentPeakIndex+i;
            if(peakIndex <0 || peakIndex>=logMzPeaksSize) continue;
            if(abs(logMz - logMzPeaks[peakIndex].logMz) > tol) continue;
            addPeakToPeakGroup(logMzPeaks[peakIndex], charge, 0, mzToMassOffset, tol, peakGroup, peakMassBins, i>0);
        }
    }

    void addIsotopicPeaksToPeakGroup(PeakGroup &peakGroup, vector<Size> &peakMassBins,
                          int &currentPeakIndex, double &mzToMassOffset, vector<LogMzPeak> &logMzPeaks,
                          int &charge, double &tol, int &isoOff, int &maxi){
        double mz = logMzPeaks[currentPeakIndex].orgPeak->getMZ() - Constants::PROTON_MASS_U;
        int logMzPeaksSize = (int)logMzPeaks.size();

        for (int d = -1; d <= 1; d += 2) { // negative then positive direction.
            double logMz = logMzPeaks[currentPeakIndex].logMz;
            int peakIndex = currentPeakIndex + d;
            for (int i = 1; peakIndex>=0 && peakIndex < logMzPeaksSize ; i++) {
                double isotopeLogMz = logMz + Constants::C13C12_MASSDIFF_U * i * d / charge / mz;
                //double isotopeLogMzMin = logMz + (Constants::C13C12_MASSDIFF_U - 0.001) * d / charge / mz;

                bool isotopePeakPresent = false;
                //double ci = -1.0;
                while (peakIndex >= 0 && peakIndex < logMzPeaksSize) {
                    double& logMzForIsotopePeak = logMzPeaks[peakIndex].logMz;
                    peakIndex += d;
                    if (logMzForIsotopePeak < isotopeLogMz - tol){ if (d < 0){ break;} else{ continue;}}
                    if (logMzForIsotopePeak > isotopeLogMz + tol){ if (d < 0){ continue;}else { break;}}
                    isotopePeakPresent = true;
                    addPeakToPeakGroup(logMzPeaks[peakIndex-d], charge, i * d, mzToMassOffset, tol, peakGroup, peakMassBins, d>0);
                    //if(ci<logMzPeaks[peakIndex-d].orgPeak->getIntensity()){
                    //    ci = logMzPeaks[peakIndex-d].orgPeak->getIntensity();
                        //logMz = logMzPeaks[peakIndex-d].logMz;
                    //}
                }
                if (!isotopePeakPresent) break;
                if (d < 0) {
                    isoOff = -i > isoOff ? isoOff : -i;
                }else{
                    maxi = i>maxi?i:maxi;
                }
            }
        }
    }

    void addPeakToPeakGroup(LogMzPeak &logMzPeak, int &charge, int isotopeIndex, double &mzToMassOffset,
             double &tol, PeakGroup &peakGroup, vector<Size> &peakMassBins, bool updatePeakMassBins = true){
        LogMzPeak p(*logMzPeak.orgPeak, charge, isotopeIndex);
        if(updatePeakMassBins)
            peakMassBins.push_back((Size)((p.logMz +mzToMassOffset) / tol));
        peakGroup.push_back(p);
    }

    boost::dynamic_bitset<> getMzBins(vector<LogMzPeak> &logMzPeaks, double tol) {
        Size binNumber = (Size)((logMzPeaks[logMzPeaks.size() - 1].logMz  - logMzPeaks[0].logMz) / tol + 1);
        boost::dynamic_bitset<> mzBins(binNumber);

        for (auto &p : logMzPeaks) {
            Size bi = (Size) ((p.logMz - logMzPeaks[0].logMz) / tol);
            if (bi >= binNumber-1) break;
            mzBins[bi] = true;
            mzBins[bi+1] = true;
        }
        //mzBins |= mzBins << 1;
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
        int binScoreThreshold = param.minContinuousChargeCount;
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
    int getMassBins(boost::dynamic_bitset<> &massBins, boost::dynamic_bitset<> &mzBins,
                     double &minBinLogMass, double &minLogMz,
                     double *filter, double *harmonicFilter,  double *harmonicFilter2, Size &binNumber, const Parameter& param) {

        long *binOffsets = new long[param.chargeRange];
        long *harmonicBinOffsets = new long[param.chargeRange];
        long *harmonicBinOffsets2 = new long[param.chargeRange];

        for (int i = 0; i < param.chargeRange; i++) {
            binOffsets[i] = (long)((minLogMz-filter[i]-minBinLogMass) / param.tolerance);
            harmonicBinOffsets[i] = (long)((minLogMz-harmonicFilter[i]-minBinLogMass) / param.tolerance);
            harmonicBinOffsets2[i] = (long)((minLogMz-harmonicFilter2[i]-minBinLogMass) / param.tolerance);
        }

        int chargeRange = getInitialMassBinsUsingMinChargeCount(massBins, mzBins, binOffsets, binNumber, minBinLogMass, param);
        return getFinalMassBinsUsingMinContinuousChargeCount(massBins, mzBins, chargeRange, binOffsets,harmonicBinOffsets, harmonicBinOffsets2, param);
    }

    int getInitialMassBinsUsingMinChargeCount(boost::dynamic_bitset<> &massBins, boost::dynamic_bitset<> &mzBins,
                             long *binOffsets, Size &binNumber, double &minBinLogMass, const Parameter &param){
        auto setBinIndex = mzBins.find_first();
        int chargeRange = param.chargeRange;
        int minChargeCount = param.minChargeCount;
        int minContinuousChargeCount = param.minContinuousChargeCount;
        long biThreshold = (long) binNumber;
        double binThresholdSmallMass = (log(15000.0)-minBinLogMass)/param.tolerance;
        //auto firstSetBinIndex = 1;
        //long lastSetBinIndex = (long)mzBins.size()-1;

        Byte *massBinScores = new Byte[binNumber];
        fill_n(massBinScores, binNumber, 0);

       // Byte *massBinScoreThresholds = new Byte[binNumber];
       // fill_n(massBinScoreThresholds, binNumber, 0);

        int maxChargeRange = 0;
        while (setBinIndex != mzBins.npos) {
            for (int j = 0; j < chargeRange; j++) {
                long bi = setBinIndex + binOffsets[j];
                if(bi < 0)continue;
                if(bi >= biThreshold) break;
                ++massBinScores[bi];
                auto &bs = massBinScores[bi];
                if(bs < minContinuousChargeCount) continue;
                if(bs < minChargeCount) {
                    if(bi>binThresholdSmallMass)continue;
                    if(bs < bi/binThresholdSmallMass * minChargeCount) continue;
                }
                massBins[bi] = true;
                maxChargeRange = maxChargeRange > j? maxChargeRange : j;
            }
            setBinIndex = mzBins.find_next(setBinIndex);
        }

        /*
        setBinIndex = massBins.find_first();
        while (setBinIndex != massBins.npos) {
            auto &bs = massBinScores[setBinIndex];
            if(bs >= minChargeCount) {
                setBinIndex = massBins.find_next(setBinIndex);
                continue;
            }
            int k = 0;
            while (k < chargeRange && setBinIndex > binOffsets[k]) {
                k++;
            }

            int l = chargeRange - 1;
            while (l >=0 && setBinIndex < lastSetBinIndex + binOffsets[l]) {
                l--;
            }

            int s = (k-l) * .8;
            if(bs >= s) {
                setBinIndex = massBins.find_next(setBinIndex);
                continue;
            }

            massBins[setBinIndex] = false;
            setBinIndex = massBins.find_next(setBinIndex);
        }*/


        return maxChargeRange;
    }

    int getFinalMassBinsUsingMinContinuousChargeCount(boost::dynamic_bitset<> &massBins, boost::dynamic_bitset<> &mzBins,
                                                    int& chargeRange, long *binOffsets, long *harmonicBinOffsets, long *harmonicBinOffsets2,
                                                    const Parameter &param){
        auto setBinIndex = massBins.find_first();
        //int chargeRange = param.chargeRange;
        int minContinuousCharge = param.minContinuousChargeCount;

        //Size threshold = (Size)max(.0, (log(param.minMass) - minBinLogMass)/param.tolerance);
        int maxChargeRange = 0;
        while (setBinIndex != massBins.npos) {
            Byte continuousChargeCntr = 0;
            //Byte maxContinuousChargeCntr = 0;
            massBins[setBinIndex] = false;
            for (int j = 0; j <= chargeRange; j++) {
                long bi = setBinIndex - binOffsets[j];
                if(bi<0) break;
                bool mzBinSet = (bi < (long)mzBins.size()) && mzBins[bi];
                if(!mzBinSet){
                    continuousChargeCntr = 0;
                    continue;
                }

                long hbi = setBinIndex - harmonicBinOffsets[j];
                long hbi2 = setBinIndex - harmonicBinOffsets2[j];

                bool harmonicMzBinClear = hbi<1 || hbi >= (long)mzBins.size()-1 || !(mzBins[hbi]||mzBins[hbi-1]||mzBins[hbi+1]);//||mzBins[hbi-1]||mzBins[hbi+1]
                if(!harmonicMzBinClear){
                    //massBins[setBinIndex] = false;
                    //continuousChargeCntr = 0;
                    break;
                }

                harmonicMzBinClear = hbi2<1 || hbi2 >= (long)mzBins.size()-1 || !(mzBins[hbi2]||mzBins[hbi2-1]||mzBins[hbi2+1]);//||mzBins[hbi2-1]||mzBins[hbi2+1]
                if(!harmonicMzBinClear){
                    //massBins[setBinIndex] = false;
                    //continuousChargeCntr = 0;
                    break;
                }

                bool set = ++continuousChargeCntr >= minContinuousCharge;
                if(!set) continue;
                massBins[setBinIndex] = true;
                //break;
                maxChargeRange = maxChargeRange > j? maxChargeRange : j;
                //maxContinuousChargeCntr = maxContinuousChargeCntr > continuousChargeCntr? maxContinuousChargeCntr : continuousChargeCntr;
            }
            setBinIndex = massBins.find_next(setBinIndex);
        }
        return maxChargeRange;
    }

    void scoreAndFilterPeakGroups(vector<PeakGroup> &peakGroups,
                                  PrecalcularedAveragine &averagines,
                                  const Parameter &param){
        vector<double> intensities;
        intensities.reserve(peakGroups.size());
        for (auto pg = peakGroups.begin(); pg != peakGroups.end(); ) {
            auto perChargeIntensities = new double[param.chargeRange];
            auto perIsotopeIntensities = new double[param.maxIsotopeCount];
            updatePerChargeIsotopeIntensities(perChargeIntensities, perIsotopeIntensities, *pg, param);
            if(!isIntensitiesQualified(perChargeIntensities, perIsotopeIntensities, param)){
                pg = peakGroups.erase(pg);
                continue;
            }

            pg->chargeDistributionScore = getChargeDistributionScore(perChargeIntensities, param);
            if(pg->chargeDistributionScore<param.chargeDistributionScoreThreshold){
                pg = peakGroups.erase(pg);
                continue;
            }

            double monoIsotopeMass;
            pg->isotopeCosineScore = getIsotopeCosine(monoIsotopeMass, *pg, perIsotopeIntensities, averagines);

            if (pg->isotopeCosineScore <= param.isotopeCosineThreshold){
                pg = peakGroups.erase(pg);
                continue;
            }
            pg->monoisotopicMass = monoIsotopeMass;
            pg->intensity = accumulate(perChargeIntensities, perChargeIntensities+param.chargeRange,.0);
            intensities.push_back(pg->intensity);
            ++pg;
        }
        filterPeakGroupsByIntensity(peakGroups, intensities, param);
    }

    void updatePerChargeIsotopeIntensities(double *perChargeIntensities, double *perIsotopeIntensities, PeakGroup &pg, const Parameter &param){
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

    bool isIntensitiesQualified(double *perChargeIntensities, double *perIsotopeIntensities, const Parameter &param){
        int setIntensityCounter = 0;
        int maxSetIntensityCounter = 0;
        for (int i = 0; i < param.chargeRange; i++) {
            if (perChargeIntensities[i] <= 0) {
                setIntensityCounter = 0;
                continue;
            }
            setIntensityCounter++;
            maxSetIntensityCounter =
                    maxSetIntensityCounter > setIntensityCounter ? maxSetIntensityCounter : setIntensityCounter;
        }
        if (maxSetIntensityCounter < param.minContinuousChargeCount) return false;

        setIntensityCounter = 0;
        maxSetIntensityCounter = 0;
        for (int i = 0; i < param.maxIsotopeCount; i++) {
            if (perIsotopeIntensities[i] <= 0) {
                setIntensityCounter = 0;
                continue;
            }
            setIntensityCounter++;
            maxSetIntensityCounter =
                    maxSetIntensityCounter > setIntensityCounter ? maxSetIntensityCounter : setIntensityCounter;
        }
        return maxSetIntensityCounter > param.minContinuousIsotopeCount;
    }

    double getIsotopeCosine(double &monoIsotopeMass, PeakGroup &pg, double *perIsotopeIntensities,
                            PrecalcularedAveragine &averagines){
        double maxIntensityForMonoIsotopeMass = -1;
        auto iso = averagines.get(pg.peaks[0].getMass());

        int isoSize = (int)iso.size();
        int maxIsotopeIndex=0, minIsotopeIndex = isoSize;

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
        //cout<<isoDiff<<endl;
        for (auto &p : pg.peaks) {
            maxIsotopeIndex = p.isotopeIndex < maxIsotopeIndex ? maxIsotopeIndex : p.isotopeIndex;
            minIsotopeIndex = p.isotopeIndex < minIsotopeIndex ? p.isotopeIndex : minIsotopeIndex;
            if (p.isotopeIndex != 0) continue;
            double intensity = p.orgPeak->getIntensity();
            if (maxIntensityForMonoIsotopeMass > intensity) continue;
            maxIntensityForMonoIsotopeMass = intensity;
            monoIsotopeMass = p.getMass() - p.isotopeIndex*isoDiff;
        }

        int mostAbundantIndex = 0;
        auto mostAbundant = iso.getMostAbundant();
        for(int i=0;i<isoSize;i++){
            if(iso[i] == mostAbundant){
                mostAbundantIndex = i;
                break;
            }
        }

        int offset = 0;
        double maxCosine = -1;
        for (int f = -mostAbundantIndex + 1; f <= 1; f++) {
            if(minIsotopeIndex < f) continue;
            if(maxIsotopeIndex - f > isoSize)continue;

            auto cos = getCosine(perIsotopeIntensities, iso, isoSize, f);

            if (maxCosine < cos) {
                maxCosine = cos;
                offset = f;
            }
        }
        for (auto &p : pg.peaks) {
            p.isotopeIndex -= offset;
        }
        monoIsotopeMass += offset * isoDiff;
        return maxCosine;
    }


    int getChargeDistributionScore(double *perChargeIntensities, const Parameter &param){
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
            int d1 = k<=maxIntensityIndex? 0 : -1;
            int d2 = k<=maxIntensityIndex? -1 : 0;
            double& int1 = perChargeIntensities[k + d1];
            double& int2 = perChargeIntensities[k + d2];
            if (int2 <= 0) continue;
            score += int1 >= int2 ? 1 : -1;
        }
        return score;
    }


    double getCosine(double *a, IsotopeDistribution &b, int& isoSize, int offset = 0) {
        double n = 0, d1 = 0, d2 = 0;
        for (int i = 0; i < isoSize; i++) {
            int j = i + offset;
            double bInt = b[i].getIntensity();
            d2 += bInt * bInt;
            if (j<0 || j >= isoSize) continue;
            n += a[j] * bInt;
            d1 += a[j] * a[j];

        }
        double d = sqrt(d1 * d2);
        if (d <= 0) return 0;
        return n / d;
    }


    void filterPeakGroupsByIntensity(vector<PeakGroup> &peakGroups, vector<double> &intensities, const Parameter &param){
        if(param.maxMassCount<0 || intensities.size() <= (Size)param.maxMassCount) return;
        Size mc = (Size)param.maxMassCount;
        sort(intensities.begin(), intensities.end());
        auto threshold = intensities[intensities.size() - mc];
        for (auto pg = peakGroups.begin(); pg != peakGroups.end(); ) {
            if(peakGroups.size() <= mc) break;
            if(pg->intensity < threshold) {
                pg = peakGroups.erase(pg);
                continue;
            }
            ++pg;
        }
    }



    void findNominalMassFeatures(vector<PeakGroup> &peakGroups, int &id, fstream& fs, const Parameter & param){

        map<int, vector<double>> massMap;
        map<int, vector<String>> writeMap;

        double delta = param.maxRTDelta;
        /*if(delta <=0){
            auto end = peakGroups[peakGroups.size()-1].spec->getRT();
            auto begin = peakGroups[0].spec->getRT();
            delta = (end-begin)/10.0;
        }*/
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
                if(v[1] - v[0] > 0 && v[1]-v[0]<delta) {
                    stringstream s;
                    s <<"\t" << param.fileName << "\t" << nm << "\t" << fixed << setprecision(5) << v[0] << "\t"
                      << v[1]<<"\t"<<v[1]-v[0] << "\t" << v[2] << "\t" << v[3]<< "\t" << v[4];
                    auto sit = writeMap.find(nm);
                    if (sit == writeMap.end()) {
                        vector<String> vs;
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
    }
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
