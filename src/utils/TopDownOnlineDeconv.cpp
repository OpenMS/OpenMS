//
// Created by JihyungKim on 07.12.18.
//

#include <OpenMS/config.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <queue>
#include "boost/dynamic_bitset.hpp"
#include <iostream>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <unordered_set>

using namespace OpenMS;
using namespace std;

class FlashDeconv:
        public TOPPBase {

public:
    FlashDeconv() :
            TOPPBase("FlashDeconv", "Real-time Deconvolution for Non-redundant MS2 acquisition with top down data",
                     false)\
 {}

    typedef struct Parameter{
        int minCharge;
        int maxCharge;
        double minMass;
        double maxMass;
        double tolerance; // up to here: ordinary user accessible parameters

        double intensityThreshold;
        int chargeRange;
        int minContinuousCharge;
        int minChargeCount;
        int minIsotopeCount;
        int maxMassCount;
        int maxIsotopeIndex;
        double isotopeCosineThreshold;
        double chargeDistributionScoreThreshold; // geek accessible parameters

        Parameter() : intensityThreshold(500), chargeRange(25), minCharge(5), maxCharge(minCharge + chargeRange), minContinuousCharge(3),
            minChargeCount(minContinuousCharge * 2), minIsotopeCount(4), maxMassCount(50), minMass(500.0), maxMass(50000.0),
            tolerance(5e-6), isotopeCosineThreshold(.5), chargeDistributionScoreThreshold(.0), maxIsotopeIndex(50) {}
    } Parameter;

    typedef struct LogMzPeak {
        Peak1D *orgPeak;
        double logMz;
        int charge;
        int isotopeIndex;
        double score;

        LogMzPeak() : orgPeak(nullptr), logMz(-10000), charge(0), isotopeIndex(0), score(0) {}
        explicit LogMzPeak(Peak1D &peak) : orgPeak(&peak), logMz(getLogMz(peak.getMZ())), charge(0), isotopeIndex(0),
                                           score(0) {}
        LogMzPeak(Peak1D &peak, int c, int i) : orgPeak(&peak), logMz(getLogMz(peak.getMZ())), charge(c),
                                                isotopeIndex(i), score(0) {}
        double getMass() {
            return exp(logMz) * charge;
        }
        double getMonoIsotopeMass(){
            return getMass() - isotopeIndex * Constants::C13C12_MASSDIFF_U;
        }
        bool operator<(const LogMzPeak &a) {
            return logMz < a.logMz;
        }
    } LogMzPeak;

    typedef struct PeakGroup{
        vector<LogMzPeak> peaks;

        void push_back(LogMzPeak & p){
            peaks.push_back(p);
        }

        void reserve(int n){
            peaks.reserve(n);
        }

    } PeakGroup;

protected:

    static double getLogMz(double mz) {
        return log(mz - Constants::PROTON_MASS_U);
    }


    void registerOptionsAndFlags_() override {
        registerInputFile_("in", "<file>", "", "Input file.");
        setValidFormats_("in", ListUtils::create<String>("mzML"));
    }

    ExitCodes main_(int, const char **) override {
        //-------------------------------------------------------------
        // parsing parameters
        //-------------------------------------------------------------
//        String infilePath = getStringOption_("in");

        String infileDir = "/Users/kyowonjeong/Documents/A4B/mzml/MS1only/yeast/";
        String infilePath = "/Users/kyowonjeong/Documents/A4B/mzml/MS1only/180523_Cytocrome_C_MS2_HCD_MS1only.mzML";
        infilePath = "/Users/kyowonjeong/Documents/A4B/mzml/MS1only/180523_Myoglobin_MS2_HCD_MS1only.mzML";
        String outfilePath = "/Users/kyowonjeong/Documents/A4B/matlab/myo.m";
        // just for quick use

        //-------------------------------------------------------------
        // reading input
        //-------------------------------------------------------------

        MzMLFile mzml;
        mzml.setLogType(log_type_);

        Parameter param;
        int specCntr = 0, qspecCntr = 0, massCntr = 0;

        double elapsed_secs = 0;
        for (int r = 1; r <= 2; r++) {
            ostringstream st;
            st << "/Users/kyowonjeong/Documents/A4B/matlab/yeast" << r << ".m";
            outfilePath = st.str();

            fstream fs;
            fs.open(outfilePath, fstream::out);
            fs << "m=[";
            for (int f = 1; f <= 12; f++) {
                infilePath = infileDir + "f"+f+"r" + r+".mzML";
                cout << "%file name : " << infilePath << endl;
                MSExperiment map;
                mzml.load(infilePath, map);
                cout << "%Loaded consensus maps" << endl;
                clock_t begin = clock();
                onlineDeconvolution(map, param, fs, specCntr, qspecCntr, massCntr);
                clock_t end = clock();
                elapsed_secs += double(end - begin) / CLOCKS_PER_SEC;
                std::cout << massCntr << " masses in "<< qspecCntr << " MS1 spectra deconvoluted so far" << endl;
            }
            fs << "];";
            fs.close();
        }

        std::cout << "%" << elapsed_secs << " seconds elapsed for " << specCntr << " MS1 spectra" << endl;
        std::cout << "%" << elapsed_secs / specCntr * 1000 << " msec per spectrum" << std::endl;

        return EXECUTION_OK;
    }

    vector<LogMzPeak> GetLogMzPeaks(MSSpectrum &spec, const Parameter &param){
        vector<LogMzPeak> logMzPeaks;
        logMzPeaks.reserve(spec.size());
        for (auto &peak : spec) {
            if (peak.getIntensity() <= param.intensityThreshold) continue;
            LogMzPeak logMzPeak(peak);
            logMzPeaks.push_back(logMzPeak);
        }
        return logMzPeaks;
    }

    double getIsotopeCosine(double &monoIsotopeMass, PeakGroup &pg, double *perIsotopeIntensities, CoarseIsotopePatternGenerator *generator, const Parameter &param){
        double maxIntensityForMonoIsotopeMass = -1;
        auto iso = generator->estimateFromPeptideWeight(pg.peaks[0].getMass());
        iso.trimRight(.1*iso.getMostAbundant().getIntensity());

        for (auto &p : pg.peaks) {
            if (p.isotopeIndex > iso.size()) continue;
            int index = p.charge - param.minCharge;
            double intensity = p.orgPeak->getIntensity();
            if (maxIntensityForMonoIsotopeMass > intensity) continue;
            maxIntensityForMonoIsotopeMass = intensity;
            monoIsotopeMass = p.getMonoIsotopeMass();
        }

        int offset = 0;
        double maxCosine = -1;
        for (int f = -3; f <= 3; f++) {
            auto cos = getCosine(perIsotopeIntensities, iso, f);

            if (maxCosine < cos) {
                maxCosine = cos;
                offset = f;
            }
        }

        monoIsotopeMass += offset * Constants::C13C12_MASSDIFF_U;
        return maxCosine;
    }

    void updatePerChargeIsotopeIntensities(double *perChargeIntensities, double *perIsotopeIntensities, PeakGroup &pg, const Parameter &param){
        fill_n(perChargeIntensities, param.chargeRange, 0);
        fill_n(perIsotopeIntensities, param.maxIsotopeIndex + 1, 0); // change it to be flexible

        for (auto &p : pg.peaks) {
            if (p.isotopeIndex > param.maxIsotopeIndex) continue;
            int index = p.charge - param.minCharge;
            double intensity = p.orgPeak->getIntensity();
            perChargeIntensities[index] += intensity;
            perIsotopeIntensities[p.isotopeIndex] += intensity;
        }
    }

    bool checkIntensities(double *perChargeIntensities, double *perIsotopeIntensities, const Parameter &param){
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
        if (maxSetIntensityCounter < param.minContinuousCharge) return false;

        setIntensityCounter = 0;
        maxSetIntensityCounter = 0;
        for (int i = 0; i < param.maxIsotopeIndex; i++) {
            if (perIsotopeIntensities[i] <= 0) {
                setIntensityCounter = 0;
                continue;
            }
            setIntensityCounter++;
            maxSetIntensityCounter =
                    maxSetIntensityCounter > setIntensityCounter ? maxSetIntensityCounter : setIntensityCounter;
        }
        if (maxSetIntensityCounter < param.minIsotopeCount) return false;
        return true;
    }

    double getChargeDistributionScore(double *perChargeIntensities, const Parameter &param){
        int maxIntensityIndex;
        double maxIntensity = -1;
        for (int i = 0; i < param.chargeRange; i++) {
            if (maxIntensity < perChargeIntensities[i]) {
                maxIntensity = perChargeIntensities[i];
                maxIntensityIndex = i;
            }
        }

        double score = 0;
        for (int k = 1; k < param.chargeRange; k++) {
            int d1 = k<=maxIntensityIndex? 0 : -1;
            int d2 = k<=maxIntensityIndex? -1 : 0;
            double int1 = perChargeIntensities[k + d1];
            double int2 = perChargeIntensities[k + d2];
            if (int2 <= 0) continue;
            score += int1 >= int2 ? 1 : -1;
        }
        return score;
    }

    std::map<double, double> getIntensityMassMap(vector<PeakGroup> &peakGroups, CoarseIsotopePatternGenerator *generator, const Parameter &param){
        std::map<double, double> intensityMassMap;

        for (auto &pg : peakGroups) {
            auto perChargeIntensities = new double[param.chargeRange];
            auto perIsotopeIntensities = new double[param.maxIsotopeIndex + 1];
            updatePerChargeIsotopeIntensities(perChargeIntensities, perIsotopeIntensities, pg, param);
            if(!checkIntensities(perChargeIntensities, perIsotopeIntensities, param)) continue;

            double chargeDistributionScore = getChargeDistributionScore(perChargeIntensities, param);
            if(chargeDistributionScore<=param.chargeDistributionScoreThreshold) continue;
            double monoIsotopeMass;
            double isotopeCosine = getIsotopeCosine(monoIsotopeMass, pg, perIsotopeIntensities, generator, param);
            if (isotopeCosine <= param.isotopeCosineThreshold) continue;

            double totalIntensity = accumulate(perChargeIntensities, perChargeIntensities+param.chargeRange,.0);
            intensityMassMap[-totalIntensity] = monoIsotopeMass;
        }
        return intensityMassMap;
    }

    void onlineDeconvolution(MSExperiment &map, Parameter &param, fstream &fs, int &specCntr, int &qspecCntr, int &massCntr) {

        double filter[param.chargeRange];
        double harmonicFilter[param.chargeRange];

        auto generator = new CoarseIsotopePatternGenerator();
        auto maxIso = generator->estimateFromPeptideWeight(param.maxMass);
        maxIso.trimRight(.1*maxIso.getMostAbundant().getIntensity());
        param.maxIsotopeIndex = min(param.maxIsotopeIndex, (int)maxIso.size());
        generator->setMaxIsotope(param.maxIsotopeIndex);

        for (int i = 0; i < param.chargeRange; i++) {
            filter[i] = log(1.0 / (i + param.minCharge)); // should be descending, and negative!
            harmonicFilter[i] = log(1.0 / (i + .5 + param.minCharge));
        }

        for (auto it = map.begin(); it != map.end(); ++it) {
            if (it->getMSLevel() != 1) continue;
            specCntr++;
            double rt = it->getRT();

            auto logMzPeaks = GetLogMzPeaks(*it, param);
            auto peakGroups = findPeakGroups(logMzPeaks, filter, harmonicFilter, param);

            auto intensityMassMap = getIntensityMassMap(peakGroups, generator, param);
            if (intensityMassMap.empty()) continue;

            qspecCntr++;

            int massCntrForThis = 0;
            for (auto iter = intensityMassMap.begin(); iter != intensityMassMap.end(); ++iter) {
                auto intensity = -iter->first;
                auto m = iter->second;
                fs << massCntrForThis << " " << qspecCntr << " " << specCntr << " " << it->getRT() << " "  <<intensity<<" "   << setprecision(15)
                   << m << " " << round(m * 0.999497) << endl;
                massCntrForThis++;
                if (massCntrForThis > param.maxMassCount) break;
            }
            massCntr += massCntrForThis;

        }
        return;
    }

    double getCosine(double *a, IsotopeDistribution b, int offset = 0) {
        double n = 0, d1 = 0, d2 = 0;
        //auto bMax = b.getMostAbundant();
        for (int i = 0; i < b.size(); i++) {
            int j = i + offset;
            if (j < 0 || j >= b.size()) continue;
            double bInt = b[i].getIntensity();
            //if (b[i].getMZ() > bMax.getMZ() && bInt / bMax.getIntensity() < rightThreshold) break;
            d2 += bInt * bInt;

            n += a[j] * bInt;
            d1 += a[j] * a[j];

        }
        double d = sqrt(d1 * d2);
        if (d <= 0) return 0;
        return n / d;
    }

    boost::dynamic_bitset<> getMzBins(vector<LogMzPeak> &logMzPeaks, double tol) {
        int binNumber = (logMzPeaks[logMzPeaks.size() - 1].logMz  - logMzPeaks[0].logMz) / tol + 2;
        boost::dynamic_bitset<> mzBins(binNumber);

        for (auto &p : logMzPeaks) {
            int bi = (int) ((p.logMz - logMzPeaks[0].logMz) / tol);
            if (bi >= binNumber) break;
            mzBins[bi] = true;
        }
        mzBins |= mzBins << 1;
        //mzBins |= mzBins >> 1;
        return mzBins;
    }

    void getMassBins(boost::dynamic_bitset<> &massBins, boost::dynamic_bitset<> &mzBins,
            Byte *massBinScores, double minBinLogMass, double minLogMz,
                     double *filter, double *harmonicFilter,  int binNumber, const Parameter& param) {
        fill_n(massBinScores, binNumber, 0);

        int minChargeCount = param.minChargeCount;
        int minContinuousCharge = param.minContinuousCharge;

        int binOffsets[param.chargeRange];
        int harmonicBinOffsets[param.chargeRange];

        for (int i = 0; i < param.chargeRange; i++) {
            binOffsets[i] = (minLogMz-filter[i]-minBinLogMass) / param.tolerance;
            harmonicBinOffsets[i] = (minLogMz-harmonicFilter[i]-minBinLogMass) / param.tolerance;
        }

        unsigned long setBinIndex = mzBins.find_first();

        while (setBinIndex != mzBins.npos) {
            for (int j = 0; j < param.chargeRange; j++) {
                int bi = setBinIndex + binOffsets[j];
                if (bi<0) continue;
                if (bi >= binNumber) break;
                massBins[bi] = ++massBinScores[bi] >= minChargeCount;
            }
            setBinIndex = mzBins.find_next(setBinIndex);
        }

        fill_n(massBinScores, binNumber, 0); // re-init
        setBinIndex = massBins.find_first();
        double threshold = (log(param.minMass) - minBinLogMass)/param.tolerance;

        while (setBinIndex != massBins.npos) {
            if(setBinIndex <threshold){
                massBins[setBinIndex] = false;
                setBinIndex = massBins.find_next(setBinIndex);
                continue;
            }

            Byte continuousChargeCntr = 0;
            bool take = false;
            for (int j = 0; j < param.chargeRange; j++) {
                int bi = setBinIndex - binOffsets[j];
                if (bi < 0) break;
                int hbi = setBinIndex - harmonicBinOffsets[j];
                if (mzBins[bi] && (hbi < 0 || !mzBins[hbi])) {
                    take = ++continuousChargeCntr >= minContinuousCharge;
                    if(take){
                        Byte score = massBinScores[setBinIndex];
                        massBinScores[setBinIndex] = score < continuousChargeCntr ? continuousChargeCntr  : score;
                    }
                } else continuousChargeCntr = 0;
            }
            massBins[setBinIndex] = take;
            setBinIndex = massBins.find_next(setBinIndex);
        }
    }

    double calculateBinScoreThreshold(boost::dynamic_bitset<> &massBins, Byte *massBinScores, const Parameter &param){
        int *binScoreDist = new int[param.chargeRange + 1];
        fill_n(binScoreDist, param.chargeRange + 1, 0);

        auto setBinIndex = massBins.find_first();
        while (setBinIndex != massBins.npos) {
            binScoreDist[massBinScores[setBinIndex]]++;
            setBinIndex = massBins.find_next(setBinIndex);
        }

        int binScoreThreshold = param.minContinuousCharge;
        int sum = 0;
        for (int i = param.chargeRange; i >= 0; i--) {
            sum += binScoreDist[i];
            if (sum >= param.maxMassCount * 2) {
                binScoreThreshold = i > binScoreThreshold ? i : binScoreThreshold;
                break;
            }
        }
        return binScoreThreshold;
    }


    void getPeakGroups(vector<PeakGroup> &peakGroups, boost::dynamic_bitset<> &massBins, Byte *massBinScores,
                          int binScoreThreshold, vector<LogMzPeak> &logMzPeaks, double minBinLogMass, double *filter, const Parameter& param) {

        int *currentPeakIndex = new int[param.chargeRange];
        fill_n(currentPeakIndex, param.chargeRange, 0);

        auto setBinIndex = massBins.find_first();
        while (setBinIndex != massBins.npos) {
            if (massBinScores[setBinIndex] >= binScoreThreshold) {

                int isoOff = 0, maxi = 0;
                PeakGroup peakGroup;
                vector<int> peakMassBins;
                peakMassBins.reserve(logMzPeaks.size());
                peakGroup.reserve(logMzPeaks.size());

                for (int j = 0; j < param.chargeRange; j++) {
                    int charge = j + param.minCharge;
                    while (currentPeakIndex[j] < logMzPeaks.size()) {
                        double logMz = logMzPeaks[currentPeakIndex[j]].logMz;
                        double mz = logMzPeaks[currentPeakIndex[j]].orgPeak->getMZ() - Constants::PROTON_MASS_U;

                        int bi = (logMz - filter[j] - minBinLogMass) / param.tolerance;
                        if (bi > setBinIndex) break;
                        if (setBinIndex == bi) {
                            LogMzPeak p(*logMzPeaks[currentPeakIndex[j]].orgPeak, charge, 0);
                            peakGroup.push_back(p);
                            peakMassBins.push_back(bi);

                            if (currentPeakIndex[j] > 0 ) {
                                if (logMz - logMzPeaks[currentPeakIndex[j] - 1].logMz < param.tolerance) {
                                    LogMzPeak p(*logMzPeaks[currentPeakIndex[j] - 1].orgPeak, charge, 0);
                                    //peakMassBins.push_back((p.logMz - filter[j] - minLogMass) / tol);
                                    peakGroup.push_back(p);
                                }
                            }

                            if (currentPeakIndex[j] < logMzPeaks.size()-1) {
                                if (logMzPeaks[currentPeakIndex[j] + 1].logMz - logMz < param.tolerance) {
                                    LogMzPeak p(*logMzPeaks[currentPeakIndex[j] + 1].orgPeak, charge, 0);
                                    peakMassBins.push_back((p.logMz - filter[j] - minBinLogMass) / param.tolerance);
                                    peakGroup.push_back(p);
                                }
                            }
                            //continue;
                            for (int d = -1; d <= 1; d += 2) { // negative then positive direction.
                                int currentPeakIndexForIsotopes = currentPeakIndex[j] + d;
                                for (int i = 1; currentPeakIndexForIsotopes>=0 && currentPeakIndexForIsotopes < logMzPeaks.size() ; i++) {
                                    double logMzIsotope = logMz + Constants::C13C12_MASSDIFF_U * i * d / charge / mz; //log(mz + Constants::C13C12_MASSDIFF_U * i * d / charge);

                                    bool isotopePeakPresent = false;
                                    while (currentPeakIndexForIsotopes >= 0 &&
                                           currentPeakIndexForIsotopes < logMzPeaks.size()) {
                                        double logMzForIsotopePeak = logMzPeaks[currentPeakIndexForIsotopes].logMz;
                                        currentPeakIndexForIsotopes += d;
                                        if (logMzForIsotopePeak < logMzIsotope - param.tolerance) if (d < 0) break; else continue;
                                        if (logMzForIsotopePeak > logMzIsotope + param.tolerance) if (d < 0) continue; else break;

                                        isotopePeakPresent = true;

                                        LogMzPeak p(*logMzPeaks[currentPeakIndexForIsotopes - d].orgPeak, charge,
                                                    i * d);
                                        if(d>0){
                                            peakMassBins.push_back((p.logMz - filter[j]- minBinLogMass) / param.tolerance);
                                        }
                                       // if(i - isoOff < param.maxIsotopeIndex)
                                        peakGroup.push_back(p);
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
                        currentPeakIndex[j]++;
                    }
                }

                if(maxi + isoOff + 1 >= param.minIsotopeCount){
                    for(auto& i :peakMassBins){
                        massBins[i] = false;
                    }
                    for (LogMzPeak &p : peakGroup.peaks) {
                        p.isotopeIndex = p.isotopeIndex - isoOff;
                    }
                    peakGroups.push_back(peakGroup);
                }
            }
            setBinIndex = massBins.find_next(setBinIndex);
            //cout<<setBinIndex<<" "<<massBins.count()<<endl;
        }
    }

    vector<PeakGroup> findPeakGroups(vector<LogMzPeak> &logMzPeaks, double *filter, double *hFilter, const Parameter &param){
        vector<PeakGroup> peakGroups;
        double maxBinLogMass = std::min(logMzPeaks[logMzPeaks.size() - 1].logMz - filter[param.chargeRange - param.minChargeCount], log(param.maxMass));
        double minBinLogMass = logMzPeaks[0].logMz - filter[param.minChargeCount];

        peakGroups.reserve(param.maxMassCount * 100);
        int massBinNumber = (maxBinLogMass - minBinLogMass) / param.tolerance + 1;
        boost::dynamic_bitset<> mzBins = getMzBins(logMzPeaks, param.tolerance);
        boost::dynamic_bitset<> massBins(massBinNumber);
        Byte *massBinScores = new Byte[massBinNumber];
        getMassBins(massBins, mzBins, massBinScores, minBinLogMass, logMzPeaks[0].logMz, filter, hFilter,  massBinNumber, param);
        double binScoreThreshold = calculateBinScoreThreshold(massBins, massBinScores, param);
        getPeakGroups(peakGroups, massBins, massBinScores, binScoreThreshold, logMzPeaks, minBinLogMass, filter, param);
        return peakGroups;
    }
};

    int main(int argc, const char **argv) {
        FlashDeconv tool;
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