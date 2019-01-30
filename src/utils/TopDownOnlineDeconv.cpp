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

class TopDownRealTimeDeconvolution:
        public TOPPBase {

public:
    TopDownRealTimeDeconvolution() :
            TOPPBase("TopDownRealTimeDeconvolution", "Real-time Deconvolution for Non-redundant MS2 acquisition with top down data",
                     false)\
 {}

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
                onlineDeconvolution(map, fs, specCntr, qspecCntr, massCntr);
                clock_t end = clock();
                elapsed_secs += double(end - begin) / CLOCKS_PER_SEC;
                std::cout << "%" << qspecCntr << " MS1 spectra deconvoluted so far" << endl;
            }
            fs << "];";
            fs.close();
        }

        std::cout << "%" << elapsed_secs << " seconds elapsed for " << specCntr << " MS1 spectra" << endl;
        std::cout << "%" << elapsed_secs / specCntr * 1000 << " msec per spectrum" << std::endl;

        return EXECUTION_OK;
    }

    void onlineDeconvolution(MSExperiment &map, fstream &fs, int &specCntr, int &qspecCntr, int &massCntr) {
        double intensityThreshold = 1000;// the lower the more diverse molecules
        int filterSize = 25;
        int minCharge = 5;
        int minContinuousChargePeak = 3;//
        int minChargePeak = filterSize * .2;
        int minIsotopeCount = 4; //
        int maxMassPerSpectrum = 30;//
        int maxIsotopeIndex = 15; // make it variable
        double minMass = 500.0;
        double maxMass = 50000.0;

        double tolerance = 5e-6; // 5 ppm
        double isotopeCosineThreshold = .5;

        double filter[filterSize];
        double harmonicFilter[filterSize];

        auto generator = new CoarseIsotopePatternGenerator();
        generator->setMaxIsotope(maxIsotopeIndex);

        for (int i = 0; i < filterSize; i++) {
            filter[i] = log(1.0 / (i + minCharge)); // should be descending, and negative!
            harmonicFilter[i] = log(1.0 / (i + .5 + minCharge));
        }

        for (auto it = map.begin(); it != map.end(); ++it) {
            if (it->getMSLevel() != 1) continue;
            specCntr++;

            double rt = it->getRT();

            vector<LogMzPeak> logMzPeaks;
            logMzPeaks.reserve(it->size());
            for (auto &peak : *it) {
                if (peak.getIntensity() <= intensityThreshold) continue;
                LogMzPeak logMzPeak(peak);
                logMzPeaks.push_back(logMzPeak);
            }

            auto peakGroups =
                    findPeakGroups(logMzPeaks, filter, harmonicFilter, filterSize, minCharge, maxIsotopeIndex, minMass, maxMass,
                                   tolerance, minContinuousChargePeak, minChargePeak, 1,
                                   maxMassPerSpectrum * 2);
            if (peakGroups.empty())
                continue;

            bool qualified = false;
            std::map<double, double> intensityMassMap; // intensity -> mass

            for (auto &pg : peakGroups) {
                auto perChargeIntensities = new double[filterSize];
                auto perIsotopeIntensities = new double[maxIsotopeIndex + 1];

                fill_n(perChargeIntensities, filterSize, 0);
                fill_n(perIsotopeIntensities, maxIsotopeIndex + 1, 0); // change it to be flexible

                double monoIsotopeMass;
                double maxIntensityForMonoIsotopeMass = -1;
                double mostAbundantMass;
                double maxIntensityForMostAbundantMass = -1;

                for (auto &p : pg) {
                    if (p.isotopeIndex > maxIsotopeIndex) continue;
                    int index = p.charge - minCharge;
                    double intensity = p.orgPeak->getIntensity();
                    perChargeIntensities[index] += intensity;
                    perIsotopeIntensities[p.isotopeIndex] += intensity;

                    if (maxIntensityForMostAbundantMass <= intensity) {
                        maxIntensityForMostAbundantMass = intensity;
                        mostAbundantMass = p.getMass();
                    }


                    if (maxIntensityForMonoIsotopeMass > intensity) continue;
                    maxIntensityForMonoIsotopeMass = intensity;
                    monoIsotopeMass = p.getMonoIsotopeMass();

                }

                int setIntensityCounter = 0;
                int maxSetIntensityCounter = 0;
                for (int i = 0; i < filterSize; i++) {
                    if (perChargeIntensities[i] <= 0) {
                        setIntensityCounter = 0;
                        continue;
                    }
                    setIntensityCounter++;
                    maxSetIntensityCounter =
                            maxSetIntensityCounter > setIntensityCounter ? maxSetIntensityCounter : setIntensityCounter;
                }
                if (maxSetIntensityCounter < minContinuousChargePeak) continue;

                setIntensityCounter = 0;
                maxSetIntensityCounter = 0;
                for (int i = 0; i < maxIsotopeIndex; i++) {
                    if (perIsotopeIntensities[i] <= 0) {
                        setIntensityCounter = 0;
                        continue;
                    }
                    setIntensityCounter++;
                    maxSetIntensityCounter =
                            maxSetIntensityCounter > setIntensityCounter ? maxSetIntensityCounter : setIntensityCounter;
                }
                if (maxSetIntensityCounter < minIsotopeCount) continue;

                int maxIntensityIndex;
                double maxIntensity = -1;
                double totalIntensity = 0;
                for (int i = 0; i < filterSize; i++) {
                    totalIntensity += perChargeIntensities[i];
                    if (maxIntensity < perChargeIntensities[i]) {
                        maxIntensity = perChargeIntensities[i];
                        maxIntensityIndex = i;
                    }
                }

                double score = 0;
                for (int k = 1; k < filterSize; k++) {
                    int d1 = k<=maxIntensityIndex? 0 : -1;
                    int d2 = k<=maxIntensityIndex? -1 : 0;
                    double int1 = perChargeIntensities[k + d1];
                    double int2 = perChargeIntensities[k + d2];
                    if (int2 <= 0) continue;
                    score += int1 >= int2 ? 1 : -1;
                }

                if(score<=0) continue;


                auto avgIso = generator->estimateFromPeptideWeight(mostAbundantMass);

                int offset = 0;
                double maxCosine = -1;
                for (int f = -3; f <= 3; f++) {
                    auto cos = getCosine(perIsotopeIntensities, avgIso, .01, f);

                    if (maxCosine < cos) {
                        maxCosine = cos;
                        offset = f;
                    }
                }

                if (maxCosine < isotopeCosineThreshold) {
                    continue;
                }

                monoIsotopeMass += offset * Constants::C13C12_MASSDIFF_U;

                intensityMassMap[-totalIntensity] = monoIsotopeMass;
                qualified = true;


                /*if(monoIsotopeMass<1.81e4&&monoIsotopeMass>1.78e4){
                    for (int i = 0; i < maxIsotopeIndex; i++) {
                        cout<<perIsotopeIntensities[i]<<",";
                    }
                    cout<<endl;
                    for (int i = 0; i < filterSize; i++) {
                        cout<<perChargeIntensities[i]<<",";
                    }
                    cout<<endl<<endl;
                }*/


            }
            if (!qualified) continue;

            qspecCntr++;

            int massCntr = 0;
            for (auto iter = intensityMassMap.begin(); iter != intensityMassMap.end(); ++iter) {
                auto intensity = -iter->first;
                auto m = iter->second;
                fs << massCntr++ << " " << qspecCntr << " " << specCntr << " " << it->getRT() << " "  <<intensity<<" "   << setprecision(15)
                   << m << " " << round(m * 0.999497) << endl;
                massCntr++;
                if (massCntr > maxMassPerSpectrum) break;
                //cout<<massCntr++<<" "<<qspecCntr << " " << specCntr << " " << it->getRT() <<" " << setprecision(15)<< m << " " << round(m*0.999497) << " " << intensity << endl;
            }
            //cout<<endl;
        }
        return;
    }


    double getCosine(double *a, IsotopeDistribution b, double rightThreshold, int offset = 0) {
        double n = 0, d1 = 0, d2 = 0;
        auto bMax = b.getMostAbundant();
        for (int i = 0; i < b.size(); i++) {
            int j = i + offset;
            if (j < 0 || j >= b.size()) continue;
            double bInt = b[i].getIntensity();
            if (b[i].getMZ() > bMax.getMZ() && bInt / bMax.getIntensity() < rightThreshold) break;
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
        return mzBins;
    }

    void getMassBins(boost::dynamic_bitset<> &massBins, boost::dynamic_bitset<> &mzBins, Byte *massBinScores,
                     double minLogMass, double minLogMz, double minMass,
                     double *filter, double *harmonicFilter, int filterSize, int binNumber, double tol,
                     int minContinuousChargePeak, int minChargePeak) {
        fill_n(massBinScores, binNumber, 0);

        int binOffsets[filterSize];
        int harmonicBinOffsets[filterSize];

        for (int i = 0; i < filterSize; i++) {
            binOffsets[i] = (minLogMz-filter[i]-minLogMass) / tol;
            harmonicBinOffsets[i] = (minLogMz-harmonicFilter[i]-minLogMass) / tol;
        }

        unsigned long setBinIndex = mzBins.find_first();

        while (setBinIndex != mzBins.npos) {
            for (int j = 0; j < filterSize; j++) {
                int bi = setBinIndex + binOffsets[j];
                if (bi<0) continue;
                if (bi >= binNumber) break;
                massBinScores[bi]++;
                if (massBinScores[bi] >= minChargePeak) {
                    massBins[bi] = true;
                }
            }
            setBinIndex = mzBins.find_next(setBinIndex);
        }

        fill_n(massBinScores, binNumber, 0);
        setBinIndex = massBins.find_first();
        double logMinMass = log(minMass);

        while (setBinIndex != massBins.npos) {
            //auto m = exp(minLogMass + setBinIndex * tol);
            if(minLogMass + setBinIndex * tol<logMinMass){
                massBins[setBinIndex] = false;
                setBinIndex = massBins.find_next(setBinIndex);
                continue;
            }

            Byte continuousChargeCntr = 0;
            bool take = false;
            for (int j = 0; j < filterSize; j++) {
                int bi = setBinIndex - binOffsets[j];
                if (bi < 0) break;
                int hbi = setBinIndex - harmonicBinOffsets[j];
                if (mzBins[bi] && (hbi < 0 || !mzBins[hbi])) {
                    continuousChargeCntr++;
                    take = continuousChargeCntr >= minContinuousChargePeak;
                    massBinScores[setBinIndex] =
                            massBinScores[setBinIndex] < continuousChargeCntr ? continuousChargeCntr
                                                                              : massBinScores[setBinIndex];
                } else continuousChargeCntr = 0;
            }
            massBins[setBinIndex] = take;
            setBinIndex = massBins.find_next(setBinIndex);
        }

    }

    /*
    void filterMassBinsUsingIsotopeCriteria(boost::dynamic_bitset<> &massBins, Byte *massBinScores, int maxIsotopeIndex,
                                            int minIsotopeCount, int binNumber, double tol, double minLogMass, double minMass) {

        unsigned long setBinIndex = massBins.find_first();

        while (setBinIndex != massBins.npos) {
            //auto allIsoSet = true;
            auto m = exp(minLogMass + setBinIndex * tol);
            if(m<minMass){
                massBins[setBinIndex] = false;
                setBinIndex = massBins.find_next(setBinIndex);
                continue;
            }

            int *isoBins = new int[maxIsotopeIndex];
            fill_n(isoBins, maxIsotopeIndex, 0);

            for (int i = 1; i <= maxIsotopeIndex; i++) {
                int bin = setBinIndex + Constants::C13C12_MASSDIFF_U * i / m / tol;
                if (bin >= binNumber - 1) {
                    allIsoSet = false;
                    break;
                }
                bool isoSet = massBins[bin] | massBins[bin + 1];
                if (i < minIsotopeCount) allIsoSet &= isoSet;//
                if (!isoSet) break;
                isoBins[i - 1] = bin;
            }

            if (allIsoSet) {
                for (int i = 1; i <= maxIsotopeIndex; i++) {
                    int bin = isoBins[i - 1];
                    if (bin == 0) break;
                    auto score = massBinScores[setBinIndex];
                    score = score < massBinScores[bin] ? massBinScores[bin] : score;
                    score = score < massBinScores[bin + 1] ? massBinScores[bin] : score;
                    massBinScores[setBinIndex] = score;
                    massBins[bin] = false;
                    if (bin < binNumber - 1) massBins[bin + 1] = false;
                }
                if (setBinIndex > 0) massBins[setBinIndex - 1] = false;
            } else {
                massBins[setBinIndex] = false;
            }
            setBinIndex = massBins.find_next(setBinIndex);
        }
    }*/

    double calculateBinScoreThreshold(boost::dynamic_bitset<> &massBins, Byte *massBinScores,
                                      int minChargePeakCount, int maxMassCountPerSpectrum, int filterSize) {
        int *binScoreDist = new int[filterSize + 1];
        fill_n(binScoreDist, filterSize + 1, 0);

        auto setBinIndex = massBins.find_first();
        while (setBinIndex != massBins.npos) {
            binScoreDist[massBinScores[setBinIndex]]++;
            setBinIndex = massBins.find_next(setBinIndex);
        }

        int binScoreThreshold = minChargePeakCount, tsum = 0;
        for (int i = filterSize; i >= 0; i--) {
            tsum += binScoreDist[i];
            if (tsum >= maxMassCountPerSpectrum) {
                binScoreThreshold = i > binScoreThreshold ? i : binScoreThreshold;
                break;
            }
        }
        return binScoreThreshold;
    }


    void updatePeakGroups(vector<vector<LogMzPeak>> &peakGroups, boost::dynamic_bitset<> &massBins, Byte *massBinScores,
                          int binScoreThreshold, vector<LogMzPeak> &logMzPeaks, int filterSize, int minCharge,
                          int minIsotopeCount, int maxIsotopeIndex, double minLogMass, double *filter, double tol) {

        int *currentPeakIndex = new int[filterSize];
        fill_n(currentPeakIndex, filterSize, 0);

        auto setBinIndex = massBins.find_first();
        while (setBinIndex != massBins.npos) {
            if (massBinScores[setBinIndex] >= binScoreThreshold) {

                int isoOff = 0, maxi = 0;
                vector<LogMzPeak> peakGroup;
                vector<int> peakMassBins;
                peakMassBins.reserve(logMzPeaks.size());
                peakGroup.reserve(logMzPeaks.size());

                for (int j = 0; j < filterSize; j++) {
                    int charge = j + minCharge;
                    while (currentPeakIndex[j] < logMzPeaks.size()) {
                        double logMz = logMzPeaks[currentPeakIndex[j]].logMz;
                        double mz = logMzPeaks[currentPeakIndex[j]].orgPeak->getMZ() - Constants::PROTON_MASS_U;

                        int bi = (logMz - filter[j] - minLogMass) / tol;
                        if (bi > setBinIndex) break;
                        if (setBinIndex == bi) {
                            LogMzPeak p(*logMzPeaks[currentPeakIndex[j]].orgPeak, charge, 0);
                            peakGroup.push_back(p);
                            peakMassBins.push_back(bi);

                            if (currentPeakIndex[j] > 0 ) {
                                if (logMz - logMzPeaks[currentPeakIndex[j] - 1].logMz < tol) {
                                    LogMzPeak p(*logMzPeaks[currentPeakIndex[j] - 1].orgPeak, charge, 0);
                                    //peakMassBins.push_back((p.logMz - filter[j] - minLogMass) / tol);
                                    peakGroup.push_back(p);
                                }
                            }

                            if (currentPeakIndex[j] < logMzPeaks.size()-1) {
                                if (logMzPeaks[currentPeakIndex[j] + 1].logMz - logMz < tol) {
                                    LogMzPeak p(*logMzPeaks[currentPeakIndex[j] + 1].orgPeak, charge, 0);
                                    peakMassBins.push_back((p.logMz - filter[j] - minLogMass) / tol);
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
                                        if (logMzForIsotopePeak < logMzIsotope - tol) if (d < 0) break; else continue;
                                        if (logMzForIsotopePeak > logMzIsotope + tol) if (d < 0) continue; else break;

                                        isotopePeakPresent = true;

                                        LogMzPeak p(*logMzPeaks[currentPeakIndexForIsotopes - d].orgPeak, charge,
                                                    i * d);
                                        if(d>0){
                                            peakMassBins.push_back((p.logMz - filter[j]- minLogMass) / tol);
                                        }
                                        if(i - isoOff < maxIsotopeIndex) peakGroup.push_back(p);
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

                //cout<<isoOff<<" "<<maxi<<endl;
                if(maxi + isoOff + 1 <minIsotopeCount) continue;

                for(auto& i :peakMassBins){
                    massBins[i] = false;
                }

                for (LogMzPeak &p : peakGroup) {
                    p.isotopeIndex = p.isotopeIndex - isoOff;
                }
                /*if(isoOff==-15){
                    for (LogMzPeak &p : peakGroup) {
                        cout<<p.getMass()<<" "<<p.orgPeak->getMZ()<<" "<<p.charge<<" "<<p.isotopeIndex<<endl;
                    }
                    cout<<endl;
                }*/
                peakGroups.push_back(peakGroup);
            }
            setBinIndex = massBins.find_next(setBinIndex);
        }
    }

    vector<vector<LogMzPeak>>
    findPeakGroups(vector<LogMzPeak> &logMzPeaks, double *filter, double *harmonicFilter, int filterSize, int minCharge,
                   int maxIsotopeIndex, double minMass, double maxMass, double tol, int minContinuousChargePeak, int minChargePeak,
                   int minIsotopeCount, int maxMassCountPerSpectrum) {
        vector<vector<LogMzPeak>> peakGroups;
        // minContinuousChargePeak, minChargePeak should be < filterSize
        double maxLogMass = std::min(logMzPeaks[logMzPeaks.size() - 1].logMz - filter[filterSize - minChargePeak], log(maxMass));
        double minLogMass = logMzPeaks[0].logMz - filter[minChargePeak];

        peakGroups.reserve(maxMassCountPerSpectrum * 10);
        //cout<<0;

        int massBinNumber = (maxLogMass - minLogMass) / tol + 1;
        boost::dynamic_bitset<> mzBins = getMzBins(logMzPeaks, tol);
        boost::dynamic_bitset<> massBins(massBinNumber);
        Byte *massBinScores = new Byte[massBinNumber];

        //cout<<1;

        getMassBins(massBins, mzBins, massBinScores,
                    minLogMass, logMzPeaks[0].logMz, minMass,
                    filter, harmonicFilter, filterSize, massBinNumber, tol,
                    minContinuousChargePeak, minChargePeak);
        //filterMassBinsUsingIsotopeCriteria(massBins, massBinScores, maxIsotopeIndex, minIsotopeCount, massBinNumber, tol,
         //                                  minLogMass, minMass);
        //cout<<2;

        double binScoreThreshold = calculateBinScoreThreshold(massBins, massBinScores, minContinuousChargePeak,
                                                              maxMassCountPerSpectrum, filterSize);
        //cout<<3;

        updatePeakGroups(peakGroups, massBins, massBinScores, binScoreThreshold, logMzPeaks, filterSize, minCharge,
                         minIsotopeCount, maxIsotopeIndex, minLogMass, filter, tol);
        //cout<<4<<endl;

        return peakGroups;
    }
};

    int main(int argc, const char **argv) {
        TopDownRealTimeDeconvolution tool;
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