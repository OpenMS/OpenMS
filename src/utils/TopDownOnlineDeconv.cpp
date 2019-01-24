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

using namespace OpenMS;
using namespace std;

class TDOnlineDeconv:
        public TOPPBase
{

public:
    TDOnlineDeconv():
            TOPPBase("TDOnlineDeconv", "Online Deconvolution for Smart MS2 acquisition with top down data", false)\
    {}

    typedef struct LogMzPeak{
        Peak1D *orgPeak;
        double logMz;
        int charge;
        int isotopeIndex;
        double score;

        LogMzPeak(): orgPeak(nullptr), logMz(-10000), charge(0), isotopeIndex(0), score(0) {}
        LogMzPeak(Peak1D &peak): orgPeak(&peak), logMz(getLogMz(peak.getMZ())), charge(0), isotopeIndex(0), score(0) {}

        double getMass(){
            return exp(logMz) * charge;
        }
        bool operator<(const LogMzPeak &a){
            return logMz < a.logMz;
        }
    };


protected:

    static double getLogMz(double mz){
        return log(mz -  Constants::PROTON_MASS_U);
    }


    void registerOptionsAndFlags_() override {
        registerInputFile_("in", "<file>", "", "Input file.");
        setValidFormats_("in", ListUtils::create<String>("mzML"));
    }
    ExitCodes main_(int, const char **) override
    {
        //-------------------------------------------------------------
        // parsing parameters
        //-------------------------------------------------------------
//        String infilePath = getStringOption_("in");

        String infilePath = "/Users/kyowonjeong/Documents/A4B/mzml/MS1only/05-26-17_B7A_yeast_td_fract12_rep2_MS1only.mzML";
        //infilePath = "/Users/kyowonjeong/Documents/A4B/mzml/180523_Cytocrome_C_MS2_HCD.mzML";
        //infilePath = "/Users/kyowonjeong/Documents/A4B/mzml/180523_Myoglobin_MS2_HCD.mzML";
        cout << "file name : " << infilePath << endl;
        // just for quick use

        //-------------------------------------------------------------
        // reading input
        //-------------------------------------------------------------

        MzMLFile mzml;
        mzml.setLogType(log_type_);

        // load input
        MSExperiment map;
        mzml.load(infilePath, map);
        cout << "Loaded consensus maps" << endl;
        clock_t begin = clock();
        int cntr = onlineDeconvolution(map);
        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        std::cout << elapsed_secs << " seconds elapsed for " <<cntr << " MS1 spectra" <<  endl;
        std::cout << elapsed_secs/cntr*1000 << " msec per spectrum"<< std::endl;

        return EXECUTION_OK;
    }

    int onlineDeconvolution(MSExperiment &map){
        double threshold = 5000;//
        int filterSize = 25;//2069 601 3006 21563.2
        int minCharge = 5;
        int minPeakPerMass = 6; // 180 124
        int minIsotopeCount = 3;
        int maxMassPerSpectrum = 6;
        int maxIsotopeIndex = 15; // make it variable TODO
        double tolerance = 5e-6; // 5 ppm

        int specCntr = 0, qspecCntr = 0, massCntr=0;
        double filter[filterSize];

        CoarseIsotopePatternGenerator* generator= new CoarseIsotopePatternGenerator();
        generator->setMaxIsotope(maxIsotopeIndex);

        for(int i=0;i<filterSize;i++){
            filter[i] = log(1.0/(i + minCharge)); // should be descending, and negative!
        }

        for (auto it = map.begin(); it != map.end(); ++it) {
            if (it->getMSLevel() != 1) continue;
            specCntr++;

            double rt = it->getRT();

            vector<LogMzPeak> logMzPeaks;
            logMzPeaks.reserve(it->size());
            for (auto &peak : *it){
                if(peak.getIntensity() <= threshold) continue;
                LogMzPeak logMzPeak(peak);
                logMzPeaks.push_back(logMzPeak);
            }

            auto peakGroups =
                    findPeakGroups(logMzPeaks, filter, filterSize, minCharge, maxIsotopeIndex,
                            tolerance, minPeakPerMass, minIsotopeCount, maxMassPerSpectrum);
            if (peakGroups.empty())
                continue;
            std::map<double, double> scoreMassMap;

            bool qualified = false;
            for(auto &pg : peakGroups){
                double* perChargeIntensities = new double[filterSize];
                double* perIsotopeIntensities = new double[maxIsotopeIndex];

                fill_n(perChargeIntensities,filterSize,0);
                fill_n(perIsotopeIntensities,maxIsotopeIndex,0); // change it to be flexible

                double monoIsotopeMass;
                for(auto &p : pg){
                    int index = p.charge - minCharge;
                    perChargeIntensities[index] += p.orgPeak->getIntensity();
                    if(p.isotopeIndex<maxIsotopeIndex)
                        perIsotopeIntensities[p.isotopeIndex] += p.orgPeak->getIntensity();
                    if(p.isotopeIndex == 0) monoIsotopeMass = p.getMass();
                }

                auto avgIso = generator->estimateFromPeptideWeight(monoIsotopeMass);
                avgIso.renormalize();
                avgIso.trimRight(0.01);//
                int offset;
                double maxCosine = -1;
                for(int f = -3;f<=3;f++){
                    auto cos = getCosine(perIsotopeIntensities, avgIso, f);

                    if(maxCosine < cos){
                        maxCosine = cos;
                        offset = f;
                    }
                }

                if(maxCosine < .5){
                   continue;
                }

                //cout<<offset<< " " << maxCosine << endl;
                monoIsotopeMass += offset * Constants::C13C12_MASSDIFF_U;

                //
                int maxIntensityIndex;
                int firstNonZeroIndex = -1;
                int lastNonZeroIndex = filterSize;
                double maxIntensity = -1;

                for(int i=0;i<filterSize;i++) {
                    if(perChargeIntensities[i]!=0){
                        lastNonZeroIndex = i;
                        if(firstNonZeroIndex<0 ) firstNonZeroIndex=i;
                    }

                    if (maxIntensity < perChargeIntensities[i]) {
                        maxIntensity = perChargeIntensities[i];
                        maxIntensityIndex = i;
                    }
                }

                double score = 0;
                for(int k=maxIntensityIndex;k<lastNonZeroIndex;k++){
                    if(k>=filterSize-1) break;
                    score += getLogLikelihoodRatioScore(perChargeIntensities[k], perChargeIntensities[k+1]);
                }
                for(int k=maxIntensityIndex;k>firstNonZeroIndex;k--){
                    if(k<=0) break;
                    score += getLogLikelihoodRatioScore(perChargeIntensities[k], perChargeIntensities[k-1]);
                }


                //continue;
                if(score > 0 && !pg.empty()){
                    cout<<massCntr++<<" "<<qspecCntr << " " << specCntr << " " << monoIsotopeMass << endl;
                    qualified = true;
                }
            }
            if(!qualified) continue;

           qspecCntr++;

           // for (int i=0; i<20 && i<scoreMassMap.size() ; ++iter, i++){
           //     cout << qspecCntr<<" "<<specCntr<<" "<<iter->first << " " << iter->second << " " << peakGroups.size() << " " << scoreMassMap.size()<< endl;
           // }
           // cout << endl;
        }
        return specCntr;
    }

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

    double getCosine(double* a, IsotopeDistribution b, int offset = 0){
        double n=0, d1=0, d2=0;
        for(int i=0;i<b.size();i++){
            int j = i + offset;
            if(j<0 || j>=b.size()) continue;
            n += a[j] * b[i].getIntensity();
            d1 += a[j] * a[j];
            d2 += b[i].getIntensity() * b[i].getIntensity();
        }
        double d = sqrt(d1 * d2);
        if(d <=0) return 0;
        return n/d;
    }

    boost::dynamic_bitset<> getMzBins(vector<LogMzPeak> &logMzPeaks, int binNumber, double tol){
        boost::dynamic_bitset<> mzBins(binNumber);

        for (auto &p : logMzPeaks) {
            int bi = (int)((p.logMz - logMzPeaks[0].logMz) / tol);
            if (bi >= binNumber) break;
            mzBins[bi] = true;
        }
        mzBins |= mzBins<<1;
        return mzBins;
    }

    void getMassBins(boost::dynamic_bitset<>& massBins, boost::dynamic_bitset<>& mzBins, Byte* massBinScores,
            double *filter, int filterSize, int binNumber, double tol, int minChargePeakCount){
        fill_n(massBinScores, binNumber,0);

        int binOffsets[filterSize];
        for(int i=0;i<filterSize;i++){
            binOffsets[i] = (filter[0] - filter[i]) / tol;
        }

        unsigned long setBinIndex = mzBins.find_first();

        while(setBinIndex != mzBins.npos){
            for(int j=0;j<filterSize;j++) {
                int bi = setBinIndex + binOffsets[j];
                if (bi >= binNumber) break;
                massBinScores[bi]++;

                if (massBinScores[bi] >= minChargePeakCount)
                    massBins[bi] = true;
            }
            setBinIndex = mzBins.find_next(setBinIndex);
        }

        fill_n(massBinScores, binNumber,0);
        setBinIndex = massBins.find_first();
        while(setBinIndex != massBins.npos){
            Byte continuousChargeCntr = 0;
            bool take = false;
            for(int j=0;j<filterSize;j++) {
                int bi = setBinIndex - binOffsets[j];
                if(bi<0) break;
                if(mzBins[bi]){
                    continuousChargeCntr++;
                    take = continuousChargeCntr >= minChargePeakCount;
                    massBinScores[setBinIndex] = massBinScores[setBinIndex] < continuousChargeCntr? continuousChargeCntr : massBinScores[setBinIndex];
                }else continuousChargeCntr = 0;
            }
            massBins[setBinIndex] = take;
            setBinIndex = massBins.find_next(setBinIndex);
        }
    }

    void filterMassBinsUsingIsotopeCriteria(boost::dynamic_bitset<>& massBins, Byte* massBinScores, int maxIsotopeIndex,
            int minIsotopeCount, int binNumber, double tol, double min){

        unsigned long setBinIndex = massBins.find_first();

        while(setBinIndex != massBins.npos){
            auto allIsoSet = true;
            auto m = exp(min + setBinIndex * tol);
            int* isoBins = new int[maxIsotopeIndex-1];
            fill_n(isoBins, maxIsotopeIndex-1,0);

            for(int i=1;i<maxIsotopeIndex;i++) {
                int bin = setBinIndex + Constants::C13C12_MASSDIFF_U*i/m/tol;
                if(bin>=binNumber - 1){
                    allIsoSet = false;
                    break;
                }
                bool isoSet = massBins[bin] | massBins[bin+1];
                if(i<minIsotopeCount) allIsoSet &= isoSet;//
                if(!isoSet) break;
                isoBins[i-1] = bin;
            }

            if(allIsoSet){
                for(int i=1;i<maxIsotopeIndex;i++) {
                    int bin = isoBins[i-1];
                    if(bin == 0) break;
                    auto score = massBinScores[setBinIndex];
                    score = score < massBinScores[bin]? massBinScores[bin] : score;
                    score = score < massBinScores[bin+1]? massBinScores[bin] : score;
                    massBinScores[setBinIndex] = score;
                    massBins[bin] = false;
                    if(bin < binNumber - 1) massBins[bin+1] = false;
                }
                if(setBinIndex > 0) massBins[setBinIndex-1] = false;
            }else{
                massBins[setBinIndex] = false;
            }
            setBinIndex = massBins.find_next(setBinIndex);
        }
    }

    double calculateBinScoreThreshold(boost::dynamic_bitset<>& massBins, Byte* massBinScores,
            int minChargePeakCount, int maxMassCountPerSpectrum, int filterSize){
        int* binScoreDist = new int[filterSize + 1];
        fill_n(binScoreDist,filterSize + 1,0);

        auto setBinIndex = massBins.find_first();
        while(setBinIndex != massBins.npos){
            binScoreDist[massBinScores[setBinIndex]]++;
            setBinIndex = massBins.find_next(setBinIndex);
        }

        int binScoreThreshold = minChargePeakCount, tsum=0;
        for(int i=filterSize;i>=0;i--){
            tsum+=binScoreDist[i];
            if(tsum >= maxMassCountPerSpectrum){
                binScoreThreshold = i > binScoreThreshold ? i : binScoreThreshold;
                break;
            }
        }
        return binScoreThreshold;
    }

    vector<vector<LogMzPeak>> findPeakGroups(vector<LogMzPeak> &logMzPeaks, double *filter, int filterSize, int minCharge,
                                            int maxIsotopeIndex, double tol, int minChargePeakCount, int minIsotopeCount,
                                            int maxMassCountPerSpectrum){
        vector<vector<LogMzPeak>> peakGroups;
        // minChargePeakCount should be < filterSize
        double min = logMzPeaks[0].logMz - filter[0]; // never fix it..
        double max = logMzPeaks[logMzPeaks.size()-1].logMz - filter[filterSize-minChargePeakCount];
        peakGroups.reserve(maxMassCountPerSpectrum*10);

        int binNumber = (max-min)/tol + 1;
        boost::dynamic_bitset<> mzBins = getMzBins(logMzPeaks, binNumber, tol);

        boost::dynamic_bitset<> massBins(binNumber);
        Byte* massBinScores = new Byte[binNumber];

        getMassBins(massBins, mzBins, massBinScores, filter, filterSize, binNumber, tol, minChargePeakCount);
        filterMassBinsUsingIsotopeCriteria(massBins, massBinScores, maxIsotopeIndex, minIsotopeCount, binNumber, tol, min);
        double binScoreThreshold = calculateBinScoreThreshold(massBins, massBinScores, minChargePeakCount, maxMassCountPerSpectrum, filterSize);

        int* currentPeakIndex = new int[filterSize];
        fill_n(currentPeakIndex,filterSize,0);

        auto setBinIndex = massBins.find_first();
        while(setBinIndex != massBins.npos){
            if(massBinScores[setBinIndex]<binScoreThreshold){
                setBinIndex = massBins.find_next(setBinIndex);
                continue;
            }
            int isoOff = 0;
            vector<LogMzPeak> peakGroup;
            peakGroup.reserve(logMzPeaks.size());

            for(int j=0;j<filterSize;j++) {
                int charge = j + minCharge;
                while(currentPeakIndex[j] < logMzPeaks.size()) {
                    double logMz = logMzPeaks[currentPeakIndex[j]].logMz;
                    int bi = (logMz - min - filter[j]) / tol;
                    if (bi>setBinIndex) break;
                    if (setBinIndex == bi) {
                        LogMzPeak p(*logMzPeaks[currentPeakIndex[j]].orgPeak);
                        p.charge = charge;
                        p.isotopeIndex = 0;
                        peakGroup.push_back(p);

                        if(currentPeakIndex[j]>0){
                            if(logMz - logMzPeaks[currentPeakIndex[j]-1].logMz < tol){
                                LogMzPeak p(*logMzPeaks[currentPeakIndex[j]-1].orgPeak);
                                p.charge = charge;
                                p.isotopeIndex = 0;
                                peakGroup.push_back(p);
                            }
                        }
                        int currentPeakIndexForIsotopes = currentPeakIndex[j] - 1;

                        for(int i=-1;i>-maxIsotopeIndex;i--) {
                            double logMzIsotope = log(exp(logMz) + Constants::C13C12_MASSDIFF_U*i/charge);

                            bool isotopePeakPresent = false;
                            while(currentPeakIndexForIsotopes >= 0){
                                double logMzForIsotope = logMzPeaks[currentPeakIndexForIsotopes].logMz;
                                currentPeakIndexForIsotopes--;
                                if(logMzForIsotope < logMzIsotope - tol) break;
                                if(logMzForIsotope > logMzIsotope + tol) continue;

                                isotopePeakPresent = true;

                                LogMzPeak p(*logMzPeaks[currentPeakIndexForIsotopes+1].orgPeak);
                                p.charge = charge;
                                p.isotopeIndex = i;
                                peakGroup.push_back(p);
                            }
                            if(!isotopePeakPresent) break;
                            isoOff = i > isoOff? isoOff : i;
                        }
                        currentPeakIndexForIsotopes = currentPeakIndex[j] + 1;

                        for(int i=1;i<maxIsotopeIndex + isoOff;i++) {
                            double logMzIsotope = log(exp(logMz) + Constants::C13C12_MASSDIFF_U*i/charge);

                            bool isotopePeakPresent = false;
                            while(currentPeakIndexForIsotopes < logMzPeaks.size()){
                                double logMzForIsotope = logMzPeaks[currentPeakIndexForIsotopes].logMz;
                                currentPeakIndexForIsotopes++;
                                if(logMzForIsotope < logMzIsotope - tol) continue;
                                if(logMzForIsotope > logMzIsotope + tol) break;

                                isotopePeakPresent = true;

                                LogMzPeak p(*logMzPeaks[currentPeakIndexForIsotopes-1].orgPeak);
                                p.charge = charge;
                                p.isotopeIndex = i;
                                peakGroup.push_back(p);
                            }
                            if(!isotopePeakPresent) break;
                        }
                    }
                    currentPeakIndex[j]++;
                }
            }

            vector<LogMzPeak> peakGroupIsotopeIndexAdjusted;
            peakGroupIsotopeIndexAdjusted.reserve(peakGroup.size());

            for(auto p : peakGroup){
                p.isotopeIndex = p.isotopeIndex - isoOff;
                if(p.isotopeIndex >= maxIsotopeIndex) continue;
                peakGroupIsotopeIndexAdjusted.push_back(p);
            }

            peakGroups.push_back(peakGroupIsotopeIndexAdjusted);
            setBinIndex = massBins.find_next(setBinIndex);
        }
        return peakGroups;
    }

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
    double getLogLikelihoodRatioScore(double int1, double int2){
        double sScore = getLogLikelihood(int1, int2, false);
        double nScore = getLogLikelihood(int1, int2, true);
       // cout<<int1 << " " << int2 << " " << pow(10, sScore) << " " << pow(10, nScore) << " " << sScore - nScore << endl;
        return sScore - nScore;
    }

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
};

int main(int argc, const char** argv)
{
    TDOnlineDeconv tool;
    return tool.main(argc, argv);
}
