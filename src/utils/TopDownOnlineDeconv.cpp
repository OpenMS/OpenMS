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

class TopDownOnlineDeconvolution:
        public TOPPBase
{

public:
    TopDownOnlineDeconvolution():
            TOPPBase("TopDownOnlineDeconvolution", "Online Deconvolution for Smart MS2 acquisition with top down data", false)\
    {}

    typedef struct LogMzPeak{
        Peak1D *orgPeak;
        double logMz;
        int charge;
        int isotopeIndex;
        double score;

        LogMzPeak(): orgPeak(nullptr), logMz(-10000), charge(0), isotopeIndex(0), score(0) {}
        explicit LogMzPeak(Peak1D &peak): orgPeak(&peak), logMz(getLogMz(peak.getMZ())), charge(0), isotopeIndex(0), score(0) {}
        LogMzPeak(Peak1D &peak, int c, int i): orgPeak(&peak), logMz(getLogMz(peak.getMZ())), charge(c), isotopeIndex(i), score(0) {}

        double getMass(){
            return exp(logMz) * charge;
        }
        bool operator<(const LogMzPeak &a){
            return logMz < a.logMz;
        }
    } LogMzPeak;


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
        for(int r=1;r<=2;r++){
            ostringstream st;
            st <<  "/Users/kyowonjeong/Documents/A4B/matlab/yeast" << r << ".m";
            outfilePath = st.str();

            fstream fs;
            fs.open(outfilePath, fstream::out);
            fs << "m=[";
            for(int f=1;f<=12;f++){
                infilePath = infileDir + "f"+f+"r" + r+".mzML";
                cout << "%file name : " << infilePath << endl;
                MSExperiment map;
                mzml.load(infilePath, map);
                cout << "%Loaded consensus maps" << endl;
                clock_t begin = clock();
                onlineDeconvolution(map, fs, specCntr, qspecCntr, massCntr);
                clock_t end = clock();
                elapsed_secs += double(end - begin) / CLOCKS_PER_SEC;
                std::cout <<"%"<< qspecCntr << " MS1 spectra deconvoluted so far" <<  endl;
            }
            fs << "];";
            fs.close();
        }

        std::cout <<"%"<< elapsed_secs << " seconds elapsed for " <<specCntr << " MS1 spectra" <<  endl;
        std::cout <<"%"<< elapsed_secs/specCntr*1000 << " msec per spectrum"<< std::endl;

        return EXECUTION_OK;
    }

    void onlineDeconvolution(MSExperiment &map, fstream& fs, int& specCntr,int& qspecCntr,int& massCntr){
        double threshold = 1000;// the lower the more diverse molecules
        int filterSize = 25;
        int minCharge = 5;
        int minContinuousChargePeak = 3;//the lower the more spectrum deconved the running time sacrifice is huge..
        int minChargePeak = filterSize*.2;
        int minIsotopeCount = 3; //
        int maxMassPerSpectrum = 30;//  it should be checked at the end...
        int maxIsotopeIndex = 15; // make it variable
        double tolerance = 5e-6; // 5 ppm

        double filter[filterSize];
        double harmonicFilter[filterSize];

        auto generator= new CoarseIsotopePatternGenerator();
        generator->setMaxIsotope(maxIsotopeIndex);

        for(int i=0;i<filterSize;i++){
            filter[i] = log(1.0/(i + minCharge)); // should be descending, and negative!
            harmonicFilter[i] = log(1.0/(i+.5+minCharge));
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
                    findPeakGroups(logMzPeaks, filter, harmonicFilter, filterSize, minCharge, maxIsotopeIndex,
                            tolerance, minContinuousChargePeak, minChargePeak, minIsotopeCount, maxMassPerSpectrum);
            if (peakGroups.empty())
                continue;
            std::map<double, double> scoreMassMap;

            bool qualified = false;
            std::unordered_set<double> massSet;

            for(auto &pg : peakGroups){
                auto perChargeIntensities = new double[filterSize];
                auto perIsotopeIntensities = new double[maxIsotopeIndex + 1];

                fill_n(perChargeIntensities,filterSize,0);
                fill_n(perIsotopeIntensities,maxIsotopeIndex + 1,0); // change it to be flexible

                double monoIsotopeMass;
                double maxIntensityForMonoIsotopeMass = -1;
                double mostAbundantMass;
                double maxIntensityForMostAbundantMass = -1;

                for(auto &p : pg){
                    if(p.isotopeIndex>maxIsotopeIndex) continue;
                    int index = p.charge - minCharge;
                    double intensity = p.orgPeak->getIntensity();
                    perChargeIntensities[index] += intensity;
                    perIsotopeIntensities[p.isotopeIndex] += intensity;

                    if(maxIntensityForMostAbundantMass <= intensity) {
                        maxIntensityForMostAbundantMass = intensity;
                        mostAbundantMass = p.getMass();
                    }

                    if(p.isotopeIndex == 0) {
                        if(maxIntensityForMonoIsotopeMass > intensity) continue;
                        maxIntensityForMonoIsotopeMass = intensity;
                        monoIsotopeMass = p.getMass();
                    }
                }

                int setIntensityCounter = 0;
                int maxSetIntensityCounter = 0;
                for(int i=0;i<filterSize;i++){
                    if(perChargeIntensities[i]<=0){
                        setIntensityCounter = 0;
                        continue;
                    }
                    setIntensityCounter++;
                    maxSetIntensityCounter = maxSetIntensityCounter > setIntensityCounter? maxSetIntensityCounter : setIntensityCounter;
                }
                if(maxSetIntensityCounter < minContinuousChargePeak) continue;

                /*
                setIntensityCounter = 0;
                maxSetIntensityCounter = 0;
                for(int i=0;i<maxIsotopeIndex;i++){
                    //if(perIsotopeIntensities[i]<=0){
                        //setIntensityCounter = 0;
                        //continue;
                    //}
                    setIntensityCounter++;
                    maxSetIntensityCounter = maxSetIntensityCounter > setIntensityCounter? maxSetIntensityCounter : setIntensityCounter;

                }
                if(maxSetIntensityCounter < minIsotopeCount) continue;
*/
                auto avgIso = generator->estimateFromPeptideWeight(mostAbundantMass);

                int offset = 0;
                double maxCosine = -1;
                for(int f = -3;f<=3;f++){
                    auto cos = getCosine(perIsotopeIntensities, avgIso, .01, f);

                    if(maxCosine < cos){
                        maxCosine = cos;
                        offset = f;
                    }
                }

                if(maxCosine < .5){
                   continue;
                }

                monoIsotopeMass += offset * Constants::C13C12_MASSDIFF_U;

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
                if(score > 0 && !pg.empty()){
                    massSet.insert(monoIsotopeMass);
                    // cout<<massCntr++<<" "<<qspecCntr << " " << specCntr << " " <<score<< " " << monoIsotopeMass << endl;
                    qualified = true;
                }
            }
            if(!qualified) continue;

            qspecCntr++;
            /*std::unordered_set<int> nominalMasseSet;
            for(auto m :massSet){
                nominalMasseSet.insert(round(m*0.999497)); //float INTEGER_MASS_SCALER = 0.999497f;
            }
            for(auto m :massSet) {
                auto nm = round(m*2*0.999497);
                if(nominalMasseSet.count(nm)) massSet.erase(m);
            }*/
            for(auto m :massSet) {
                fs<<massCntr++<<" "<<qspecCntr << " " << specCntr << " " << it->getRT() <<" " << setprecision(15)<< m << " " << round(m*0.999497) << endl;
                //cout<<massCntr++<<" "<<qspecCntr << " " << specCntr << " " << it->getRT() <<" " << setprecision(15)<< m << " " << round(m*0.999497) << endl;
            }

           // for (int i=0; i<20 && i<scoreMassMap.size() ; ++iter, i++){
           //     cout << qspecCntr<<" "<<specCntr<<" "<<iter->first << " " << iter->second << " " << peakGroups.size() << " " << scoreMassMap.size()<< endl;
           // }
           // cout << endl;
        }
        return;
    }


    double getCosine(double* a, IsotopeDistribution b, double rightThreshold, int offset = 0){
        double n=0, d1=0, d2=0;
        auto bMax = b.getMostAbundant();
        for(int i=0;i<b.size();i++){
            int j = i + offset;
            if(j<0 || j>=b.size()) continue;
            double bInt = b[i].getIntensity();
            if(b[i].getMZ()>bMax.getMZ() && bInt/bMax.getIntensity() < rightThreshold) break;
            n += a[j] * bInt;
            d1 += a[j] * a[j];
            d2 += bInt * bInt;
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
            double *filter, double*harmonicFilter, int filterSize, int binNumber, double tol, int minContinuousChargePeak, int minChargePeak){
        fill_n(massBinScores, binNumber,0);

        int binOffsets[filterSize];
        int harmonicBinOffsets[filterSize];

        for(int i=0;i<filterSize;i++){
            binOffsets[i] = (filter[0] - filter[i]) / tol;
            harmonicBinOffsets[i] = (filter[0] - harmonicFilter[i]) / tol;
        }

        unsigned long setBinIndex = mzBins.find_first();

        while(setBinIndex != mzBins.npos){
            for(int j=0;j<filterSize;j++) {
                int bi = setBinIndex + binOffsets[j];
                if (bi >= binNumber) break;
               // int hbi = setBinIndex + harmonicBinOffsets[j];
               // if(hbi>=binNumber - 1 || !mzBins[hbi])
               massBinScores[bi]++;
                if (massBinScores[bi] >= minChargePeak) {
                    massBins[bi] = true;
                }
            }
            setBinIndex = mzBins.find_next(setBinIndex);
        }

        fill_n(massBinScores, binNumber,0);
        setBinIndex = massBins.find_first();
        while(setBinIndex != massBins.npos){
        //for(int setBinIndex : setBins){
            Byte continuousChargeCntr = 0;
            bool take = false;
            for(int j=0;j<filterSize;j++) {
                int bi = setBinIndex - binOffsets[j];
                if(bi<0) break;
                int hbi = setBinIndex - harmonicBinOffsets[j];
                if(mzBins[bi] && (hbi<0 || !mzBins[hbi])){
                    continuousChargeCntr++;
                    take = continuousChargeCntr >= minContinuousChargePeak;
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
            int* isoBins = new int[maxIsotopeIndex];
            fill_n(isoBins, maxIsotopeIndex,0);

            for(int i=1;i<=maxIsotopeIndex;i++) {
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
                for(int i=1;i<=maxIsotopeIndex;i++) {
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


    void updatePeakGroups(vector<vector<LogMzPeak>>& peakGroups, boost::dynamic_bitset<>& massBins, Byte* massBinScores,
                                    int binScoreThreshold, vector<LogMzPeak> &logMzPeaks, int filterSize, int minCharge,
                                    int maxIsotopeIndex, double min, double* filter, double tol){

        int* currentPeakIndex = new int[filterSize];
        fill_n(currentPeakIndex,filterSize,0);

        auto setBinIndex = massBins.find_first();
        while(setBinIndex != massBins.npos){
            if(massBinScores[setBinIndex]>=binScoreThreshold) {

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
                            LogMzPeak p(*logMzPeaks[currentPeakIndex[j]].orgPeak, charge, 0);
                            peakGroup.push_back(p);

                            if(currentPeakIndex[j]>0){
                                if(logMz - logMzPeaks[currentPeakIndex[j]-1].logMz < tol){
                                    LogMzPeak p(*logMzPeaks[currentPeakIndex[j]-1].orgPeak, charge, 0);
                                    massBins[(p.logMz - min - filter[j]) / tol] = false;
                                    peakGroup.push_back(p);
                                }
                            }

                            for(int d=-1;d<=1;d+=2){ // negative then positive direction.
                                int currentPeakIndexForIsotopes = currentPeakIndex[j] + d;
                                for(int i=1;i<=maxIsotopeIndex + (d<0? 0 : isoOff);i++) {
                                    double logMzIsotope = log(exp(logMz) + Constants::C13C12_MASSDIFF_U*i*d/charge);

                                    bool isotopePeakPresent = false;
                                    while(currentPeakIndexForIsotopes >= 0 && currentPeakIndexForIsotopes < logMzPeaks.size()){
                                        double logMzForIsotope = logMzPeaks[currentPeakIndexForIsotopes].logMz;
                                        currentPeakIndexForIsotopes += d;
                                        if(logMzForIsotope < logMzIsotope - tol) if(d<0) break; else continue;
                                        if(logMzForIsotope > logMzIsotope + tol) if(d<0) continue; else break;

                                        isotopePeakPresent = true;

                                        LogMzPeak p(*logMzPeaks[currentPeakIndexForIsotopes-d].orgPeak, charge, i*d);
                                        if(d > 0) massBins[(p.logMz - min - filter[j]) / tol] = false;
                                        peakGroup.push_back(p);
                                    }
                                    if(!isotopePeakPresent) break;
                                    if(d < 0) isoOff = -i > isoOff? isoOff : -i;
                                }
                            }
                        }
                        currentPeakIndex[j]++;
                    }
                }

                for(LogMzPeak &p : peakGroup){
                    p.isotopeIndex = p.isotopeIndex - isoOff;
                }

                peakGroups.push_back(peakGroup);
            }
            setBinIndex = massBins.find_next(setBinIndex);
        }
    }

    vector<vector<LogMzPeak>> findPeakGroups(vector<LogMzPeak> &logMzPeaks, double *filter, double* harmonicFilter, int filterSize, int minCharge,
                                            int maxIsotopeIndex, double tol, int minContinuousChargePeak, int minChargePeak, int minIsotopeCount,
                                            int maxMassCountPerSpectrum){
        vector<vector<LogMzPeak>> peakGroups;
        // minContinuousChargePeak should be < filterSize
        double min = logMzPeaks[0].logMz - filter[0]; // never fix it..
        double max = logMzPeaks[logMzPeaks.size()-1].logMz - filter[filterSize-minChargePeak];
        peakGroups.reserve(maxMassCountPerSpectrum*10);

        int binNumber = (max-min)/tol + 1;
        boost::dynamic_bitset<> mzBins = getMzBins(logMzPeaks, binNumber, tol);
        boost::dynamic_bitset<> massBins(binNumber);
        Byte* massBinScores = new Byte[binNumber];

        getMassBins(massBins, mzBins, massBinScores, filter, harmonicFilter, filterSize, binNumber, tol, minContinuousChargePeak, minChargePeak);

        filterMassBinsUsingIsotopeCriteria(massBins, massBinScores, maxIsotopeIndex, minIsotopeCount, binNumber, tol, min);
        double binScoreThreshold = calculateBinScoreThreshold(massBins, massBinScores, minContinuousChargePeak, maxMassCountPerSpectrum, filterSize);

        updatePeakGroups(peakGroups, massBins, massBinScores, binScoreThreshold, logMzPeaks, filterSize, minCharge,
                maxIsotopeIndex, min, filter, tol);

        return peakGroups;
    }

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
    TopDownOnlineDeconvolution tool;
    return tool.main(argc, argv);
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