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

using namespace OpenMS;
using namespace std;

class TDOnlineDeconv:
        public TOPPBase
{

public:
    TDOnlineDeconv():
            TOPPBase("TDOnlineDeconv", "Online Deconvolution for Smart MS2 acquisition with top down data", false)\
    {}

    typedef pair<double, pair<Size, Size>> ppi;

    typedef struct LogMzPeak{
        Peak1D *orgPeak;
        double logmz;
        int charge;
        double score;

        LogMzPeak(): orgPeak(nullptr), logmz(-10000), charge(0){}
        LogMzPeak(Peak1D &peak): orgPeak(&peak), logmz(getLogMz(peak.getMZ())), charge(0){}

        double getMass(){
            return exp(logmz) * charge;
        }
        bool operator<(const LogMzPeak &a){
            return logmz < a.logmz;
        }
    };


protected:

    static double getLogMz(double mz){
        const double protonMass = 1.0072764668;// 1.0072765;
        return log(mz - protonMass);
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
      //  String infilePath = "/Users/kyowonjeong/Documents/A4B/mzml/05-26-17_B7A_yeast_td_fract12_rep1.mzML";
        String infilePath = "/Users/kyowonjeong/Documents/A4B/mzml/180523_Cytocrome_C_MS2_HCD.mzML";
     //   String infilePath = "/Users/kyowonjeong/Documents/A4B/mzml/180523_Myoglobin_MS2_HCD.mzML";
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
        int cntr = onlineDeonvolution(map);
        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        std::cout << elapsed_secs << " seconds elapsed" << std::endl;
        std::cout << elapsed_secs/cntr*1000 << " msec per spectrum"<< std::endl;

        return EXECUTION_OK;
    }

    int onlineDeonvolution(MSExperiment &map){
        double threshold = 5000;
        int filterSize = 15;
        int minCharge = 5;
        int minChargePeakCount = filterSize*0.3; // inclusive
        int massCount = 20; // TODO - count isotope cluster a single mass! and reduce this to like 5..
        double tolerance = 5e-6; // 5 ppm

        int specCntr = 0;
        double filter[filterSize];

        for(int i=0;i<filterSize;i++){
            filter[i] = log(1.0/(i + minCharge)); // should be descending, and negative!
        }

        ofstream mfile;
        mfile.open ("/Users/kyowonjeong/Documents/A4B/matlab/bins.m");


        for (auto it = map.begin(); it != map.end(); ++it) {
            if (it->getMSLevel() != 1) continue;
           // if(it->size()<1000) continue;
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
                    findPeakGroups(logMzPeaks, filter, filterSize, minCharge, tolerance, minChargePeakCount, massCount, mfile, specCntr);
            if (peakGroups.empty()) continue;
            std::map<double, double> scoreMassMap;

            for(auto &pg : peakGroups){
                double* intensities = new double[filterSize];
                int maxIntensityIndex;
                double maxIntensity = -1;

                fill_n(intensities,filterSize,0);

                for(auto &p : pg){
                    int index = p.charge - minCharge;
                    double intensity = p.orgPeak->getIntensity();
                    intensities[index] += intensity;
                    if(maxIntensity < intensities[index]){
                        maxIntensity = intensities[index];
                        maxIntensityIndex = index;
                    }
                }
                double score = 0;
                for(int k=maxIntensityIndex;k<filterSize-1;k++){
                    score += getLogLikelihoodRatioScore(intensities[k], intensities[k+1]);
                }
                for(int k=maxIntensityIndex;k>0;k--){
                    score += getLogLikelihoodRatioScore(intensities[k], intensities[k-1]);
                }
                if(score > 10 && !pg.empty()) scoreMassMap[score] = pg[pg.size()/2].getMass();
            }

            if(scoreMassMap.empty()) continue;
            auto iter = scoreMassMap.rbegin();

            for (int i=0; i<20 && i<scoreMassMap.size() ; ++iter, i++){
                cout << iter->first << " " << iter->second << " " << peakGroups.size() << " " << scoreMassMap.size()<< endl;
            }
            cout << endl;
        }
        mfile.close();
        return specCntr;
    }

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
    }

    vector<vector<LogMzPeak>> findPeakGroups(vector<LogMzPeak> &spectrum, double *filter, int filterSize, int minCharge,
                                             double tol, int minChargePeakCount, int massCount, ofstream& mfile, int specCntr){
        vector<vector<LogMzPeak>> peakGroups;
        // minChargePeakCount should be < filterSize
        double min = spectrum[0].logmz - filter[0]; // never fix it..
        double max = spectrum[spectrum.size()-1].logmz - filter[filterSize-minChargePeakCount];
        peakGroups.reserve(massCount*100);

        int binNumber = (max-min)/tol;
        if(binNumber <=0) return peakGroups;
        boost::dynamic_bitset<> mzBins(binNumber + 1);

        for (auto &p : spectrum) {
            int bi = (p.logmz - spectrum[0].logmz) / tol;
            if (bi >= binNumber) break;
            mzBins[bi] = true;
        }
        mzBins = mzBins>>1 | mzBins;

        if(mzBins.count() == 0) return peakGroups;

        boost::dynamic_bitset<> massBins(binNumber + 1);

        Byte* massBinScores = new Byte[binNumber + 1];
        fill_n(massBinScores, binNumber + 1,0);

        int binf[filterSize];
        for(int i=0;i<filterSize;i++){
            binf[i] = (filter[0] - filter[i]) / tol;
        }

        int index = mzBins.find_first();

        while(index != mzBins.npos){
            for(int j=0;j<filterSize;j++) {
                int bi = (index + binf[j]);
                if (bi >= binNumber) break;
                massBinScores[bi]++;

                if (massBinScores[bi] >= minChargePeakCount)
                    massBins[bi] = true;
            }
            index = mzBins.find_next(index);
        }

        if(massBins.count() == 0) return peakGroups;

        double isom = 1.003355;
        int isoMinCntr = 4;

        index = massBins.find_first();
        while(index != massBins.npos){
            auto isopresent = true;
            auto m = exp(min + index * tol);
            for(int i=1;i<isoMinCntr;i++) {
                int bin = index + isom*i/m/tol;
                isopresent &= massBins[bin];
            }

            if(isopresent){
// meaning of bin??    remove massBins here?
            }else{
                massBins[index] = false;
            }
            index = massBins.find_next(index);
        }

        if(massBins.count() == 0) return peakGroups;

        int* currentPeakIndex = new int[filterSize];
        fill_n(currentPeakIndex,filterSize,0);

        int* binScoreDist = new int[filterSize];
        fill_n(binScoreDist,filterSize,0);

        index = massBins.find_first();
        while(index != massBins.npos){
            int j = massBinScores[index] >= filterSize? filterSize - 1 : massBinScores[index];
            binScoreDist[j]++;
            index = massBins.find_next(index);
        }

        if(massBins.count() == 0) return peakGroups;

        int binScoreThreshold = minChargePeakCount, tsum=0;
        for(int i=filterSize;i>=1;i--){
            tsum+=binScoreDist[i-1];
            if(tsum >= massCount){
                binScoreThreshold = i + 1 > binScoreThreshold ? i + 1 : binScoreThreshold;
                binScoreThreshold = binScoreThreshold > filterSize ? filterSize : binScoreThreshold;
                break;
            }
        }
        /*if(specCntr==602) {
            mfile << "mass" << specCntr << "=[";
            index = mzBins.find_first();

            while(index != mzBins.npos){
                for(int j=0;j<filterSize;j++) {
                    int bi = (index + binf[j]);
                    if (bi >= binNumber) break;

                    mfile << setprecision(15) << exp(spectrum[0].logmz - filter[0] + bi * tol) << "," << j << ";";

                }
                index = mzBins.find_next(index);
            }


            mfile << "];" << endl;

            mfile << "spec" << specCntr<<"=[";
            for (auto &p : spectrum) {
                mfile<<p.orgPeak->getMZ()<<","<<p.orgPeak->getIntensity()<<","<<setprecision(15)<<p.logmz<<";";
            }
            mfile << "];" << endl;
        }*/
        index = massBins.find_first();
        while(index != massBins.npos){
            if(massBinScores[index]<binScoreThreshold){
                index = massBins.find_next(index);
                continue;
            }
            map<int, int> toselect; // index - > charge

            for(int j=0;j<filterSize;j++) {
                while(currentPeakIndex[j] < spectrum.size()) {
                    double logMz = spectrum[currentPeakIndex[j]].logmz;
                    int bi = (logMz - min - filter[j]) / tol;
                    if (bi>index) break;
                    if (index == bi) {
                        toselect[currentPeakIndex[j]] = j + minCharge;
                        if(currentPeakIndex[j]>0){
                            if(logMz - spectrum[currentPeakIndex[j]-1].logmz < tol){
                                toselect[currentPeakIndex[j]-1] = j + minCharge;
                            }
                        }
                    }
                    currentPeakIndex[j]++;
                }
            }
            vector<LogMzPeak> logpeaks;
            logpeaks.reserve(toselect.size());

            for(auto iter = toselect.begin(); iter != toselect.end(); ++ iter){
                auto t = iter->first;
                spectrum[t].charge = iter->second; // charge should be multiple...
                logpeaks.push_back(spectrum[t]);
            }
            peakGroups.push_back(logpeaks);
            index = massBins.find_next(index);
        }
       // cout<<binScoreThreshold<<endl;
        //cout<<peakGroups.size()<<endl;
        return peakGroups;
    }


    void sortMatrix(vector<vector<LogMzPeak>> &matrix, vector<LogMzPeak> &result){
        priority_queue< ppi, vector<ppi>, greater<ppi> > pq;

        for (Size i=0; i< matrix.size(); i++){
            pq.push({matrix[i][0].logmz, {i, 0}});
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
                pq.push({ matrix[i][j + 1].logmz, { i, j + 1 } });
        }
    }

    double getLogLikelihoodRatioScore(double int1, double int2){
        return getLogLikelihood(int1, int2, false) - getLogLikelihood(int1, int2, true);
    }

    double getLogLikelihood(double int1, double int2, bool isH0){
        double tmp = 1e-4;
        double ret;
        if(int1<=0){
            if(int2<=0){
                if(isH0) ret = 1-tmp;
                else ret = 1 - 2*tmp;
            }else{
                if(isH0) ret = tmp;
                else ret = 2*tmp;
            }
        }else{
            if(int2<=0){
                if(isH0) ret = 1-tmp;
                else ret = .1;
            }else if(int1<int2){
                if(isH0) ret = tmp/2;
                else ret = .05;
            }else {
                if (isH0) ret = tmp / 2;
                else ret = .85;
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
