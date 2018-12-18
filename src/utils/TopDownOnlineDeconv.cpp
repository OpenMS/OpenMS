//
// Created by JihyungKim on 07.12.18.
//

#include <OpenMS/config.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <queue>

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
        const double protonMass = 1.0072765;
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
        String infilePath = "/Users/kyowonjeong/Documents/A4B/mzml/05-26-17_B7A_yeast_td_fract12_rep1.mzML";
      //  String infilePath = "/Users/kyowonjeong/Documents/A4B/mzml/180523_Cytocrome_C_MS2_HCD.mzML";
      //  String infilePath = "/Users/kyowonjeong/Documents/A4B/mzml/180523_Myoglobin_MS2_HCD.mzML";
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
        int filterSize = 25;
        int minCharge = 5;
        int minChargePeakCount = filterSize*.5; // inclusive
        int massCount = 20;
        double tolerance = 5e-6;

        int specCntr = 0;
        double filter[filterSize];

        for(int i=0;i<filterSize;i++){
            filter[i] = log(1.0/(i + minCharge)); // should be descending, and negative!
        }

        for (auto it = map.begin(); it != map.end(); ++it) {
            if (it->getMSLevel() != 1) continue;
            specCntr++;

            double rt = it->getRT();
            //auto filteredSpec = *it;//filterSpectrum(*it, 1, 2, threshold);

            vector<LogMzPeak> logMzPeaks;
            logMzPeaks.reserve(it->size());
            for (auto &peak : *it){
                if(peak.getIntensity() <= threshold) continue;
                LogMzPeak logmzpeak(peak);
                logMzPeaks.push_back(logmzpeak);
            }

            //continue;
            auto peakGroups =
                    findPeakGroups(logMzPeaks, filter, filterSize, minCharge, tolerance, minChargePeakCount, massCount);
            if (peakGroups.empty()) continue;

            std::map<double, double> scoreMassMap;

            for(auto &pg : peakGroups){
                double intensities[filterSize];
                int maxIntensityIndex;
                double maxIntensity = -1;

                for(int k=0;k<filterSize;k++){
                    intensities[k] = 0;
                }
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
                                             double tol, int minChargePeakCount, int massCount){
        vector<vector<LogMzPeak>> peakGroups;
        // minChargePeakCount should be < filterSize
        double min = spectrum[0].logmz - filter[minChargePeakCount-1];
        double max = spectrum[spectrum.size()-1].logmz - filter[filterSize-minChargePeakCount];
        peakGroups.reserve(massCount*100);

        int binNumber = (max-min)/tol;
        if(binNumber <=0) return peakGroups;

        Byte bins[binNumber + 1];

        for(int i=0;i<binNumber+1;i++){
            bins[i] = 0;
        }

        double binf[filterSize];
        for(int i=0;i<filterSize;i++){
            binf[i] = filter[i]/tol;
        }
        for(auto &p : spectrum) {
            double logMz = p.logmz;
            double fbi = (logMz - min)/tol;
            for(int j=0;j<filterSize;j++){
                int bi = fbi - binf[j];
                if(bi<1) continue;
                if(bi+1>binNumber) continue;
                bins[bi-1]++;
                bins[bi]++;
            }
        }

        int currentPeakIndex[filterSize];
        int binScoreDist[filterSize];
        for(int i=0;i<filterSize;i++){
            currentPeakIndex[i] = binScoreDist[i] = 0;
        }

        for(int i=0;i<binNumber+1;i++) {
            if(bins[i] < minChargePeakCount) continue;
            int j = bins[i] >= filterSize? filterSize - 1 : bins[i];
            binScoreDist[j]++;
        }

        int binScoreThreshold = minChargePeakCount, tsum=0;
        for(int i=filterSize-1;i>=0;i--){
            tsum+=binScoreDist[i];
            if(tsum >= massCount){
                binScoreThreshold = i > binScoreThreshold ? i  : binScoreThreshold;
                break;
            }
        }

        for(int i=0;i<binNumber+1;i++){
            if(bins[i]<binScoreThreshold) continue;
            map<int, int> toselect; // index - > charge

            for(int j=0;j<filterSize;j++) {
                while(currentPeakIndex[j] < spectrum.size()) {
                    double logMz = spectrum[currentPeakIndex[j]].logmz;
                    int bi = (logMz - min) / tol - binf[j];
                    if (bi>i) break;
                    if (i == bi) {
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
        }
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
