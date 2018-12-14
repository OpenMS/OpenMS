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
        double protonMass = 1.0072765;

        LogMzPeak(): orgPeak(nullptr), logmz(-10000), charge(0){}
        LogMzPeak(Peak1D &peak): orgPeak(&peak), logmz(log(peak.getMZ() - protonMass)), charge(0){}
        LogMzPeak(const LogMzPeak &other, double shift, int charge): orgPeak(other.orgPeak), logmz(shift + other.logmz), charge(charge){}

        double getMass(){
            return exp(logmz);
        }
        bool operator<(const LogMzPeak &a){
            return logmz < a.logmz;
        }
    };


protected:
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
       // String infilePath = "/Users/kyowonjeong/Documents/A4B/mzml/05-26-17_B7A_yeast_td_fract12_rep1.mzML";
        String infilePath = "/Users/kyowonjeong/Documents/A4B/mzml/180523_Cytocrome_C_MS2_HCD.mzML";
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
        double threshold = 5e3;
        int filterSize = 25;
        int minCharge = 5;
        int minPeakNumForScoring = 3; // inclusive

        int specCntr = 0;
        double filter[filterSize];
        double tolerance = 1e-5;
        int maxPeakCntr = 2000;

        for(int i=0;i<filterSize;i++){
            filter[i] = log(1.0/(i+minCharge));
        }

        for (auto it = map.begin(); it != map.end(); ++it) {
            if (it->getMSLevel() != 1) continue;
            specCntr++;


            double rt = it->getRT();
            double currentThreshold = threshold;

            if(it->size() > maxPeakCntr){
                vector<Peak1D> allPeaks;
                allPeaks.reserve(it->size());
                for (auto &peak : (*it) ) {
                    allPeaks.push_back(peak);
                }
                while(allPeaks.size() > maxPeakCntr) {
                    vector<Peak1D> tmpPeaks;
                    tmpPeaks.reserve(allPeaks.size());
                    for (auto &peak : allPeaks){
                        if (peak.getIntensity() < currentThreshold) continue;
                        tmpPeaks.push_back(peak);
                    }
                    if (tmpPeaks.size() <= maxPeakCntr) {
                        break;
                    }
                    allPeaks = tmpPeaks;
                    currentThreshold = currentThreshold * 1.2;
                }
            }

            auto filteredSpec = filterSpectrum(*it, 2, 1, threshold);
           // auto filteredSpec = *it;

            //cout<<filteredSpec.size()<<endl;
            vector<LogMzPeak> logMzPeaks;
            logMzPeaks.reserve(filteredSpec.size());

            for (auto &peak : filteredSpec ){
               // if (peak.getIntensity() < currentThreshold) continue;

                LogMzPeak logmzpeak(peak);
                logMzPeaks.push_back(logmzpeak);
            }

            vector<vector<LogMzPeak>> matrix(filterSize); // initialization

            for (int i = 0; i < filterSize; i++) {
                matrix[i].reserve(logMzPeaks.size());
                for (auto &logMzPeak : logMzPeaks) {
                    LogMzPeak p(logMzPeak, -filter[i], i+minCharge);
                    matrix[i].push_back(p);
                }
            }

            vector<LogMzPeak> result;
            result.reserve(filterSize * logMzPeaks.size());
            sortMatrix(matrix, result);


            vector<LogMzPeak> groupedPeaks;
            groupedPeaks.reserve(100);
            groupedPeaks.push_back(result[0]);

            std::map<double, double> scoreMassMap;

            for (int i=1; i<result.size();i++){
                auto peak = result[i];
                if (peak.logmz - groupedPeaks[0].logmz < tolerance * 2) {
                    groupedPeaks.push_back(peak);
                }else{
                    if(groupedPeaks.size()>=minPeakNumForScoring){
                        LogMzPeak peaksToScore[filterSize];
                        int maxIntensityIndex;
                        double maxIntensity = -1;

                        for(int k=0;k<filterSize;k++){
                            peaksToScore[k] = LogMzPeak();
                        }
                        for(auto &gp : groupedPeaks){
                            int index = gp.charge - minCharge;
                            double intensity = gp.orgPeak->getIntensity();
                            peaksToScore[index] = gp;
                            if(maxIntensity < intensity){
                                maxIntensity = intensity;
                                maxIntensityIndex = index;
                            }
                        }
                        double score = 0;
                        for(int k=maxIntensityIndex;k<filterSize-1;k++){
                            double int1 = peaksToScore[k].orgPeak == nullptr? 0 : peaksToScore[k].orgPeak->getIntensity();
                            double int2 = peaksToScore[k+1].orgPeak == nullptr? 0 : peaksToScore[k+1].orgPeak->getIntensity();
                            score += getLogLikelihoodRatioScore(int1, int2);
                        }
                        for(int k=maxIntensityIndex;k>0;k--){
                            double int1 = peaksToScore[k].orgPeak == nullptr? 0 : peaksToScore[k].orgPeak->getIntensity();
                            double int2 = peaksToScore[k-1].orgPeak == nullptr? 0 : peaksToScore[k-1].orgPeak->getIntensity();
                            score += getLogLikelihoodRatioScore(int1, int2);
                        }
                        if(score > 10 && !groupedPeaks.empty()) scoreMassMap[score] = groupedPeaks[groupedPeaks.size()/2].getMass();
                        //if(score > log10(2)) cout<<groupedPeaks[groupedPeaks.size()/2].getMass() << "  " << score<<endl;
                    }
                    groupedPeaks.clear();
                    groupedPeaks.push_back(peak);
                }
            }

            if(scoreMassMap.empty()) continue;
          /*  auto iter = scoreMassMap.rbegin();
            for (int i=0; i<5 && i<scoreMassMap.size() ; ++iter, i++){
                cout << iter->first << " " << iter->second << endl;
            }
            cout << endl;*/
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
        //make_heap(begin(intensityHeap), end(intensityHeap));
        vector<Peak1D> initFiltered;
        initFiltered.reserve(spectrum.size());
        for(auto p : spectrum) if(p.getIntensity() > th) initFiltered.push_back(p);

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
           // cout<< spectrum[wsIndex].getMZ() << " " << mz << " " << spectrum[weIndex].getMZ()<< " " << changed << endl;
            if(p.getIntensity() > median * factor)
                filtered.push_back(p);

            prevMedian = median;
        }
        return filtered;
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
