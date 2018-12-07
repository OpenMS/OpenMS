//
// Created by JihyungKim on 07.12.18.
//

#include <OpenMS/config.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/MzMLFile.h>


using namespace OpenMS;
using namespace std;

class TDOnlineDeconv:
        public TOPPBase
{

public:
    TDOnlineDeconv():
            TOPPBase("TDOnlineDeconv", "Online Deconvolution for Smart MS2 acquisition with top down data", false)\
    {}

    struct {
        bool operator()(Peak1D a, Peak1D b) const
        {
            return a.getMZ() < b.getMZ();
        }
    } peakMzCompare;



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
        String infilePath = "/Users/kyowonjeong/Documents/A4B/mzml/05-26-17_B7A_yeast_td_fract12_rep1.mzML";

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
        onlineDeonvolution(map);
        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        std::cout << elapsed_secs << " seconds elapsed" << std::endl;
        std::cout << elapsed_secs/map.size()*1000 << " msec per spectrum"<< std::endl;

        return EXECUTION_OK;
    }


    void onlineDeonvolution(MSExperiment map){

        double threshold = 5e3;
        int filterSize = 25;
        int minCharge = 5;

        for (PeakMap::ConstIterator it = map.begin(); it != map.end(); ++it) {
            // check if this is a MS1 survey scan
            if (it->getMSLevel() != 1) continue;

            // RT value
            double rt = it->getRT();
            double filter[filterSize];

            for(int i=0;i<filterSize;i++){
                filter[i] = log(1.0/(i+minCharge));
            }

            Size peak_idx = 0;
            vector<Peak1D> peaks;
            peaks.reserve(it->size());

            for (; peak_idx < it->size(); ++peak_idx) {
                 Peak1D peak = (*it)[peak_idx];
                if(peak.getIntensity() < threshold) continue;

                Peak1D logmzpeak;
                 logmzpeak.setIntensity(peak.getIntensity());
                 logmzpeak.setMZ(log(peak.getMZ()));
                 peaks.push_back(logmzpeak);
//                file << (*it)[peak_idx].getMZ() << "," << rt << "," << (*it)[peak_idx].getIntensity() << ";";
            }

            vector<Peak1D> result;
            result.reserve(peaks.size() * filterSize);

            for(int i=0;i<filterSize;i++) {
                for(int j=0;j<peaks.size();j++){
                    Peak1D p;
                    p.setIntensity(peaks[j].getIntensity());
                    p.setMZ(peaks[j].getMZ() - filter[i]);
                    result.push_back(p);
                }
            }
           // cout<<result[100].getMZ()<<endl;
           // sort(result.begin(), result.end(), peakMzCompare);
           // cout<<result[100].getMZ()<<endl;
            //cout << result.size() << endl;
            double before;
            int cntr = 0;
            for(std::vector<Peak1D>::iterator iter= result.begin(); iter!=result.end(); ++iter ){
                Peak1D current = *iter;
                if(current.getMZ() - before < 1e-5){
                    cntr++;
                }
                before = current.getMZ();
            }

        }
    }
};

int main(int argc, const char** argv)
{
    TDOnlineDeconv tool;
    return tool.main(argc, argv);
}
