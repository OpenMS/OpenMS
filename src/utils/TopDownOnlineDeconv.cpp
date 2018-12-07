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
        String infilePath = "/Users/jeek/Documents/Thermo/05-26-17_B7A_yeast_td_fract12_rep1.mzML";
//        String infilePath = "/Users/jeek/Documents/A4B/1st Set/TD/180523_Cytocrome_C_MS2_HCD.mzML";

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
        onlineDeconvolution(map);
        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        std::cout << elapsed_secs << " seconds elapsed" << std::endl;
        std::cout << elapsed_secs/map.size()*1000 << " msec per spectrum"<< std::endl;

        return EXECUTION_OK;
    }

    void onlineDeconvolution(MSExperiment &map) {

        double threshold = 5e3;
        int filterSize = 25;
        int minCharge = 5;

        double filter[filterSize];
        for (int i = 0; i < filterSize; i++) {
            filter[i] = log(1.0 / (i + minCharge));
        }

        for (PeakMap::ConstIterator it = map.begin(); it != map.end(); ++it) {
            // check if this is a MS1 survey scan
            if (it->getMSLevel() != 1) continue;

            // RT value
//            double rt = it->getRT();

            vector<Peak1D> peaks;
            peaks.reserve(it->size());

            // create vector with log(MZ) & filtered intensity
            for (auto &peak : (*it) ) {
                if (peak.getIntensity() < threshold) continue;

                Peak1D logmzpeak;
                logmzpeak.setIntensity(peak.getIntensity());
                logmzpeak.setMZ(log(peak.getMZ()));
                peaks.push_back(logmzpeak);
            }

            vector<vector<Peak1D>> matrix(filterSize); // initialization

            for (int i = 0; i < filterSize; i++) {
                matrix[i].reserve(peaks.size());
                for (auto &peak : peaks) {
                    Peak1D p;
                    p.setIntensity(peak.getIntensity());
                    p.setMZ(peak.getMZ() - filter[i]);
                    matrix[i].push_back(p);
                }
            }

            vector<Peak1D> result;
            result.reserve(filterSize * peaks.size());
            sortMatrix(matrix, result);

            double before;
            int cntr = 0;
            for (std::vector<Peak1D>::iterator iter = result.begin(); iter != result.end(); ++iter) {
                Peak1D current = *iter;
                if (current.getMZ() - before < 1e-5) {
                    cntr++;
                }
                before = current.getMZ();
            }
        }
    }
    void sortMatrix(vector<vector<Peak1D>> &matrix, vector<Peak1D> &result){
        priority_queue< ppi, vector<ppi>, greater<ppi> > pq;

        for (Size i=0; i< matrix.size(); i++){
            pq.push({matrix[i][0].getMZ(), {i, 0}});
        }

        while (pq.empty() == false) {
            ppi curr = pq.top();
            pq.pop();

            // i ==> Array Number (filter index)
            // j ==> Index in the array number (peak index)
            Size i = curr.second.first;
            Size j = curr.second.second;

            result.push_back(matrix[i][j]);

            // The next element belongs to same array as current.
            if (j + 1 < matrix[i].size())
                pq.push({ matrix[i][j + 1].getMZ(), { i, j + 1 } });
        }
    }
};

int main(int argc, const char** argv)
{
    TDOnlineDeconv tool;
    return tool.main(argc, argv);
}
