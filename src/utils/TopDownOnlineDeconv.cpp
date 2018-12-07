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
        String infilePath = "/Users/jeek/VirtualBox VMs/shared/05-26-17_B7A_yeast_td_fract12_rep1.mzML";

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

        onlineDeonvolution(map);

        return EXECUTION_OK;
    }

    void onlineDeonvolution(MSExperiment map){

        for (PeakMap::ConstIterator it = map.begin(); it != map.end(); ++it) {
            // check if this is a MS1 survey scan
            if (it->getMSLevel() != 1) continue;

            // RT value
            double rt = it->getRT();

            Size peak_idx = 0;
            for (; peak_idx < it->size(); ++peak_idx) {
//                file << (*it)[peak_idx].getMZ() << "," << rt << "," << (*it)[peak_idx].getIntensity() << ";";
            }

        }
    }
};

int main(int argc, const char** argv)
{
    TDOnlineDeconv tool;
    return tool.main(argc, argv);
}
