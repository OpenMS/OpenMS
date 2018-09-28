//
// Created by JihyungKim on 25.09.18.
//

#include <OpenMS/config.h>

#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/CachedMzML.h>
#include <OpenMS/FORMAT/HANDLERS/CachedMzMLHandler.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>
#include <OpenMS/FORMAT/IndexedMzMLFileLoader.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>


using namespace OpenMS;
using namespace std;

class TDfeatureDeconvolution :
    public TOPPBase
{
public:
    TDfeatureDeconvolution():
        TOPPBase("TDfeatureDeconvolution", "Feature deconvolution with top down data", false)
    {
    }

protected:
    // this function will be used to register the tool parameters
    // it gets automatically called on tool execution
    void registerOptionsAndFlags_() override
    {
        registerInputFile_("in", "<file>", "", "Input file for deconvolution.");
        setValidFormats_("in", ListUtils::create<String>("mzML"));
    }

    // the main_ function is called after all parameters are read
    ExitCodes main_(int, const char **) override
    {
        //-------------------------------------------------------------
        // parsing parameters
        //-------------------------------------------------------------
        String infilePath = getStringOption_("in");

        cout << "file name : " << infilePath << endl;

        //-------------------------------------------------------------
        // reading input
        //-------------------------------------------------------------
        MzMLFile mzml;
        mzml.setLogType(log_type_);

        // load input
        MSExperiment map;
        mzml.load(infilePath, map);
        cout << "Loaded consensus maps" << endl;

        //-------------------------------------------------------------
        // calculations
        //-------------------------------------------------------------
        // mass trace part
        cout << map.size() << endl;




        //-------------------------------------------------------------
        // writing output
        //-------------------------------------------------------------
    }
};

int main(int argc, const char** argv)
{
    TDfeatureDeconvolution tool;
    return tool.main(argc, argv);
}
