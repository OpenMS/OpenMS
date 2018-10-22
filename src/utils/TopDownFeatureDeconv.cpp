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
#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/KERNEL/Peak2D.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>

#include <OpenMS/CONCEPT/Types.h>

#include <OpenMS/FILTERING/DATAREDUCTION/MassTraceDetection.h>



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

        vector<MassTrace> m_traces;

        //-------------------------------------------------------------
        // set parameters
        //-------------------------------------------------------------

        Param common_param = getParam_().copy("algorithm:common:", true);
        writeDebug_("Common parameters passed to sub-algorithms (mtd and ffm)", common_param, 3);

        Param mtd_param = getParam_().copy("algorithm:mtd:", true);
        writeDebug_("Parameters passed to MassTraceDetection", mtd_param, 3);


        //-------------------------------------------------------------
        // configure and run mass trace detection
        //-------------------------------------------------------------

        MassTraceDetection mtdet;
        mtd_param.insert("", common_param);
        mtd_param.remove("chrom_fwhm");
        mtd_param.setValue("mass_error_ppm", 10.0, "Allowed mass deviation (in ppm).");
        //mtd_param.setValue("max_trace_length", 300.0, "");
        mtdet.setParameters(mtd_param);

        mtdet.run(map, m_traces);  // m_traces : output of this function

        //-------------------------------------------------------------
        // calculate interval
        //-------------------------------------------------------------
        generateInterval(m_traces);


        //m_traces[0][0].getPosition()

        //cout << m_traces.size() << endl;

        //-------------------------------------------------------------
        // writing output
        //-------------------------------------------------------------
        // write masstrace results in file
        /*
        ofstream mfile;
        mfile.open ("/Users/jeek/Documents/A4B/matlab/masstrace.m");

        mfile << "traces={";
        for(int i=0;i<m_traces.size();i++){
            //if(m_traces[i].getMaxIntensity(false) < 1e5)continue;
            ostringstream os;
            for(int j=0;j<m_traces[i].getSize();j++) {
                os << m_traces[i][j].getRT() << "," << m_traces[i][j].getMZ() << "," << m_traces[i][j].getIntensity() << endl;
            }
            String str = os.str();
            if(str.empty()) continue;
            mfile << "[" << str <<  "];" << endl;

        }
        mfile << "};" << endl;
        mfile.close();
        */
    }

    void generateInterval(std::vector<MassTrace>& m_traces)
    {
        // input : MS1 (map), m_traces
        // output : interval set (RTs, RTe) -> m_traces?
        // param : deltaRT, N
        float rw_wndw_size = 0.01; // delatRT
        int N = 5;

        cout << m_traces[0].findMaxByIntPeak() << endl;
        cout << m_traces[0][m_traces[0].findMaxByIntPeak()] << endl;
        for(int j = 0; j<m_traces[0].getSize(); j++){
            cout << m_traces[0][j].getRT() << "," << m_traces[0][j].getMZ() << "," << m_traces[0][j].getIntensity() << endl;
        }

        // find apexes


//    for (const auto &trace : m_traces) {
//        trace.getCentroidRT()
//    }

    }
};

int main(int argc, const char** argv)
{
    TDfeatureDeconvolution tool;
    return tool.main(argc, argv);
}
