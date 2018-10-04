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
        /*
         * defaults_.setValue("mass_error_ppm", 20.0, "Allowed mass deviation (in ppm).");
    defaults_.setValue("noise_threshold_int", 10.0, "Intensity threshold below which peaks are removed as noise.");
    defaults_.setValue("chrom_peak_snr", 3.0, "Minimum intensity above noise_threshold_int (signal-to-noise) a peak should have to be considered an apex.");

    defaults_.setValue("reestimate_mt_sd", "true", "Enables dynamic re-estimation of m/z variance during mass trace collection stage.");
    defaults_.setValidStrings("reestimate_mt_sd", ListUtils::create<String>("true,false"));

    defaults_.setValue("quant_method", String(MassTrace::names_of_quantmethod[0]), "Method of quantification for mass traces. For LC data 'area' is recommended, 'median' for direct injection data.");
    defaults_.setValidStrings("quant_method", std::vector<String>(MassTrace::names_of_quantmethod, MassTrace::names_of_quantmethod +(int)MassTrace::SIZE_OF_MT_QUANTMETHOD));

    // advanced parameters
    defaults_.setValue("trace_termination_criterion", "outlier", "Termination criterion for the extension of mass traces. In 'outlier' mode, trace extension cancels if a predefined number of consecutive outliers are found (see trace_termination_outliers parameter). In 'sample_rate' mode, trace extension in both directions stops if ratio of found peaks versus visited spectra falls below the 'min_sample_rate' threshold.", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("trace_termination_criterion", ListUtils::create<String>("outlier,sample_rate"));
    defaults_.setValue("trace_termination_outliers", 5, "Mass trace extension in one direction cancels if this number of consecutive spectra with no detectable peaks is reached.", ListUtils::create<String>("advanced"));

    defaults_.setValue("min_sample_rate", 0.5, "Minimum fraction of scans along the mass trace that must contain a peak.", ListUtils::create<String>("advanced"));
    defaults_.setValue("min_trace_length", 5.0, "Minimum expected length of a mass trace (in seconds).", ListUtils::create<String>("advanced"));
    defaults_.setValue("max_trace_length", -1.0, "Maximum expected length of a mass trace (in seconds). Set to a negative value to disable maximal length check during mass trace detection.", ListUtils::create<String>("advanced"));
         * */

       // MSExperiment::ConstAreaIterator start = map.areaBeginConst(10, 11, 1000,1100);
       // MSExperiment::ConstAreaIterator end = map.areaBeginConst(20, 21, 1000,1100);

        //mtdet.run(start, end, m_traces);
        mtdet.run(map, m_traces);

        ofstream mfile;
        mfile.open ("/Users/kyowonjeong/Documents/A4B/matlab/masstrace.m");
       // myfile << "Writing this to a file.\n";

        mfile << "traces={";
        for(int i=0;i<m_traces.size();i++){
            //if(m_traces[i].getMaxIntensity(false) < 1e5)continue;
            ostringstream os;
            for(int j=0;j<m_traces[i].getSize();j++) {
                //if(m_traces[i][j].getRT() > 320 || m_traces[i][j].getRT() < 280) continue;
                //if(m_traces[i][j].getMZ() > 880 || m_traces[i][j].getMZ() < 870) continue;

                os << m_traces[i][j].getRT() << "," << m_traces[i][j].getMZ() << "," << m_traces[i][j].getIntensity() << endl;
            }


            String str = os.str();
            if(str.empty()) continue;
            mfile << "[" << str <<  "];" << endl;

        }
        mfile << "};" << endl;

        mfile.close();


        //m_traces[0][0].getPosition()

        //cout << m_traces.size() << endl;

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
