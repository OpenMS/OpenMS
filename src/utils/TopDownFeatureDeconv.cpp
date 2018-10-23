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
#include <OpenMS/FILTERING/DATAREDUCTION/ElutionPeakDetection.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace OpenMS;
using namespace std;

class TDfeatureDeconvolution :
    public TOPPBase, public ProgressLogger
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
        // set parameters
        //-------------------------------------------------------------
        vector<MassTrace> m_traces;
        Param common_param = getParam_().copy("algorithm:common:", true);
        writeDebug_("Common parameters passed to sub-algorithms (mtd and ffm)", common_param, 3);

        Param mtd_param = getParam_().copy("algorithm:mtd:", true);
        writeDebug_("Parameters passed to MassTraceDetection", mtd_param, 3);

        Param epd_param = getParam_().copy("algorithm:epd:", true);
        writeDebug_("Parameters passed to ElutionPeakDetection", epd_param, 3);

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
        // configure and run elution peak detection
        //-------------------------------------------------------------

//        std::vector<MassTrace> m_traces_final;

        std::vector<MassTrace> splitted_mtraces;
        epd_param.remove("enabled"); // artificially added above
        epd_param.insert("", common_param);
        ElutionPeakDetection epdet;
        epdet.setParameters(epd_param);
        // fill mass traces with smoothed data as well --> code exists in FeatureFinderMetabo.cpp (not using here)
        // some filtering should work here, but removed for now. (22.10.18)

        vector<pair<MassTrace, vector<MassTrace>>> mt_submt;
        detectPeaks(epdet, m_traces, mt_submt);

        //
        cout << "mass traces" << endl;
        cout << m_traces.size() << endl;
        cout << m_traces[0][m_traces[0].findMaxByIntPeak()] << endl;
        int cnt = 0;
        for (auto &m_trace : m_traces) {
            cnt += m_trace.getSize();
        }
        cout << cnt << endl;

        cout << "elution detection" << endl;
        cout << mt_submt.size() << endl;\
        cout << mt_submt[0].first[mt_submt[0].first.findMaxByIntPeak()] << endl;
        cout << mt_submt[0].second.size() << endl;

        cnt = 0;
        vector<int> maxINmt;
        for (auto &i : mt_submt) {
            maxINmt.push_back(i.second.size());
            for (auto &j : i.second) {
                cnt += j.getSize();
            }
        }
        cout << (1.0 * accumulate(maxINmt.begin(), maxINmt.end(), 0LL) / maxINmt.size()) << endl;
        cout << cnt << endl;



        //-------------------------------------------------------------
        // calculate interval
        //-------------------------------------------------------------
//        generateInterval(m_traces);


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

        return EXECUTION_OK;
    }

    void generateInterval(std::vector<MassTrace>& m_traces)
    {
        // input : MS1 (map), m_traces
        // output : interval set (RTs, RTe) -> m_traces?
        // param : deltaRT, N
        float rw_wndw_size = 0.01; // delatRT
        int N = 5;

//        cout << m_traces[0].findMaxByIntPeak() << endl;
//        cout << m_traces[0][m_traces[0].findMaxByIntPeak()] << endl;
//        for(int j = 0; j<m_traces[0].getSize(); j++){
//            cout << m_traces[0][j].getRT() << "," << m_traces[0][j].getMZ() << "," << m_traces[0][j].getIntensity() << endl;
//        }

        // find apexes


//    for (const auto &trace : m_traces) {
//        trace.getCentroidRT()
//    }

    }


    // this is from 'elutionPeakDetection'
    void detectPeaks(ElutionPeakDetection& epdet, vector<MassTrace>& mt_vec, vector<pair<MassTrace, vector<MassTrace>>>& out_mtraces)
    {
        out_mtraces.clear();
        vector<MassTrace> single_mtraces;

        this->startProgress(0, mt_vec.size(), "elution peak detection");
        Size progress(0);
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif

        SignedSize mt_size = mt_vec.size();
        for (SignedSize i = 0; i < mt_size; ++i)
        {
            IF_MASTERTHREAD this->setProgress(progress);

    #ifdef _OPENMP
    #pragma omp atomic
    #endif
            ++progress;

            // push_back to 'single_mtraces' is protected, so threading is ok
            epdet.detectPeaks(mt_vec[i], single_mtraces);
            if( !single_mtraces.empty() )
                out_mtraces.emplace_back( make_pair(mt_vec[i], single_mtraces) );
        }

        this->endProgress();
    }

};

int main(int argc, const char** argv)
{
    TDfeatureDeconvolution tool;
    return tool.main(argc, argv);
}
