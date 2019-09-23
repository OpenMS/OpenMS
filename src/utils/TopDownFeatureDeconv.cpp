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

    typedef struct MassTraceApex{
        double MZ, RT;
        double FWHM;
        Size mt_idx, et_idx;

        MassTraceApex(): MZ(0.0), RT(0.0), FWHM(0.0), mt_idx(0), et_idx(0){}

        static double cluster_similartity(const MassTraceApex& a, const MassTraceApex& b){
            return abs( a.MZ-b.MZ ); // just simple example
        }
    }MTApex;

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
        Param common_param = getParam_().copy("algorithm:common:", true);
        writeDebug_("Common parameters passed to sub-algorithms (mtd and ffm)", common_param, 3);

        Param mtd_param = getParam_().copy("algorithm:mtd:", true);
        writeDebug_("Parameters passed to MassTraceDetection", mtd_param, 3);

        Param epd_param = getParam_().copy("algorithm:epd:", true);
        writeDebug_("Parameters passed to ElutionPeakDetection", epd_param, 3);

        //-------------------------------------------------------------
        // configure and run mass trace detection
        //-------------------------------------------------------------
        vector<MassTrace> m_traces;
        MassTraceDetection mtdet;
        mtd_param.insert("", common_param);
        mtd_param.remove("chrom_fwhm");
        mtd_param.setValue("mass_error_ppm", 10.0, "Allowed mass deviation (in ppm).");
        //mtd_param.setValue("max_trace_length", 300.0, "");
        mtdet.setParameters(mtd_param);

        mtdet.run(map, m_traces);  // m_traces : output of this function

        //-------------------------------------------------------------
        // configure and run elution peak detection (find apexes)
        //-------------------------------------------------------------
        std::vector<vector<MassTrace>> elutionPeaks;
        epd_param.remove("enabled"); // artificially added above
        epd_param.insert("", common_param);
        ElutionPeakDetection epdet;
        epdet.setParameters(epd_param);
        // fill mass traces with smoothed data as well --> code exists in FeatureFinderMetabo.cpp (not using here)
        // some filtering should work here, but removed for now. (22.10.18)

        vector<MTApex> apexes; // pair of elution peaks and original mass trace
        detectPeaks(epdet, m_traces, elutionPeaks, apexes);

        //-------------------------------------------------------------
        // generate elution peak clusters (close_apexes)
        //-------------------------------------------------------------
        // sort apexes w/ RT
        sort(apexes.begin(), apexes.end(), [this](MTApex &a, MTApex &b){
            return this->sortApexes_byRT(a,b);
        });
        std::vector<vector<MTApex>> elu_clusters;
        generateElutionPeakClusters(apexes, elu_clusters);

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

    void print_apexes(vector<MTApex> &m)
    {
        cout << "m_starts : " << m.size() << endl;
        for (auto&x : m) {
            cout << x.RT << "," << x.MZ << "," << x.FWHM << endl;
        }
    }

    bool sortApexes_byRT(const MTApex &a, const MTApex &b)
    {
        return ( a.RT < b.RT);
    }

    void generateElutionPeakClusters(std::vector<MTApex>& apexes, std::vector<vector<MTApex>> &elucluster)
    {
        // input : the set of apexes
        // output : the set of elution peak clusters
        // param : deltaRT, N

        elucluster.clear();

        double rw_wndw_size = 0.1; // delatRT
        unsigned long N = 6;

        std::vector<MTApex> closeApexes;
        int apexesSize = apexes.size();
        int pntr=0; // end of the window
        for (int i=0;  i<apexesSize; i++) { // i : start of the window
//            cout << i << " : " << apexes[i].RT << "," << apexes[i].MZ << endl;
            if(closeApexes.empty()) closeApexes.push_back(apexes[i]);

            // collect closeApexes
            if( pntr == i ) pntr++;
            double stdRT = apexes[i].RT;
            while( pntr<apexesSize ){
                if( apexes[pntr].RT-stdRT > rw_wndw_size ) break;

                closeApexes.push_back(apexes[pntr]);
                pntr++;
            }

            if( closeApexes.size() < N ) {
                closeApexes.erase(closeApexes.begin());
                continue;
            }
            // otherwise, calculate cluster here
            clusteringByFWHM(closeApexes, elucluster);

            //
//            print_apexes(closeApexes); // for testing
            closeApexes.erase(closeApexes.begin()); // remove i from window after clustering
        }
    }

    void clusteringByFWHM(std::vector<MTApex> &closeAps, std::vector<vector<MTApex>> &elu_clus){
        for(auto& x: closeAps){
        }
    }

    // this is from 'elutionPeakDetection'
    void detectPeaks(ElutionPeakDetection& epdet, vector<MassTrace>& mt_vec,vector<vector<MassTrace>>& elu_vec, vector<MTApex>& apexes)
    {
        // make sure all things are cleared
        elu_vec.clear();
        apexes.clear();
        vector<MassTrace> mt_vec_tmp; // copy to remember the reference

        cout << "elution peak detection" << endl;

        SignedSize mt_size = mt_vec.size();
        for (SignedSize i = 0; i < mt_size; ++i) {

            // push_back to 'single_mtraces' is protected, so threading is ok
            vector<MassTrace> single_mtraces;
            single_mtraces.clear();
            epdet.detectPeaks(mt_vec[i], single_mtraces);
            if (!single_mtraces.empty()) {
                elu_vec.push_back(single_mtraces);
                mt_vec_tmp.push_back(mt_vec[i]);
                for (Size j = 0; j < single_mtraces.size(); j++) {
                    MTApex mta = MTApex();
                    PeakType p = single_mtraces[j][single_mtraces[j].findMaxByIntPeak()];
                    mta.MZ = p.getMZ();
                    mta.RT = p.getRT();
                    mta.et_idx = j;
                    mta.mt_idx = elu_vec.size() - 1;
                    mta.FWHM = single_mtraces[j].getFWHM();
                    apexes.push_back(mta);
                }
            }
        }
        // no apex -> remove mass trace from vector
        mt_vec.clear();
        copy(mt_vec_tmp.begin(), mt_vec_tmp.end(), back_inserter(mt_vec));
    }

};

int main(int argc, const char** argv)
{
    TDfeatureDeconvolution tool;
    return tool.main(argc, argv);
}
