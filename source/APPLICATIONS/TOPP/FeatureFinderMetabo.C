// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Erhan Kenar $
// $Authors: Erhan Kenar, Holger Franken $
// --------------------------------------------------------------------------
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/FILTERING/DATAREDUCTION/MassTraceDetection.h>
#include <OpenMS/FILTERING/DATAREDUCTION/ElutionPeakDetection.h>
#include <OpenMS/FILTERING/DATAREDUCTION/FeatureFindingMetabo.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
        @page TOPP_FeatureFinderMetabo FeatureFinderMetabo

        @brief detects mass traces. Haha.
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFeatureFinderMetabo
        : public TOPPBase
{
public:
    TOPPFeatureFinderMetabo()
        : TOPPBase("FeatureFinderMetabo", "Detects mass traces in LC-MS data.")
    {
    }

protected:

    void registerOptionsAndFlags_()
    {
        registerInputFile_("in","<file>", "", "input centroided mzML file");
        setValidFormats_("in",StringList::create("mzML"));
        registerOutputFile_("out", "<file>", "", "output featureXML file with mass traces");
        setValidFormats_("out",StringList::create("featureXML"));

        addEmptyLine_();
        addText_("Parameters for the mass trace detection algorithm can be given in the 'algorithm' part of INI file.");
        registerSubsection_("algorithm","Algorithm parameters section");
    }

    Param getSubsectionDefaults_(const String& /*section*/) const
    {
        Param combined;
        Param p_com;
        p_com.setValue("mass_error_ppm", 20.0, "Allowed mass error deviation in ppm");
        p_com.setValue("chrom_fwhm" , 3.0 , "Lower bound for a chromatographic peak's FWHM (in seconds)");

        combined.insert("common:", p_com);

        Param p_mtd = MassTraceDetection().getDefaults();
        p_mtd.remove("mass_error_ppm");
        p_mtd.remove("chrom_fwhm");

        combined.insert("mtd:", p_mtd);

        Param p_epd = ElutionPeakDetection().getDefaults();
        p_epd.setValue("enabled", "true", "Do post-filtering of detected mass traces?");
        p_epd.setValidStrings("enabled", StringList::create("true,false"));
        combined.insert("epd:", p_epd);

        Param p_ffm = FeatureFindingMetabo().getDefaults();
        p_ffm.remove("mass_error_ppm");
        p_ffm.remove("chrom_fwhm");

        combined.insert("ffm:", p_ffm);

        return combined;
    }

    ExitCodes main_(int , const char**)
    {

        //-------------------------------------------------------------
        // parameter handling
        //-------------------------------------------------------------

        String in = getStringOption_("in");
        String out = getStringOption_("out");

        //-------------------------------------------------------------
        // loading input
        //-------------------------------------------------------------
        MzMLFile mz_data_file;
        mz_data_file.setLogType(log_type_);
        MSExperiment<Peak1D> ms_peakmap;
        mz_data_file.load(in,ms_peakmap);

        if ( ms_peakmap.empty() )
        {
            LOG_WARN << "The given file does not contain any conventional peak data, but might"
                        " contain chromatograms. This tool currently cannot handle them, sorry.";
            return INCOMPATIBLE_INPUT_DATA;
        }


        FeatureMap<> ms_feat_map;
        vector<MassTrace> m_traces;

        //-------------------------------------------------------------
        // set parameters
        //-------------------------------------------------------------

        Param common_param = getParam_().copy("algorithm:common:", true);
        writeDebug_("Common parameters passed to all sub-algorithms", common_param,3);

        Param mtdet_param = getParam_().copy("algorithm:mtd:",true);
        writeDebug_("Parameters passed to MassTraceDetection", mtdet_param,3);

        Param epd_param = getParam_().copy("algorithm:epd:",true);
        writeDebug_("Parameters passed to ElutionPeakDetection", epd_param,3);

        Param ffm_param = getParam_().copy("algorithm:ffm:", true);
        writeDebug_("Parameters passed to FeatureFindingMetabo", ffm_param,3);

        //-------------------------------------------------------------
        // configure and run mass trace detection
        //-------------------------------------------------------------

        MassTraceDetection mtdet;
        mtdet_param.insert("", common_param);
        // std::cout << "errppm ffm:" << mtdet_param.getValue("mass_error_ppm") << std::endl;
        mtdet.setParameters(mtdet_param);

        mtdet.run(ms_peakmap, m_traces);


        //-------------------------------------------------------------
        // configure and run elution peak detection
        //-------------------------------------------------------------

        //        DoubleReal fwhm(mtdet.getParameters().getValue("chrom_fwhm"));
        //        DoubleReal scan_rt_diff ((ms_peakmap[ms_peakmap.size() - 1].getRT() - ms_peakmap[0].getRT())/(ms_peakmap.size()));
        //        Size min_datapoints = std::floor(fwhm/scan_rt_diff);

        ElutionPeakDetection epdet;
        epdet.setParameters(epd_param);

        std::vector<MassTrace> splitted_mtraces;
        std::vector<MassTrace> filtered_mtraces;

        epdet.detectPeaks(m_traces, splitted_mtraces);
        epdet.filterByPeakWidth(splitted_mtraces, filtered_mtraces);


        //-------------------------------------------------------------
        // configure and run feature finding
        //-------------------------------------------------------------

        FeatureFindingMetabo ffmet;
        ffm_param.insert("", common_param);

        ffmet.setParameters(ffm_param);
        // ffmet.run(splitted_mtraces, ms_feat_map);
        ffmet.run(filtered_mtraces, ms_feat_map);



        //        for (Size i = 0; i < splitted_mtraces.size(); ++i)
        //        {
        //            MassTrace tmp_mt(splitted_mtraces[i]);

        //            Feature f;
        //            f.setMetaValue(3,tmp_mt.getLabel());
        //            f.setCharge(0);
        //            f.setMZ(tmp_mt.getCentroidMZ());
        //            f.setIntensity(tmp_mt.computePeakArea());
        //            f.setRT(tmp_mt.getSmoothedMaxRT());
        //            f.setWidth(tmp_mt.estimateFWHM());
        //            f.setOverallQuality(0.0);
        //            f.getConvexHulls().push_back(tmp_mt.getConvexhull());

        //            ms_feat_map.push_back(f);
        //        }

        //-------------------------------------------------------------
        // writing output
        //-------------------------------------------------------------

        //annotate output with data processing info TODO
        // addDataProcessing_(ms_featmap, getProcessingInfo_(DataProcessing::PEAK_PICKING));

        FeatureXMLFile().store(out, ms_feat_map);

        return EXECUTION_OK;
    }
};


int main( int argc, const char** argv )
{
    TOPPFeatureFinderMetabo tool;
    return tool.main(argc, argv);
}

/// @endcond
