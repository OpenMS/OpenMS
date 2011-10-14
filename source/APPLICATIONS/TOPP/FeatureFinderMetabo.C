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
        return MassTraceDetection().getDefaults();
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

        if (ms_peakmap.size()==0)
        {
            LOG_WARN << "The given file does not contain any conventional peak data, but might"
                    " contain chromatograms. This tool currently cannot handle them, sorry.";
            return INCOMPATIBLE_INPUT_DATA;
        }


        //-------------------------------------------------------------
        // set parameters and start extraction
        //-------------------------------------------------------------
        FeatureMap<> ms_feat_map;
        vector<MassTrace> m_traces;

        Param mt_ext_param = getParam_().copy("algorithm:",true);
        writeDebug_("Parameters passed to FeatureFinderMetabo", mt_ext_param,3);

        MassTraceDetection mt_ext;
        // mt_ext.setLogType(log_type_);
        mt_ext.setParameters(mt_ext_param);

        mt_ext.run(ms_peakmap, m_traces);

        DoubleReal fwhm(mt_ext.getParameters().getValue("chrom_fwhm"));
        DoubleReal scan_rt_diff ((ms_peakmap[ms_peakmap.size() - 1].getRT() - ms_peakmap[0].getRT())/(ms_peakmap.size()));
        Size min_datapoints = std::floor(fwhm/scan_rt_diff);

        ElutionPeakDetection ep_det;

        Param ep_det_param;
        ep_det_param.setValue("window_size", min_datapoints);
        ep_det.setParameters(ep_det_param);

        std::vector<MassTrace> splitted_mtraces;
        std::vector<MassTrace> filtered_mtraces;


        ep_det.detectPeaks(m_traces, splitted_mtraces);
        // ep_det.filterByPeakWidth(splitted_mtraces, filtered_mtraces);

        FeatureFindingMetabo ff_met;

        ff_met.run(splitted_mtraces, ms_feat_map);
        // ff_met.run(filtered_mtraces, ms_feat_map);



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
