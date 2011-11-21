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
// #include <OpenMS/FILTERING/TRANSFORMERS/TICResampling.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
        @page TOPP_MassTraceExtractor MassTraceExtractor

        @brief detects mass traces.

        <CENTER>
        <table>
        <tr>
        <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
        <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ MassTraceExtractor \f$ \longrightarrow \f$</td>
        <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerHiRes </td>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_Decharger</td>
        </tr>
        <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerWavelet </td>
        </tr>
        </table>
        </CENTER>


        Annotates mass traces in centroided LC/MS maps.
        Useful for metabolomics and top-down proteomics. 
        Use FeatureFinder<xxx> tools for peptide data.
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPMassTraceExtractor
    : public TOPPBase
{
public:
    TOPPMassTraceExtractor()
        : TOPPBase("MassTraceExtractor", "Detects mass traces in LC-MS data.")
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
        combined.insert("mtd:", MassTraceDetection().getDefaults());
        Param p_epd = ElutionPeakDetection().getDefaults();
        p_epd.setValue("enabled", "true", "Do post-filtering of detected mass traces?!");
        p_epd.setValidStrings("enabled", StringList::create("true,false"));
        combined.insert("epd:", p_epd);

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

        if (ms_peakmap.size()==0)
        {
            LOG_WARN << "The given file does not contain any conventional peak data, but might"
                    " contain chromatograms. This tool currently cannot handle them, sorry.";
            return INCOMPATIBLE_INPUT_DATA;
        }


        //                MSExperiment<Peak1D> testmap;

        //                TICResampling tic_res;

        //                tic_res.run(ms_peakmap, testmap);

        //                MzMLFile().store("resampled.mzML", testmap);

        //                return EXECUTION_OK;


        //-------------------------------------------------------------
        // set parameters and start extraction
        //-------------------------------------------------------------
        FeatureMap<> ms_feat_map;
        vector<MassTrace> m_traces;

        Param mt_ext_param = getParam_().copy("algorithm:mtd:",true);
        writeDebug_("Parameters passed to MassTraceDetection", mt_ext_param,3);

        Param epd_param = getParam_().copy("algorithm:epd:",true);
        writeDebug_("Parameters passed to ElutionPeakDetection", epd_param,3);


        MassTraceDetection mt_ext;
        // mt_ext.setLogType(log_type_);
        mt_ext.setParameters(mt_ext_param);

        mt_ext.run(ms_peakmap, m_traces);
        
        vector<MassTrace> m_traces_final = m_traces;

        bool use_epd = epd_param.getValue("enabled") == "true";

        if (use_epd)
        {
            DoubleReal fwhm(mt_ext.getParameters().getValue("chrom_fwhm"));
            DoubleReal scan_rt_diff ((ms_peakmap[ms_peakmap.size() - 1].getRT() - ms_peakmap[0].getRT())/(ms_peakmap.size()));
            Size min_datapoints = std::floor(fwhm/scan_rt_diff);

            ElutionPeakDetection ep_det;

            epd_param.remove("enabled"); // artificially added above
            epd_param.setValue("window_size", min_datapoints);

            ep_det.setParameters(epd_param);

            std::vector<MassTrace> splitted_mtraces;

            std::vector<MassTrace> filtered_mtraces;

            ep_det.detectPeaks(m_traces, splitted_mtraces);

            if (ep_det.getParameters().getValue("width_filtering") == "true")
            {
                ep_det.filterByPeakWidth(splitted_mtraces, filtered_mtraces);

                LOG_INFO << "After filtering: " << filtered_mtraces.size() << " of " << splitted_mtraces.size() << std::endl;

                splitted_mtraces = filtered_mtraces;
            }

            m_traces_final = splitted_mtraces;
        }

        //-----------------------------------------------------------
        // convert mass traces to features
        //-----------------------------------------------------------

        std::map<DoubleReal, map<DoubleReal, DoubleReal> > out_map;

        for (Size i = 0; i < m_traces_final.size(); ++i)
        {
            MassTrace tmp_mt(m_traces_final[i]);
            if (tmp_mt.getSize() == 0) continue;

            Feature f;
            f.setMetaValue(3,tmp_mt.getLabel());
            f.setCharge(0);
            f.setMZ(tmp_mt.getCentroidMZ());
            f.setIntensity(tmp_mt.computePeakArea());
            //            if (use_epd)
            //            {
            //              f.setRT(tmp_mt.getSmoothedMaxRT());
            //              f.setWidth(tmp_mt.estimateFWHM());
            //            }
            //            else
            //            {
            //              f.setRT(tmp_mt.getCentroidRT());
            //            }
            f.setRT(tmp_mt.getCentroidRT());
            f.setWidth(tmp_mt.estimateFWHM());

            f.setOverallQuality(1 - (1.0/tmp_mt.getSize()));
            f.getConvexHulls().push_back(tmp_mt.getConvexhull());

            Size mtr_idx(0);

            for (MassTrace::const_iterator c_it = m_traces_final[i].begin(); c_it != m_traces_final[i].end(); ++c_it)
            {
                DoubleReal p_int(m_traces_final[i].getSmoothedIntensities()[mtr_idx]);
                if (p_int > 0.0) {
                    out_map[c_it->getRT()][c_it->getMZ()] = p_int;
                }
                ++mtr_idx;
            }


            ms_feat_map.push_back(f);
        }
        //-------------------------------------------------------------
        // writing output
        //-------------------------------------------------------------

        MSExperiment<Peak1D> out_msexp;

        for (std::map<DoubleReal, map<DoubleReal, DoubleReal> >::const_iterator m_it = out_map.begin(); m_it != out_map.end(); ++m_it)
        {
            MSSpectrum<Peak1D> tmp_spec;
            tmp_spec.setRT(m_it->first);

            for (std::map<DoubleReal, DoubleReal>::const_iterator c_it = m_it->second.begin(); c_it != m_it->second.end(); ++c_it)
            {
                Peak1D tmp_peak;
                tmp_peak.setMZ(c_it->first);
                tmp_peak.setIntensity(c_it->second);
                tmp_spec.push_back(tmp_peak);
            }

            out_msexp.push_back(tmp_spec);
        }


        MzMLFile().store("raw_out.mzML", out_msexp);

        //annotate output with data processing info TODO
        // addDataProcessing_(ms_featmap, getProcessingInfo_(DataProcessing::PEAK_PICKING));

        FeatureXMLFile().store(out, ms_feat_map);

        return EXECUTION_OK;
    }
};


int main( int argc, const char** argv )
{
    TOPPMassTraceExtractor tool;
    return tool.main(argc, argv);
}

/// @endcond
