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

        @brief FeatureFinderMetabo assembles metabolite features from singleton mass traces.

        <CENTER>
        <table>
        <tr>
        <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
        <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ FeatureFinderMetabo \f$ \longrightarrow \f$</td>
        <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerHiRes </td>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_TextExporter</td>
        </tr>
        <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerWavelet </td>
        </tr>
        </table>
        </CENTER>

        Mass traces alone would allow for further analyzes such as metabolite ID or statistical
        evaluation. However, in general, monoisotopic mass traces are accompanied with satellite
        C13 peaks and thus may render the analysis more difficult. @ref FeatureFinderMetabo fulfills
        a further data reduction step by assembling compatible mass traces to metabolite features
        (that is, mass traces all stemming from one metabolite). To this end, multiple metabolite
        hypotheses are formulated and scored according to how well differences in RT and m/z or
        intensity ratios match to those of theoretical isotope patterns.

        <B>The command line parameters of this tool are:</B>
        @verbinclude TOPP_FeatureFinderMetabo.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFeatureFinderMetabo
        : public TOPPBase
{
public:
    TOPPFeatureFinderMetabo()
        : TOPPBase("FeatureFinderMetabo", "Assembles metabolite features from singleton mass traces.")
    {
    }

protected:

    void registerOptionsAndFlags_()
    {
        registerInputFile_("in","<file>", "", "input centroided mzML file");
        setValidFormats_("in",StringList::create("mzML"));
        registerOutputFile_("out", "<file>", "", "output featureXML file with metabolite features");
        setValidFormats_("out",StringList::create("featureXML"));

        addEmptyLine_();
        addText_("Parameters for the mass trace detection algorithm can be given in the 'algorithm' part of INI file.");
        registerSubsection_("algorithm","Algorithm parameters section");
    }

    Param getSubsectionDefaults_(const String& /*section*/) const
    {
        Param combined;
        Param p_com;
        p_com.setValue("noise_threshold_int", 10.0, "Intensity threshold below which peaks are regarded as noise.");
        p_com.setValue("chrom_peak_snr", 3.0, "Minimum signal-to-noise a mass trace should have.");

        combined.insert("common:", p_com);

        Param p_mtd = MassTraceDetection().getDefaults();
        p_mtd.remove("noise_threshold_int");
        p_mtd.remove("chrom_peak_snr");

        combined.insert("mtd:", p_mtd);

        Param p_epd = ElutionPeakDetection().getDefaults();
        p_epd.remove("noise_threshold_int");
        p_epd.remove("chrom_peak_snr");

        // p_epd.setValue("enabled", "true", "Enables/disables the chromatographic peak detection of mass traces");
        // p_epd.setValidStrings("enabled", StringList::create("true,false"));
        combined.insert("epd:", p_epd);

        Param p_ffm = FeatureFindingMetabo().getDefaults();

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
        std::vector<Int> ms_level(1, 1);
        (mz_data_file.getOptions()).setMSLevels(ms_level);
        mz_data_file.load(in, ms_peakmap);

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
        writeDebug_("Common parameters passed to subalgorithms (mtd and ffm)", common_param,3);

        Param mtd_param = getParam_().copy("algorithm:mtd:",true);
        writeDebug_("Parameters passed to MassTraceDetection", mtd_param,3);

        Param epd_param = getParam_().copy("algorithm:epd:",true);
        writeDebug_("Parameters passed to ElutionPeakDetection", epd_param,3);

        Param ffm_param = getParam_().copy("algorithm:ffm:", true);
        writeDebug_("Parameters passed to FeatureFindingMetabo", ffm_param,3);

        //-------------------------------------------------------------
        // configure and run mass trace detection
        //-------------------------------------------------------------

        MassTraceDetection mtdet;
        mtd_param.insert("", common_param);
        mtdet.setParameters(mtd_param);

        mtdet.run(ms_peakmap, m_traces);


        //-------------------------------------------------------------
        // configure and run elution peak detection
        //-------------------------------------------------------------

        // bool use_epd = epd_param.getValue("enabled").toBool();

        std::vector<MassTrace> m_traces_final = m_traces;

        // DoubleReal pw_est(epd_param.getValue("chrom_fwhm"));
        DoubleReal scan_time(std::fabs(ms_peakmap[ms_peakmap.size() - 1].getRT() - ms_peakmap[0].getRT())/ms_peakmap.size());

        ElutionPeakDetection epdet;
        // epd_param.remove("enabled"); // artificially added above
        epd_param.insert("", common_param);
        epdet.setParameters(epd_param);

        std::vector<MassTrace> splitted_mtraces;

        epdet.setScanTime(scan_time);

        epdet.detectPeaks(m_traces, splitted_mtraces);


        if (epdet.getParameters().getValue("width_filtering").toBool())
        {
            m_traces_final.clear();
            epdet.filterByPeakWidth(splitted_mtraces, m_traces_final);
        }
        else
        {
            m_traces_final = splitted_mtraces;
        }



        //-------------------------------------------------------------
        // configure and run feature finding
        //-------------------------------------------------------------

        FeatureFindingMetabo ffmet;

        ffmet.setParameters(ffm_param);
        ffmet.run(m_traces_final, ms_feat_map);

        ms_feat_map.applyMemberFunction(&UniqueIdInterface::setUniqueId);

        //-------------------------------------------------------------
        // writing output
        //-------------------------------------------------------------

        // annotate output with data processing info
        addDataProcessing_(ms_feat_map, getProcessingInfo_(DataProcessing::QUANTITATION));

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
