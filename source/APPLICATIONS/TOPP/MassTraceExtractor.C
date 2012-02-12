// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/FORMAT/PeakFileOptions.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/FILTERING/DATAREDUCTION/MassTraceDetection.h>
#include <OpenMS/FILTERING/DATAREDUCTION/ElutionPeakDetection.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
        @page TOPP_MassTraceExtractor MassTraceExtractor

        @brief MassTraceExtractor extracts mass traces from a @ref MSExperiment map and stores them into a @ref FeatureXMLFile.

        <CENTER>
        <table>
        <tr>
        <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
        <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ MassTraceExtractor \f$ \longrightarrow \f$</td>
        <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerHiRes </td>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureFinderMetabo</td>
        </tr>
        <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerWavelet </td>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_TextExporter </td>
        </tr>
        </table>
        </CENTER>


        This TOPP tool detects mass traces in centroided LC-MS maps and stores them as features in
        a @ref FeatureMap. These features may be either used directly as input for an metabolite ID approach or further
        be assembled to aggregate features according to a theoretical isotope pattern. For metabolomics experiments,
        the @ref TOPP_FeatureFinderMetabo tool offers both mass trace extraction and isotope pattern assembly.
        For proteomics data, please refer to the @ref TOPP_FeatureFinderCentroided tool.

        <B>The command line parameters of this tool are:</B>
        @verbinclude TOPP_MassTraceExtractor.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPMassTraceExtractor
    : public TOPPBase
{
public:
    TOPPMassTraceExtractor()
        : TOPPBase("MassTraceExtractor", "Detects mass traces in centroided LC-MS data.")
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
        p_com.setValue("chrom_fwhm" , 0.0 , "Allows filtering of mass traces with peak width (in seconds) less than this threshold. Disabled by default (set to 0.0).");
        combined.insert("common:", p_com);

        Param p_mtd = MassTraceDetection().getDefaults();
        p_mtd.remove("chrom_fwhm");
        combined.insert("mtd:", p_mtd);

        Param p_epd = ElutionPeakDetection().getDefaults();
        p_epd.remove("chrom_fwhm");
        p_epd.setValue("enabled", "true", "Switches on/off the detection of elution peaks");
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
        std::vector<Int> ms_level(1, 1);
        (mz_data_file.getOptions()).setMSLevels(ms_level);
        mz_data_file.load(in,ms_peakmap);

        if (ms_peakmap.size()==0)
        {
            LOG_WARN << "The given file does not contain any conventional peak data, but might"
                    " contain chromatograms. This tool currently cannot handle them, sorry.";
            return INCOMPATIBLE_INPUT_DATA;
        }


        FeatureMap<> ms_feat_map;
        vector<MassTrace> m_traces;

        //-------------------------------------------------------------
        // get params for MTD and EPD algorithms
        //-------------------------------------------------------------
        Param com_param = getParam_().copy("algorithm:common:",true);
        writeDebug_("Common parameters passed to both subalgorithms (mtd and epd)", com_param,3);

        Param mtd_param = getParam_().copy("algorithm:mtd:",true);
        writeDebug_("Parameters passed to MassTraceDetection", mtd_param,3);

        Param epd_param = getParam_().copy("algorithm:epd:",true);
        writeDebug_("Parameters passed to ElutionPeakDetection", epd_param,3);


        //-------------------------------------------------------------
        // configure and run MTD
        //-------------------------------------------------------------

        MassTraceDetection mt_ext;
        mtd_param.insert("", com_param);
        mt_ext.setParameters(mtd_param);
        mt_ext.run(ms_peakmap, m_traces);
        
        vector<MassTrace> m_traces_final = m_traces;

        bool use_epd = epd_param.getValue("enabled").toBool();

        if (use_epd)
        {
            ElutionPeakDetection ep_det;

            epd_param.remove("enabled"); // artificially added above
            epd_param.insert("", com_param);

            ep_det.setParameters(epd_param);

            std::vector<MassTrace> splitted_mtraces;

            ep_det.detectPeaks(m_traces, splitted_mtraces);

            if (ep_det.getParameters().getValue("width_filtering").toBool())
            {
                m_traces_final.clear();
                ep_det.filterByPeakWidth(splitted_mtraces, m_traces_final);

                LOG_INFO << "Notice: " << splitted_mtraces.size() - m_traces_final.size() << " of total " << splitted_mtraces.size() << " were dropped because of too low peak width." << std::endl;
            }
            else
            {
            m_traces_final = splitted_mtraces;
        }
        }

        //-----------------------------------------------------------
        // convert mass traces to features
        //-----------------------------------------------------------

        for (Size i = 0; i < m_traces_final.size(); ++i)
        {
            if (m_traces_final[i].getSize() == 0) continue;

            Feature f;
            f.setMetaValue(3, m_traces_final[i].getLabel());
            f.setCharge(0);
            f.setMZ(m_traces_final[i].getCentroidMZ());
            f.setIntensity(m_traces_final[i].computePeakArea());
            f.setRT(m_traces_final[i].getCentroidRT());
            f.setWidth(m_traces_final[i].estimateFWHM(use_epd));
            f.setOverallQuality(1 - (1.0/m_traces_final[i].getSize()));
            f.getConvexHulls().push_back(m_traces_final[i].getConvexhull());

            ms_feat_map.push_back(f);
                }

        ms_feat_map.applyMemberFunction(&UniqueIdInterface::setUniqueId);

        //-------------------------------------------------------------
        // writing output
        //-------------------------------------------------------------

        //annotate output with data processing info TODO
        addDataProcessing_(ms_feat_map, getProcessingInfo_(DataProcessing::QUANTITATION));
        ms_feat_map.setUniqueId();

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
