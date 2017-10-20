// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Erhan Kenar, Holger Franken $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FILTERING/DATAREDUCTION/MassTraceDetection.h>
#include <OpenMS/FILTERING/DATAREDUCTION/ElutionPeakDetection.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/FORMAT/FileHandler.h>

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

class TOPPMassTraceExtractor :
  public TOPPBase
{
public:
  TOPPMassTraceExtractor() :
    TOPPBase("MassTraceExtractor", "Detects mass traces in centroided LC-MS data.")
  {
  }

protected:

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "input centroided mzML file");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "output featureXML file with mass traces");
    setValidFormats_("out", ListUtils::create<String>("featureXML,consensusXML"));
    registerStringOption_("out_type", "<type>", "", "output file type -- default: determined from file extension or content", false);
    setValidStrings_("out_type", ListUtils::create<String>("featureXML,consensusXML"));

    addEmptyLine_();
    registerSubsection_("algorithm", "Algorithm parameters section");

  }

  Param getSubsectionDefaults_(const String& /*section*/) const
  {
    Param combined;
    Param p_com;
    p_com.setValue("noise_threshold_int", 10.0, "Intensity threshold below which peaks are regarded as noise.");
    p_com.setValue("chrom_peak_snr", 3.0, "Minimum signal-to-noise a mass trace should have.");
    p_com.setValue("chrom_fwhm", 5.0, "Expected chromatographic peak width (in seconds).");

    combined.insert("common:", p_com);

    Param p_mtd = MassTraceDetection().getDefaults();
    p_mtd.remove("noise_threshold_int");
    p_mtd.remove("chrom_peak_snr");

    combined.insert("mtd:", p_mtd);

    Param p_epd = ElutionPeakDetection().getDefaults();
    p_epd.remove("noise_threshold_int");
    p_epd.remove("chrom_peak_snr");
    p_epd.remove("chrom_fwhm");

    p_epd.setValue("enabled", "true", "Enables/disables the chromatographic peak detection of mass traces");
    p_epd.setValidStrings("enabled", ListUtils::create<String>("true,false"));
    combined.insert("epd:", p_epd);

    return combined;
  }

  ExitCodes main_(int, const char**)
  {

    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    String in = getStringOption_("in");
    String out = getStringOption_("out");
    FileTypes::Type out_type = FileTypes::nameToType(getStringOption_("out_type"));

    if (out_type == FileTypes::UNKNOWN)
    {
      out_type = FileHandler().getTypeByFileName(out);
    }

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------
    MzMLFile mz_data_file;
    mz_data_file.setLogType(log_type_);
    PeakMap ms_peakmap;
    std::vector<Int> ms_level(1, 1);
    (mz_data_file.getOptions()).setMSLevels(ms_level);
    mz_data_file.load(in, ms_peakmap);

    if (ms_peakmap.size() == 0)
    {
      LOG_WARN << "The given file does not contain any conventional peak data, but might"
                  " contain chromatograms. This tool currently cannot handle them, sorry.";
      return INCOMPATIBLE_INPUT_DATA;
    }

    // make sure that the spectra are sorted by m/z
    ms_peakmap.sortSpectra(true);

    //-------------------------------------------------------------
    // get params for MTD and EPD algorithms
    //-------------------------------------------------------------
    Param com_param = getParam_().copy("algorithm:common:", true);
    writeDebug_("Common parameters passed to both sub-algorithms (mtd and epd)", com_param, 3);

    Param mtd_param = getParam_().copy("algorithm:mtd:", true);
    writeDebug_("Parameters passed to MassTraceDetection", mtd_param, 3);

    Param epd_param = getParam_().copy("algorithm:epd:", true);
    writeDebug_("Parameters passed to ElutionPeakDetection", epd_param, 3);


    //-------------------------------------------------------------
    // configure and run MTD
    //-------------------------------------------------------------

    MassTraceDetection mt_ext;
    mtd_param.insert("", com_param);
    mtd_param.remove("chrom_fwhm");
    mt_ext.setParameters(mtd_param);
    vector<MassTrace> m_traces;
    mt_ext.run(ms_peakmap, m_traces);

    vector<MassTrace> m_traces_final;

    bool use_epd = epd_param.getValue("enabled").toBool();

    if (!use_epd)
    {
      swap(m_traces_final, m_traces);
    }
    else
    {
      ElutionPeakDetection ep_det;

      epd_param.remove("enabled"); // artificially added above
      epd_param.insert("", com_param);

      ep_det.setParameters(epd_param);

      std::vector<MassTrace> split_mtraces;
      // note: this step will destroy any meta data annotation (e.g. FWHM_mz_avg)
      ep_det.detectPeaks(m_traces, split_mtraces);

      if (ep_det.getParameters().getValue("width_filtering") == "auto")
      {
        m_traces_final.clear();
        ep_det.filterByPeakWidth(split_mtraces, m_traces_final);

        LOG_INFO << "Notice: " << split_mtraces.size() - m_traces_final.size()
                 << " of total " << split_mtraces.size() 
                 << " were dropped because of too low peak width." << std::endl;
      }
      else
      {
        swap(m_traces_final, split_mtraces);
      }
    }

    //-------------------------------------------------------------
    // writing consensus map output
    //-------------------------------------------------------------
    if (out_type == FileTypes::CONSENSUSXML)
    {
      ConsensusMap consensus_map;
      StringList ms_runs;
      ms_peakmap.getPrimaryMSRunPath(ms_runs);
      consensus_map.setPrimaryMSRunPath(ms_runs);

      for (Size i = 0; i < m_traces_final.size(); ++i)
      {
        if (m_traces_final[i].getSize() == 0) continue;

        ConsensusFeature fcons;
        int k = 0;
        for (MassTrace::const_iterator it = m_traces_final[i].begin(); it != m_traces_final[i].end(); ++it)
        {
          FeatureHandle fhandle;
          fhandle.setRT(it->getRT());
          fhandle.setMZ(it->getMZ());
          fhandle.setIntensity(it->getIntensity());
          fhandle.setUniqueId(++k);
          fcons.insert(fhandle);
        }

        fcons.setMetaValue(3, m_traces_final[i].getLabel());
        fcons.setCharge(0);
        fcons.setWidth(m_traces_final[i].estimateFWHM(use_epd));
        fcons.setQuality(1 - (1.0 / m_traces_final[i].getSize()));

        fcons.setRT(m_traces_final[i].getCentroidRT());
        fcons.setMZ(m_traces_final[i].getCentroidMZ());
        fcons.setIntensity(m_traces_final[i].getIntensity(false));
        consensus_map.push_back(fcons);
      }
      consensus_map.applyMemberFunction(&UniqueIdInterface::setUniqueId);
      addDataProcessing_(consensus_map, getProcessingInfo_(DataProcessing::QUANTITATION));
      consensus_map.setUniqueId();
      ConsensusXMLFile().store(out, consensus_map);

    }
    else //(out_type == FileTypes::FEATUREXML)
    {

      //-----------------------------------------------------------
      // convert mass traces to features
      //-----------------------------------------------------------

      std::vector<double> stats_sd;
      FeatureMap ms_feat_map;
      StringList ms_runs;
      ms_peakmap.getPrimaryMSRunPath(ms_runs);
      ms_feat_map.setPrimaryMSRunPath(ms_runs);
      for (Size i = 0; i < m_traces_final.size(); ++i)
      {
        if (m_traces_final[i].getSize() == 0) continue;

        m_traces_final[i].updateMeanMZ();
        m_traces_final[i].updateWeightedMZsd();

        Feature f;
        f.setMetaValue(3, m_traces_final[i].getLabel());
        f.setCharge(0);
        f.setMZ(m_traces_final[i].getCentroidMZ());
        f.setIntensity(m_traces_final[i].getIntensity(false));
        f.setRT(m_traces_final[i].getCentroidRT());
        f.setWidth(m_traces_final[i].estimateFWHM(use_epd));
        f.setOverallQuality(1 - (1.0 / m_traces_final[i].getSize()));
        f.getConvexHulls().push_back(m_traces_final[i].getConvexhull());
        double sd = m_traces_final[i].getCentroidSD();
        f.setMetaValue("SD", sd);
        f.setMetaValue("SD_ppm", sd / f.getMZ() * 1e6);
        if (m_traces_final[i].fwhm_mz_avg > 0) f.setMetaValue("FWHM_mz_avg", m_traces_final[i].fwhm_mz_avg);
        stats_sd.push_back(m_traces_final[i].getCentroidSD());
        ms_feat_map.push_back(f);
      }

      // print some stats about standard deviation of mass traces
      if (stats_sd.size() > 0)
      {
        std::sort(stats_sd.begin(), stats_sd.end());
        LOG_INFO << "Mass trace m/z s.d.\n"
                 << "    low quartile: " << stats_sd[stats_sd.size() * 1 / 4] << "\n"
                 << "          median: " << stats_sd[stats_sd.size() * 1 / 2] << "\n"
                 << "    upp quartile: " << stats_sd[stats_sd.size() * 3 / 4] << std::endl;
      }


      ms_feat_map.applyMemberFunction(&UniqueIdInterface::setUniqueId);

      //-------------------------------------------------------------
      // writing output
      //-------------------------------------------------------------

      // annotate output with data processing info TODO
      addDataProcessing_(ms_feat_map, getProcessingInfo_(DataProcessing::QUANTITATION));
      //ms_feat_map.setUniqueId();

      FeatureXMLFile().store(out, ms_feat_map);
    }

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPMassTraceExtractor tool;
  return tool.main(argc, argv);
}

/// @endcond
