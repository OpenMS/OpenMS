// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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

        If the raw data scans contain the scan polarity information, it is stored as meta value "scan_polarity" in the output file.

        <B>The command line parameters of this tool are:</B>
        @verbinclude TOPP_FeatureFinderMetabo.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFeatureFinderMetabo :
  public TOPPBase
{
public:
  TOPPFeatureFinderMetabo() :
    TOPPBase("FeatureFinderMetabo", "Assembles metabolite features from singleton mass traces.")
  {
  }

protected:

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "input centroided mzML file");
    setValidFormats_("in", StringList::create("mzML"));
    registerOutputFile_("out", "<file>", "", "output featureXML file with metabolite features");
    setValidFormats_("out", StringList::create("featureXML"));

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

    // p_epd.setValue("enabled", "true", "Enables/disables the chromatographic peak detection of mass traces");
    // p_epd.setValidStrings("enabled", StringList::create("true,false"));
    combined.insert("epd:", p_epd);

    Param p_ffm = FeatureFindingMetabo().getDefaults();
    p_ffm.remove("chrom_fwhm");

    combined.insert("ffm:", p_ffm);

    return combined;
  }

  ExitCodes main_(int, const char**)
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

    if (ms_peakmap.empty())
    {
      LOG_WARN << "The given file does not contain any conventional peak data, but might"
                  " contain chromatograms. This tool currently cannot handle them, sorry.";
      return INCOMPATIBLE_INPUT_DATA;
    }

    // make sure the spectra are sorted by m/z
    ms_peakmap.sortSpectra(true);

    vector<MassTrace> m_traces;

    //-------------------------------------------------------------
    // set parameters
    //-------------------------------------------------------------

    Param common_param = getParam_().copy("algorithm:common:", true);
    writeDebug_("Common parameters passed to subalgorithms (mtd and ffm)", common_param, 3);

    Param mtd_param = getParam_().copy("algorithm:mtd:", true);
    writeDebug_("Parameters passed to MassTraceDetection", mtd_param, 3);

    Param epd_param = getParam_().copy("algorithm:epd:", true);
    writeDebug_("Parameters passed to ElutionPeakDetection", epd_param, 3);

    Param ffm_param = getParam_().copy("algorithm:ffm:", true);
    writeDebug_("Parameters passed to FeatureFindingMetabo", ffm_param, 3);

    //-------------------------------------------------------------
    // configure and run mass trace detection
    //-------------------------------------------------------------

    MassTraceDetection mtdet;
    mtd_param.insert("", common_param);
    mtd_param.remove("chrom_fwhm");
    mtdet.setParameters(mtd_param);

    mtdet.run(ms_peakmap, m_traces);


    //-------------------------------------------------------------
    // configure and run elution peak detection
    //-------------------------------------------------------------

    // bool use_epd = epd_param.getValue("enabled").toBool();

    std::vector<MassTrace> m_traces_final = m_traces;

    // DoubleReal pw_est(epd_param.getValue("chrom_fwhm"));
    // DoubleReal scan_time(std::fabs(ms_peakmap[ms_peakmap.size() - 1].getRT() - ms_peakmap[0].getRT()) / ms_peakmap.size());

    ElutionPeakDetection epdet;
    // epd_param.remove("enabled"); // artificially added above
    epd_param.insert("", common_param);
    epdet.setParameters(epd_param);

    std::vector<MassTrace> splitted_mtraces;

    // epdet.setScanTime(scan_time);

    epdet.detectPeaks(m_traces, splitted_mtraces);


    if (epdet.getParameters().getValue("width_filtering") == "auto")
    {
      m_traces_final.clear();
      epdet.filterByPeakWidth(splitted_mtraces, m_traces_final);
    }
    else
    {
      m_traces_final = splitted_mtraces;
    }

//    std::cout << "m_traces: " << m_traces_final.size() << std::endl;

    //-------------------------------------------------------------
    // configure and run feature finding
    //-------------------------------------------------------------

    FeatureFindingMetabo ffmet;
    ffm_param.insert("", common_param);
    ffm_param.remove("noise_threshold_int");
    ffm_param.remove("chrom_peak_snr");

    FeatureMap<> ms_feat_map;
    ffmet.setParameters(ffm_param);
    ffmet.run(m_traces_final, ms_feat_map);

    ms_feat_map.sortByMZ();
    ms_feat_map.applyMemberFunction(&UniqueIdInterface::setUniqueId);

    // store ionization mode of spectra (useful for postprocessing by AccurateMassSearch tool)
    if (ms_feat_map.size() > 0)
    {
      set<IonSource::Polarity> pols;
      for (Size i=0; i < ms_peakmap.size(); ++i)
      {
        pols.insert(ms_peakmap[i].getInstrumentSettings().getPolarity());
      }
      // concat to single string
      StringList sl_pols;
      for (set<IonSource::Polarity>::const_iterator it = pols.begin();
                                                    it !=pols.end();
                                                    ++it)
      {
        sl_pols.push_back(String(IonSource::NamesOfPolarity[*it]));
      }
      ms_feat_map[0].setMetaValue("scan_polarity", sl_pols.concatenate(";"));
    }


    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    // annotate output with data processing info
    addDataProcessing_(ms_feat_map, getProcessingInfo_(DataProcessing::QUANTITATION));

    FeatureXMLFile().store(out, ms_feat_map);

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPFeatureFinderMetabo tool;
  return tool.main(argc, argv);
}

/// @endcond
