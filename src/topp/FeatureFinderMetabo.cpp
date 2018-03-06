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

  Mass traces alone would allow for further analysis such as metabolite ID or
  statistical evaluation. However, in general, monoisotopic mass traces are
  accompanied by satellite C13 peaks and thus may render the analysis more
  difficult. @ref FeatureFinderMetabo fulfills a further data reduction step by
  assembling compatible mass traces to metabolite features (that is, all mass
  traces originating from one metabolite). To this end, multiple metabolite
  hypotheses are formulated and scored according to how well differences in RT (optional),
  m/z or intensity ratios match to those of theoretical isotope patterns.

  If the raw data scans contain the scan polarity information, it is stored as
  meta value "scan_polarity" in the output file.

  Mass trace clustering can be done using either 13C distances or a linear model (Kenar et al) -- see parameter 'ffm:mz_scoring_13C'.
  Generally, for lipidomics, use 13C, since lipids contain a lot of 13C.
  For general metabolites, the linear model is usually more appropriate.
  To decide what is better, the total number of features can be used as indirect measure
  - the lower(!) the better (since more mass traces are assembled into single features).
  Detailed information is stored in the featureXML output: it contains meta-values for each feature about the 
  mass trace differences (inspectable via TOPPView). If you want this in a tabular format, use TextExporter, i.e.,
  @code
     TextExporter.exe -feature:add_metavalues 1 -in <ff_metabo.featureXML> -out <ff_metabo.csv>
  @endcode
  By default, the linear model is used.

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_FeatureFinderMetabo.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_FeatureFinderMetabo.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFeatureFinderMetabo :
  public TOPPBase
{
public:
  TOPPFeatureFinderMetabo() :
    TOPPBase("FeatureFinderMetabo", "Assembles metabolite features from centroided (LC-)MS data using the mass trace approach.")
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Centroided mzML file");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "FeatureXML file with metabolite features");
    setValidFormats_("out", ListUtils::create<String>("featureXML"));

    registerOutputFile_("out_chrom", "<file>", "", "Optional mzML file with chromatograms", false);
    setValidFormats_("out_chrom", ListUtils::create<String>("mzML"));

    addEmptyLine_();
    registerSubsection_("algorithm", "Algorithm parameters section");
  }

  Param getSubsectionDefaults_(const String& /*section*/) const override
  {
    Param combined;
    Param p_com;
    p_com.setValue("noise_threshold_int", 10.0, "Intensity threshold below which peaks are regarded as noise.");
    p_com.setValue("chrom_peak_snr", 3.0, "Minimum signal-to-noise a mass trace should have.");
    p_com.setValue("chrom_fwhm", 5.0, "Expected chromatographic peak width (in seconds).");
    combined.insert("common:", p_com);
    combined.setSectionDescription("common", "Common parameters for all other subsections");

    Param p_mtd = MassTraceDetection().getDefaults();
    p_mtd.remove("noise_threshold_int");
    p_mtd.remove("chrom_peak_snr");
    combined.insert("mtd:", p_mtd);
    combined.setSectionDescription("mtd", "Mass Trace Detection parameters");

    Param p_epd;
    p_epd.setValue("enabled", "true", "Enable splitting of isobaric mass traces by chromatographic peak detection. Disable for direct injection.");
    p_epd.setValidStrings("enabled", ListUtils::create<String>("true,false"));
    p_epd.insert("", ElutionPeakDetection().getDefaults());
    p_epd.remove("chrom_peak_snr");
    p_epd.remove("chrom_fwhm");

    combined.insert("epd:", p_epd);
    combined.insert("epd:", p_epd);
    combined.setSectionDescription("epd", "Elution Profile Detection (to separate isobaric Mass Traces by elution time).");

    Param p_ffm = FeatureFindingMetabo().getDefaults();
    p_ffm.remove("chrom_fwhm");
    p_ffm.remove("report_chromatograms");
    combined.insert("ffm:", p_ffm);
    combined.setSectionDescription("ffm", "FeatureFinder parameters (assembling mass traces to charged features)");

    return combined;
  }

  ExitCodes main_(int, const char**) override
  {

    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    String in = getStringOption_("in");
    String out = getStringOption_("out");
    String out_chrom = getStringOption_("out_chrom");

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------
    MzMLFile mz_data_file;
    mz_data_file.setLogType(log_type_);
    PeakMap ms_peakmap;
    std::vector<Int> ms_level(1, 1);
    mz_data_file.getOptions().setMSLevels(ms_level);
    mz_data_file.load(in, ms_peakmap);

    if (ms_peakmap.empty())
    {
      LOG_WARN << "The given file does not contain any conventional peak data, but might"
                  " contain chromatograms. This tool currently cannot handle them, sorry.";
      return INCOMPATIBLE_INPUT_DATA;
    }

    // determine type of spectral data (profile or centroided)
    SpectrumSettings::SpectrumType spectrum_type = ms_peakmap[0].getType();

    if (spectrum_type == SpectrumSettings::PROFILE)
    {
      if (!getFlag_("force"))
      {
        throw OpenMS::Exception::FileEmpty(__FILE__, __LINE__, __FUNCTION__,
            "Error: Profile data provided but centroided spectra expected. To enforce processing of the data set the -force flag.");
      }
    }

    // make sure the spectra are sorted by m/z
    ms_peakmap.sortSpectra(true);

    vector<MassTrace> m_traces;

    //-------------------------------------------------------------
    // set parameters
    //-------------------------------------------------------------

    Param common_param = getParam_().copy("algorithm:common:", true);
    writeDebug_("Common parameters passed to sub-algorithms (mtd and ffm)", common_param, 3);

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

    std::vector<MassTrace> m_traces_final;
    if (epd_param.getValue("enabled").toBool())
    {
      std::vector<MassTrace> splitted_mtraces;
      epd_param.remove("enabled"); // artificially added above
      epd_param.insert("", common_param);
      ElutionPeakDetection epdet;
      epdet.setParameters(epd_param);
      // fill mass traces with smoothed data as well .. bad design..
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
    }
    else // no elution peak detection
    {
      m_traces_final = m_traces;
      for (Size i = 0; i < m_traces_final.size(); ++i) // estimate FWHM, so .getIntensity() can be called later
      {
        m_traces_final[i].estimateFWHM(false);
      }
      if (ffm_param.getValue("use_smoothed_intensities").toBool())
      {
        LOG_WARN << "Without EPD, smoothing is not supported. Setting 'use_smoothed_intensities' to false!" << std::endl;
        ffm_param.setValue("use_smoothed_intensities", "false");
      }
    }

    //-------------------------------------------------------------
    // configure and run feature finding
    //-------------------------------------------------------------

    ffm_param.insert("", common_param);
    ffm_param.remove("noise_threshold_int");
    ffm_param.remove("chrom_peak_snr");
    String report_chromatograms = out_chrom.empty() ? "false" : "true";
    ffm_param.setValue("report_chromatograms", report_chromatograms);

    FeatureMap feat_map;
    StringList ms_runs;
    ms_peakmap.getPrimaryMSRunPath(ms_runs);
    feat_map.setPrimaryMSRunPath(ms_runs);

    std::vector< std::vector< OpenMS::MSChromatogram > > feat_chromatograms;

    FeatureFindingMetabo ffmet;
    ffmet.setParameters(ffm_param);
    ffmet.run(m_traces_final, feat_map, feat_chromatograms);

    Size trace_count(0);
    for (Size i = 0; i < feat_map.size(); ++i)
    {
      OPENMS_PRECONDITION(feat_map[i].metaValueExists("num_of_masstraces"),
          "MetaValue 'num_of_masstraces' missing from FFMetabo output!");
      trace_count += (Size) feat_map[i].getMetaValue("num_of_masstraces");
    }

    LOG_INFO << "-- FF-Metabo stats --\n"
             << "Input traces:    " << m_traces_final.size() << "\n"
             << "Output features: " << feat_map.size() << " (total trace count: " << trace_count << ")" << std::endl;

    if (trace_count != m_traces_final.size())
    {
        LOG_ERROR << "FF-Metabo: Internal error. Not all mass traces have been assembled to features! Aborting." << std::endl;
        return UNEXPECTED_RESULT;
    }

    // store chromatograms
    if (!out_chrom.empty())
    {
        if (feat_chromatograms.size() == feat_map.size())
        {
          MSExperiment out_exp;
            for (Size i = 0; i < feat_chromatograms.size(); ++i)
            {
                for (Size j = 0; j < feat_chromatograms[i].size(); ++j)
                {
                    out_exp.addChromatogram(feat_chromatograms[i][j]);
                }
            }
          MzMLFile().store(out_chrom, out_exp);
        }
        else
        {
            LOG_ERROR << "FF-Metabo: Internal error. The number of features (" << feat_chromatograms.size() << ") and chromatograms (" << feat_map.size() << ") are different! Aborting." << std::endl;
            return UNEXPECTED_RESULT;
        }
    }

    // store ionization mode of spectra (useful for post-processing by AccurateMassSearch tool)
    if (feat_map.size() > 0)
    {
      set<IonSource::Polarity> pols;
      for (Size i = 0; i < ms_peakmap.size(); ++i)
      {
        pols.insert(ms_peakmap[i].getInstrumentSettings().getPolarity());
      }
      // concat to single string
      StringList sl_pols;
      for (set<IonSource::Polarity>::const_iterator it = pols.begin(); it != pols.end(); ++it)
      {
        sl_pols.push_back(String(IonSource::NamesOfPolarity[*it]));
      }
      feat_map[0].setMetaValue("scan_polarity", ListUtils::concatenate(sl_pols, ";"));
    }

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    // annotate output with data processing info
    addDataProcessing_(feat_map, getProcessingInfo_(DataProcessing::QUANTITATION));

    FeatureXMLFile feature_xml_file;
    feature_xml_file.setLogType(log_type_);
    feature_xml_file.store(out, feat_map);

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPFeatureFinderMetabo tool;
  return tool.main(argc, argv);
}

/// @endcond
