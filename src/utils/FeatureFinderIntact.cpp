// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Jihyung Kim $
// $Authors: Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/FILTERING/DATAREDUCTION/MassTraceDetection.h>
#include <OpenMS/FILTERING/DATAREDUCTION/ElutionPeakDetection.h>
#include <OpenMS/FILTERING/DATAREDUCTION/FeatureFindingMetabo.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <OpenMS/ANALYSIS/QUANTITATION/FeatureFindingIntact.h>
#include <OpenMS/ANALYSIS/QUANTITATION/FLASHDeconvQuant.h>

// for testing
#include <OpenMS/FORMAT/FeatureXMLFile.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------
/**
    @page TOPP_FeatureFinderIntact TOPP_FeatureFinderIntact

    @brief TOPP_FeatureFinderIntact The intact protein feature detection for quantification (centroided).
 */

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFeatureFinderIntact :
    public TOPPBase,
    public ProgressLogger
{
public:
  TOPPFeatureFinderIntact():
    TOPPBase("FeatureFinderIntact", "The intact protein feature detection for quantification", false, {}, false), ProgressLogger()
  {
    this->setLogType(CMD);
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input file", true);
    setValidFormats_("in", ListUtils::create<String>("mzML"));

    registerOutputFile_("out", "<string>", "", "FeatureXML file with metabolite features", false);
//    setValidFormats_("out", ListUtils::create<String>("featureXML"));

    addEmptyLine_();
    registerSubsection_("algorithm", "Algorithm parameters section");
  }

  Param getSubsectionDefaults_(const String& /*section*/) const override
  {
    Param combined;

    // TODO: add mtd, epd param

    Param p_ffi = FeatureFindingIntact().getDefaults();
    combined.insert("ffi:", p_ffi);
    combined.setSectionDescription("ffi", "FeatureFinder parameters (assembling mass traces to charged features)");

    return combined;
  }

public:
  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------
    // TODO: need to update this value (came from FeatureFindingMetabo)
    Param ffi_param = getParam_().copy("algorithm:ffi:", true);

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------
    String in = getStringOption_("in");
    String out_dir = getStringOption_("out");

    MzMLFile mz_data_file;
    mz_data_file.setLogType(log_type_);
    PeakMap ms_peakmap;
    std::vector<Int> ms_level(1, 1);
    mz_data_file.getOptions().setMSLevels(ms_level);
    mz_data_file.load(in, ms_peakmap);

    if (ms_peakmap.empty())
    {
      OPENMS_LOG_WARN << "The given file does not contain any conventional peak data, but might"
                         " contain chromatograms. This tool currently cannot handle them, sorry.";
      return INCOMPATIBLE_INPUT_DATA;
    }
    OPENMS_LOG_INFO << "using " << ms_peakmap.getNrSpectra() << " MS1 spectra" << endl;

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

    //-------------------------------------------------------------
    // Mass traces detection
    //-------------------------------------------------------------
    Param p_mtd = MassTraceDetection().getDefaults();
    p_mtd.setValue("noise_threshold_int" , 0.0);
    p_mtd.setValue("chrom_peak_snr" , 0.0);
    p_mtd.setValue("mass_error_ppm", 5.0);
//    p_mtd.setValue("reestimate_mt_sd", "false");
//    p_mtd.setValue("trace_termination_criterion", "sample_rate");
//    p_mtd.setValue("min_sample_rate", 0.2);

    vector<MassTrace> m_traces;
    MassTraceDetection mtdet;
    mtdet.setParameters(p_mtd);
    mtdet.run(ms_peakmap, m_traces);
    OPENMS_LOG_INFO << "# initial input mass traces : " << m_traces.size() << endl;

    //-------------------------------------------------------------
    // Elution peak detection
    //-------------------------------------------------------------
    Param p_epd = ElutionPeakDetection().getDefaults();
    p_epd.setValue("width_filtering", "off");

    std::vector<MassTrace> m_traces_final;
    ElutionPeakDetection epdet;
    epdet.setParameters(p_epd);
    // fill mass traces with smoothed data as well .. bad design..
    epdet.detectPeaks(m_traces, m_traces_final);
    for (auto &m : m_traces)
    {
      m.estimateFWHM(false);
    }
    OPENMS_LOG_INFO << "# final input mass traces : " << m_traces_final.size() << endl;

    // for test output TODO: remove
    Size last_index = in.find_last_of(".");
//    String out = in.substr(0, last_index) + ".mt.FeatureXML";
//    write_mtraces(in, m_traces_final, out);

    //-------------------------------------------------------------
    // Feature finding
    //-------------------------------------------------------------
    FLASHDeconvQuant fdq;
    fdq.setParameters(ffi_param);
    if (!out_dir.empty())
    {
      Size dirindex = in.find_last_of("/");
      String dirpath = in.substr(0, dirindex) + "/" + out_dir + "/";

      String filename = in.substr(dirindex+1);
      Size file_extension_idx = filename.find_last_of(".");
      filename = filename.substr(0, file_extension_idx) + ".features.tsv";
      fdq.outfile_path = dirpath + filename;
    }
    else
    {
      fdq.outfile_path = in.substr(0, last_index) + ".features.tsv";
    }
    FeatureMap out_map;
    fdq.run(m_traces_final, out_map);

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    return EXECUTION_OK;
  }

  // copied from MassTraceExtractor
  void write_mtraces(String in, std::vector<MassTrace> m_traces_final, String out)
  {
    std::vector<double> stats_sd;
    FeatureMap ms_feat_map;
    ms_feat_map.setPrimaryMSRunPath({in});

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
      f.setWidth(m_traces_final[i].estimateFWHM(true));
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
      OPENMS_LOG_INFO << "Mass trace m/z s.d.\n"
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

};

int main(int argc, const char** argv)
{
  TOPPFeatureFinderIntact tool;
  return tool.main(argc, argv);
}

/// @endcond