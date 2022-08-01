// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/FILTERING/DATAREDUCTION/MassTraceDetection.h>
#include <OpenMS/FILTERING/DATAREDUCTION/ElutionPeakDetection.h>
#include <OpenMS/FILTERING/DATAREDUCTION/FeatureFindingMetabo.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvQuantAlgorithm.h>


using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------
/**
    @page TOPP_FLASHDeconvQ TOPP_FLASHDeconvQ

    @brief TOPP_FLASHDeconvQ The intact protein feature detection for quantification (centroided).
 */

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFLASHDeconvQ :
    public TOPPBase,
    public ProgressLogger
{
public:
  TOPPFLASHDeconvQ():
    TOPPBase("FLASHDeconvQ", "The intact protein feature detection for quantification", false, {}, false), ProgressLogger()
  {
    this->setLogType(CMD);
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input file (mzML)", true);
    setValidFormats_("in", ListUtils::create<String>("mzML"));

    registerOutputFile_("out", "<file>", "", "feature level quantification output tsv file", true);
    setValidFormats_("out", ListUtils::create<String>("tsv"));

    registerOutputFile_("out_feat", "<file>", "", "featureXML format feature level quantification output file", false);
    setValidFormats_("out_feat", ListUtils::create<String>("featureXML"));

    addEmptyLine_();
    registerSubsection_("algorithm", "Algorithm parameters section");
  }

  Param getSubsectionDefaults_(const String& /*section*/) const override
  {
    Param combined;

    Param p_mtd = MassTraceDetection().getDefaults();
    p_mtd.setValue("noise_threshold_int", 0.0);
    p_mtd.setValue("chrom_peak_snr", 0.0);
    p_mtd.setValue("mass_error_ppm", 5.0);
    combined.insert("mtd:", p_mtd);
    combined.setSectionDescription("mtd", "Mass Trace Detection parameters");

    Param p_epd = ElutionPeakDetection().getDefaults();
    p_epd.setValue("width_filtering", "off");
    combined.insert("epd:", p_epd);
    combined.setSectionDescription("epd", "Elution Profile Detection (to separate isobaric Mass Traces by elution time).");

    Param p_ffi = FLASHDeconvQuantAlgorithm().getDefaults();
    combined.insert("fdq:", p_ffi);
    combined.setSectionDescription("fdq", "FLASHDeconvQ parameters (assembling mass traces to charged features)");

    return combined;
  }

public:
  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------
    String in = getStringOption_("in");
    String out = getStringOption_("out");
    String out_feat = getStringOption_("out_feat");

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
    // parameter handling
    //-------------------------------------------------------------
    Param mtd_param = getParam_().copy("algorithm:mtd:", true);
    writeDebug_("Parameters passed to MassTraceDetection", mtd_param, 3);

    Param epd_param = getParam_().copy("algorithm:epd:", true);
    writeDebug_("Parameters passed to ElutionPeakDetection", epd_param, 3);

    Param fdq_param = getParam_().copy("algorithm:fdq:", true);
    writeDebug_("Parameters passed to FLASHDeconvQ", fdq_param, 3);

    //-------------------------------------------------------------
    // Mass traces detection
    //-------------------------------------------------------------
//    Param p_mtd = MassTraceDetection().getDefaults();
//    p_mtd.setValue("noise_threshold_int" , 0.0);
//    p_mtd.setValue("chrom_peak_snr" , 0.0);
//    p_mtd.setValue("mass_error_ppm", 5.0);
//    p_mtd.setValue("trace_termination_criterion", "sample_rate");
//    p_mtd.setValue("min_sample_rate", 0.2);

    vector<MassTrace> m_traces;
    MassTraceDetection mtdet;
    mtdet.setParameters(mtd_param);
    mtdet.run(ms_peakmap, m_traces);
    OPENMS_LOG_INFO << "# initial input mass traces : " << m_traces.size() << endl;

    //-------------------------------------------------------------
    // Elution peak detection
    //-------------------------------------------------------------
//    Param p_epd = ElutionPeakDetection().getDefaults();
//    p_epd.setValue("width_filtering", "off");

    std::vector<MassTrace> m_traces_final;
    ElutionPeakDetection epdet;
    epdet.setParameters(epd_param);
    // fill mass traces with smoothed data as well .. bad design..
    epdet.detectPeaks(m_traces, m_traces_final);

    OPENMS_LOG_INFO << "# final input mass traces : " << m_traces_final.size() << endl;

    //-------------------------------------------------------------
    // Feature finding
    //-------------------------------------------------------------
    FLASHDeconvQuantAlgorithm fdq;
    fdq.setParameters(fdq_param);

    fdq.inputfile_path = in;
    fdq.outfile_path = out;
    if (!out_feat.empty())
    {
      fdq.outFeatureXML = true;
    }
    FeatureMap out_map;
    fdq.run(m_traces_final, out_map);

    //-------------------------------------------------------------
    // writing featureXML output
    //-------------------------------------------------------------
    if (!out_feat.empty())
    {
      out_map.setPrimaryMSRunPath({in});
      addDataProcessing_(out_map, getProcessingInfo_(DataProcessing::QUANTITATION));

      FeatureXMLFile().store(out_feat, out_map);
    }

    return EXECUTION_OK;
  }
};

int main(int argc, const char** argv)
{
  TOPPFLASHDeconvQ tool;
  return tool.main(argc, argv);
}

/// @endcond