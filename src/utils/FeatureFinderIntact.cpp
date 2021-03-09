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
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/FILTERING/DATAREDUCTION/FeatureFindingMetabo.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <OpenMS/ANALYSIS/QUANTITATION/FeatureFindingIntact.h>


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

private:
  double local_rt_range_;
  double local_mz_range_;
  double charge_lower_bound_;
  double charge_upper_bound_;
  bool use_smoothed_intensities_;

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input file");
    setValidFormats_("in", ListUtils::create<String>("mzML"));

//    registerOutputFile_("out", "<file>", "", "FeatureXML file with metabolite features");
//    setValidFormats_("out", ListUtils::create<String>("featureXML"));
  }

public:
  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------
    // TODO: need to update this value (came from FeatureFindingMetabo)
    local_mz_range_ = 6.5 ; // MZ range where to look for isotopic mass traces (-> decides size of isotopes =(local_mz_range_ * charge))
    local_rt_range_ = 15.0 ; // RT range where to look for coeluting mass traces
    charge_lower_bound_ = 7;
    charge_upper_bound_ = 30;
    use_smoothed_intensities_ = true; // for intensity of a mass trace

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------
    String in = getStringOption_("in");

    MzMLFile mz_data_file;
    mz_data_file.setLogType(log_type_);
    PeakMap ms_peakmap;
    std::vector<Int> ms_level(1, 1);
    mz_data_file.getOptions().setMSLevels(ms_level);
    /// for test purpose : reduce in_ex
//    mz_data_file.getOptions().setMZRange(DRange<1>(DPosition<1>(1230.0), DPosition<1>(1455)));
//    mz_data_file.getOptions().setRTRange(DRange<1>(DPosition<1>(130.0), DPosition<1>(410)));
    mz_data_file.load(in, ms_peakmap);

    if (ms_peakmap.empty())
    {
      OPENMS_LOG_WARN << "The given file does not contain any conventional peak data, but might"
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

    //-------------------------------------------------------------
    // Mass traces detection
    //-------------------------------------------------------------
    vector<MassTrace> m_traces;
    MassTraceDetection mtdet;
    mtdet.run(ms_peakmap, m_traces);

    //-------------------------------------------------------------
    // Elution peak detection
    //-------------------------------------------------------------
    std::vector<MassTrace> m_traces_final;
    std::vector<MassTrace> splitted_mtraces;
    ElutionPeakDetection epdet;
    // fill mass traces with smoothed data as well .. bad design..
    epdet.detectPeaks(m_traces, splitted_mtraces);
//    epdet.filterByPeakWidth(splitted_mtraces, m_traces_final);
    m_traces_final = splitted_mtraces;
    OPENMS_LOG_INFO << "# input mass traces : " << m_traces_final.size() << endl;

    //-------------------------------------------------------------
    // Feature finding
    //-------------------------------------------------------------
//    std::vector<FeatureHypothesis> feat_hypos;
//    build_feature_hypotheses_(m_traces_final, feat_hypos);
    FeatureFindingIntact ffi;
    ffi.updateMembers_(); // TODO: change this with param handler
    FeatureMap out_map;
    ffi.run(m_traces_final, out_map);

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPFeatureFinderIntact tool;
  return tool.main(argc, argv);
}

/// @endcond