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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------
#include <OpenMS/config.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>
#include <OpenMS/FORMAT/EDTAFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/FILTERING/SMOOTHING/GaussFilter.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerCWT.h>

#include <functional>
#include <numeric>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page TOPP_EICExtractor EICExtractor

    @brief Extracts EICs from an MS experiment, in order to quantify analytes at a given position

<CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ EICExtractor \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FileConverter</td>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> statistical tools, e.g., Excel, R, ... </td>
        </tr>
    </table>
</CENTER>

  Use this instead of FeatureFinder, if you have bad features  which are not recognized (much noise etc)
    or if you want to quantify non-peptides.

    The input EDTA file specifies where to search for signal in RT and m/z.
    Retention time is in seconds [s]. A third intensity column is ignored but needs to be present.

Example (replace space separator with &lt;TAB&gt;):<br>
@code
RT m/z int
19.2 431.85 0
21.1 678.77 0
25.7 660.76 0
59.2 431.85 0
@endcode

  RT positions can also be automatically generated using the 'auto-RT' functionality, which can be enabled by the flag 'auto_rt:enabled'.
  All EDTA input lines with negative RT and some m/z values are replaced by 'n' other lines, where the m/z value is identical
  and the RT column is replaced by of the 'n' RT estimates as discovered by the auto-RT functionality.
  This allows you to specify only the expected m/z positions and let the auto-RT function handle the RT positions.
  Info: auto-RT positions are only generated once from the FIRST mzML input file. All other mzML input files are expected to
        have similar RT positions!
  To debug auto-RT, you can specify an mzML output file (see 'auto_rt:out_debug_TIC' option) which will contain four single spectra which represent:
    1. the TIC (of the first mzML input file)
    2. the smoothed version of #1
    3. the signal/noise (S/N) ratio of #2
    4. the centroided version of #2 (not #3!)
  Since you can specify the smoothing aggressiveness using 'auto_rt:FHWM' and
  the minimum S/N theshold for centroided using 'auto_rt:SNThreshold', this should give you all information needed to set the best parameters which fit your data.
  Sensible default thresholds have been chosen though, such that adaption should only be required in extreme cases.

    The intensity reported is the MAXIMUM intensity of all peaks each within the given tolerances for this row's position.

  As output, one file in text format is given. It contains the actual RT and m/z positions of the data,
  as well as RT delta (in [s]) and m/z delta (in ppm) from the expected position as specified in the EDTA file or as found by the auto-RT feature.

<pre>
  RT	  - expected RT position (in [s])
  mz    - expected m/z position
  RTobs - RT position (in [s]) of the quantified entity
  dRT	  - RT delta (in [s]) to the input RT value (as specified in input file or as computed by the auto-rt heuristic)
  mzobs - m/z position of the quantified entity
  dppm  - m/z delta (in parts-per-million) to the input m/z value (as specified in input file) intensity quantification (height of centroided peak); this is an average over multiple scans, thus does usually not correspond to the maximum peak ppm
  intensity
 </pre>

  Each input experiment gives rise to the two RT and mz columns plus additional five columns (starting from RTobs) for each input file.


  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_EICExtractor.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_EICExtractor.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


struct HeaderInfo
{
  explicit HeaderInfo(const String& filename)
  {
    header_description = "-- empty --";
    TextFile tf;
    tf.load(filename);
    String content;
    content.concatenate(tf.begin(), tf.end(), ";");

    String search = "$$ Sample Description:";
    Size pos = content.find(search);
    if (pos != std::string::npos)
    {
      pos += search.size();
      Size pos_end = content.find("$$", pos);
      if (pos_end != std::string::npos)
      {
        String tmp = content.substr(pos, pos_end - pos - 1);
        if (!tmp.trim().empty()) header_description = tmp;
        //std::cerr << "Header info is: " << header_description << std::endl;
      }
    }
  }

  String header_description;
  String filename;
};


class TOPPEICExtractor :
  public TOPPBase
{
public:
  TOPPEICExtractor() :
    TOPPBase("EICExtractor", "Extracts intensities from dedicates positions in a LC/MS map")
  {
  }

  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<file>", ListUtils::create<String>(""), "Input raw data file");
    setValidFormats_("in", ListUtils::create<String>("mzML"));

    registerInputFileList_("in_header", "<file>", ListUtils::create<String>(""), "[for Waters data only] Read additional information from _HEADER.TXT. Provide one for each raw input file.", false);
    setValidFormats_("in_header", ListUtils::create<String>("txt"));

    registerInputFile_("pos", "<file>", "", "Input config file stating where to find signal");
    setValidFormats_("pos", ListUtils::create<String>("edta"));
    registerDoubleOption_("rt_tol", "", 3, "RT tolerance in [s] for finding max peak (whole RT range around RT middle)", false, false);
    registerDoubleOption_("mz_tol", "", 10, "m/z tolerance in [ppm] for finding a peak", false, false);
    registerIntOption_("rt_collect", "", 1, "# of scans up & down in RT from highest point for ppm estimation in result", false, false);

    registerTOPPSubsection_("auto_rt", "Parameters for automatic detection of injection RT peaks (no need to specify them in 'pos' input file)");
    registerFlag_("auto_rt:enabled", "Automatically detect injection peaks from TIC and quantify all m/z x RT combinations.");
    registerDoubleOption_("auto_rt:FHWM", "<FWHM [s]>", 5, "Expected full width at half-maximum of each raw RT peak in [s]. Gaussian smoothing filter with this width is applied to TIC.", false, true);
    registerDoubleOption_("auto_rt:SNThreshold", "<S/N>", 5, "S/N threshold for a smoothed raw peak to pass peak picking. Higher thesholds will result in less peaks.", false, true);
    registerOutputFile_("auto_rt:out_debug_TIC", "<file>", "", "Optional output file (for first input) containing the smoothed TIC, S/N levels and picked RT positions", false, true);
    setValidFormats_("auto_rt:out_debug_TIC", ListUtils::create<String>("mzML"));

    registerStringOption_("out_separator", "<sep>", ",", "Separator character for output CSV file.", false, true);
    //setValidStrings_("out_separator", ListUtils::create<String>(",!\t! ", '!')); // comma not allowed as valid string

    registerOutputFile_("out", "<file>", "", "Output quantitation file (multiple columns for each input compound)");
    setValidFormats_("out", ListUtils::create<String>("csv"));
  }

  MSChromatogram toChromatogram(const MSSpectrum& in)
  {
    MSChromatogram out;
    for (Size ic = 0; ic < in.size(); ++ic)
    {
      ChromatogramPeak peak;
      peak.setMZ(in[ic].getMZ());
      peak.setIntensity(in[ic].getIntensity());
      out.push_back(peak);
    }
    out.setChromatogramType(ChromatogramSettings::SELECTED_ION_CURRENT_CHROMATOGRAM);

    return out;
  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------
    StringList in = getStringList_("in");
    String edta = getStringOption_("pos");
    String out = getStringOption_("out");
    String out_sep = getStringOption_("out_separator");
    String out_TIC_debug = getStringOption_("auto_rt:out_debug_TIC");

    StringList in_header = getStringList_("in_header");


    // number of out_debug_TIC files and input files must be identical
    /*if (out_TIC_debug.size() > 0 && in.size() != out_TIC_debug.size())
    {
        LOG_FATAL_ERROR << "Error: number of input file 'in' and auto_rt:out_debug_TIC files must be identical!" << std::endl;
        return ILLEGAL_PARAMETERS;
    }*/

    // number of header files and input files must be identical
    if (in_header.size() > 0 && in.size() != in_header.size())
    {
      LOG_FATAL_ERROR << "Error: number of input file 'in' and 'in_header' files must be identical!" << std::endl;
      return ILLEGAL_PARAMETERS;
    }

    if (!getFlag_("auto_rt:enabled") && !out_TIC_debug.empty())
    {
      LOG_FATAL_ERROR << "Error: TIC output file requested, but auto_rt is not enabled! Either do not request the file or switch on 'auto_rt:enabled'." << std::endl;
      return ILLEGAL_PARAMETERS;
    }

    double rttol = getDoubleOption_("rt_tol");
    double mztol = getDoubleOption_("mz_tol");
    Size rt_collect = getIntOption_("rt_collect");

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------
    MzMLFile mzml_file;
    mzml_file.setLogType(log_type_);
    PeakMap exp, exp_pp;

    EDTAFile ed;
    ConsensusMap cm;
    ed.load(edta, cm);

    StringList tf_single_header0, tf_single_header1, tf_single_header2; // header content, for each column

    std::vector<String> vec_single; // one line for each compound, multiple columns per experiment
    vec_single.resize(cm.size());
    for (Size fi = 0; fi < in.size(); ++fi)
    {
      // load raw data
      mzml_file.load(in[fi], exp);
      exp.sortSpectra(true);

      if (exp.empty())
      {
        LOG_WARN << "The given file does not contain any conventional peak data, but might"
                    " contain chromatograms. This tool currently cannot handle them, sorry." << std::endl;
        return INCOMPATIBLE_INPUT_DATA;
      }

      // try to detect RT peaks (only for the first input file -- all others should align!)
      // cm.size() might change in here...
      if (getFlag_("auto_rt:enabled") && fi == 0)
      {
        ConsensusMap cm_local = cm; // we might have different RT peaks for each map if 'auto_rt' is enabled
        cm.clear(false); // reset global list (about to be filled)

        // compute TIC
        MSChromatogram tic = exp.getTIC();
        MSSpectrum tics, tic_gf, tics_pp, tics_sn;
        for (Size ic = 0; ic < tic.size(); ++ic)
        { // rewrite Chromatogram to MSSpectrum (GaussFilter requires it)
          Peak1D peak;
          peak.setMZ(tic[ic].getRT());
          peak.setIntensity(tic[ic].getIntensity());
          tics.push_back(peak);
        }
        // smooth (no PP_CWT here due to efficiency reasons -- large FWHM take longer!)
        double fwhm = getDoubleOption_("auto_rt:FHWM");
        GaussFilter gf;
        Param p = gf.getParameters();
        p.setValue("gaussian_width", fwhm * 2); // wider than FWHM, just to be sure we have a fully smoothed peak. Merging two peaks is unlikely
        p.setValue("use_ppm_tolerance", "false");
        gf.setParameters(p);
        tic_gf = tics;
        gf.filter(tic_gf);
        // pick peaks
        PeakPickerHiRes pp;
        p = pp.getParameters();
        p.setValue("signal_to_noise", getDoubleOption_("auto_rt:SNThreshold"));
        pp.setParameters(p);
        pp.pick(tic_gf, tics_pp);

        if (tics_pp.size())
        {
          LOG_INFO << "Found " << tics_pp.size() << " auto-rt peaks at: ";
          for (Size ipp = 0; ipp != tics_pp.size(); ++ipp) LOG_INFO << " " << tics_pp[ipp].getMZ();
        }
        else
        {
          LOG_INFO << "Found no auto-rt peaks. Change threshold parameters!";
        }
        LOG_INFO << std::endl;

        if (!out_TIC_debug.empty()) // if debug file was given
        { // store intermediate steps for debug
          PeakMap out_debug;
          out_debug.addChromatogram(toChromatogram(tics));
          out_debug.addChromatogram(toChromatogram(tic_gf));

          SignalToNoiseEstimatorMedian<MSSpectrum> snt;
          snt.init(tics);
          for (Size is = 0; is < tics.size(); ++is)
          {
            Peak1D peak;
            peak.setMZ(tic[is].getMZ());
            peak.setIntensity(snt.getSignalToNoise(tics[is]));
            tics_sn.push_back(peak);
          }
          out_debug.addChromatogram(toChromatogram(tics_sn));

          out_debug.addChromatogram(toChromatogram(tics_pp));
          // get rid of "native-id" missing warning
          for (Size id = 0; id < out_debug.size(); ++id) out_debug[id].setNativeID(String("spectrum=") + id);

          mzml_file.store(out_TIC_debug, out_debug);
          LOG_DEBUG << "Storing debug AUTO-RT: " << out_TIC_debug << std::endl;
        }

        // add target EICs: for each m/z with no/negative RT, add all combinations of that m/z with auto-RTs
        // duplicate m/z entries will be ignored!
        // all other lines with positive RT values are copied unaffected
        //do not allow doubles
        std::set<double> mz_doubles;
        for (ConsensusMap::Iterator cit = cm_local.begin(); cit != cm_local.end(); ++cit)
        {
          if (cit->getRT() < 0)
          {
            if (mz_doubles.find(cit->getMZ()) == mz_doubles.end())
            {
              mz_doubles.insert(cit->getMZ());
            }
            else
            {
              LOG_INFO << "Found duplicate m/z entry (" << cit->getMZ() << ") for auto-rt. Skipping ..." << std::endl;
              continue;
            }

            ConsensusMap cm_RT_multiplex;
            for (MSSpectrum::ConstIterator itp = tics_pp.begin(); itp != tics_pp.end(); ++itp)
            {
              ConsensusFeature f = *cit;
              f.setRT(itp->getMZ());
              cm.push_back(f);
            }

          }
          else
          { // default feature with no auto-rt
            LOG_INFO << "copying feature with RT " << cit->getRT() << std::endl;
            cm.push_back(*cit);
          }
        }

        // resize, since we have more positions now
        vec_single.resize(cm.size());
      }


      // search for each EIC and add up
      Int not_found(0);
      Map<Size, double> quant;

      String description;
      if (fi < in_header.size())
      {
        HeaderInfo info(in_header[fi]);
        description = info.header_description;
      }

      if (fi == 0)
      { // two additional columns for first file (theoretical RT and m/z)
        tf_single_header0 << "" << "";
        tf_single_header1 << "" << "";
        tf_single_header2 << "RT" << "mz";
      }

      // 5 entries for each input file
      tf_single_header0 << File::basename(in[fi]) << "" << "" << "" << "";
      tf_single_header1 << description << "" << "" << "" << "";
      tf_single_header2 << "RTobs" << "dRT" << "mzobs" << "dppm" << "intensity";

      for (Size i = 0; i < cm.size(); ++i)
      {
        //std::cerr << "Rt" << cm[i].getRT() << "  mz: " << cm[i].getMZ() << " R " <<  cm[i].getMetaValue("rank") << "\n";

        double mz_da = mztol * cm[i].getMZ() / 1e6; // mz tolerance in Dalton
        PeakMap::ConstAreaIterator it = exp.areaBeginConst(cm[i].getRT() - rttol / 2,
                                                                  cm[i].getRT() + rttol / 2,
                                                                  cm[i].getMZ() - mz_da,
                                                                  cm[i].getMZ() + mz_da);
        Peak2D max_peak;
        max_peak.setIntensity(0);
        max_peak.setRT(cm[i].getRT());
        max_peak.setMZ(cm[i].getMZ());
        for (; it != exp.areaEndConst(); ++it)
        {
          if (max_peak.getIntensity() < it->getIntensity())
          {
            max_peak.setIntensity(it->getIntensity());
            max_peak.setRT(it.getRT());
            max_peak.setMZ(it->getMZ());
          }
        }
        double ppm = 0; // observed m/z offset

        if (max_peak.getIntensity() == 0)
        {
          ++not_found;
        }
        else
        {
          // take median for m/z found
          std::vector<double> mz;
          PeakMap::Iterator itm = exp.RTBegin(max_peak.getRT());
          SignedSize low = std::min<SignedSize>(std::distance(exp.begin(), itm), rt_collect);
          SignedSize high = std::min<SignedSize>(std::distance(itm, exp.end()) - 1, rt_collect);
          PeakMap::AreaIterator itt = exp.areaBegin((itm - low)->getRT() - 0.01, (itm + high)->getRT() + 0.01, cm[i].getMZ() - mz_da, cm[i].getMZ() + mz_da);
          for (; itt != exp.areaEnd(); ++itt)
          {
            mz.push_back(itt->getMZ());
            //std::cerr << "ppm: " << itt.getRT() << " " <<  itt->getMZ() << " " << itt->getIntensity() << std::endl;
          }

          if ((SignedSize)mz.size() > (low + high + 1)) LOG_WARN << "Compound " << i << " has overlapping peaks [" << mz.size() << "/" << low + high + 1 << "]" << std::endl;

          if (!mz.empty())
          {
            double avg_mz = std::accumulate(mz.begin(), mz.end(), 0.0) / double(mz.size());
            //std::cerr << "avg: " << avg_mz << "\n";
            ppm = (avg_mz - cm[i].getMZ()) / cm[i].getMZ() * 1e6;
          }

        }

        // appending the second column set requires separator
        String append_sep = (fi == 0 ? "" : out_sep);

        vec_single[i] += append_sep; // new line
        if (fi == 0)
        {
          vec_single[i] += String(cm[i].getRT()) + out_sep +
                           String(cm[i].getMZ()) + out_sep;
        }
        vec_single[i] += String(max_peak.getRT()) + out_sep +
                         String(max_peak.getRT() - cm[i].getRT()) + out_sep +
                         String(max_peak.getMZ()) + out_sep +
                         String(ppm)  + out_sep +
                         String(max_peak.getIntensity());
      }

      if (not_found) LOG_INFO << "Missing peaks for " << not_found << " compounds in file '" << in[fi] << "'.\n";
    }

    //-------------------------------------------------------------
    // create header
    //-------------------------------------------------------------
    vec_single.insert(vec_single.begin(), ListUtils::concatenate(tf_single_header2, out_sep));
    vec_single.insert(vec_single.begin(), ListUtils::concatenate(tf_single_header1, out_sep));
    vec_single.insert(vec_single.begin(), ListUtils::concatenate(tf_single_header0, out_sep));

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    TextFile tf;
    for (std::vector<String>::iterator v_it = vec_single.begin(); v_it != vec_single.end(); ++v_it)
    {
      tf.addLine(*v_it);
    }
    tf.store(out);

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPEICExtractor tool;
  return tool.main(argc, argv);
}

/// @endcond
