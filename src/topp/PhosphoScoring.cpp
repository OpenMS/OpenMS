// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: David Wojnar $
// $Authors: David Wojnar, Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/ANALYSIS/ID/AScore.h>

using namespace OpenMS;
using namespace std;

/**
  @page TOPP_PhosphoScoring PhosphoScoring

  @brief Tool to score phosphorylation sites of peptides.

  <CENTER>
    <table>
      <tr>
        <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
        <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ PhosphoScoring \f$ \longrightarrow \f$</td>
        <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
      </tr>
      <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_MascotAdapter (or other ID engines) </td>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_PeptideIndexer </td>
      </tr>
    </table>
  </CENTER>

  @experimental This TOPP tool is not well tested and not all features might be working correctly.

  This tool performs phosphorylation analysis and site localization.
  Input files are a MS/MS an spectrum file as well as the corresponding identification file.
  Firstly, the peptide identifications are mapped onto the spectra.
  Secondly, the tools uses at the moment an implementation of the Ascore according to Beausoleil <em>et al.<em> in order to localize the most probable phosphorylation sites.

  For details, see:\n
  Beausoleil <em>et al.</em>: <a href="http://dx.doi.org/10.1038/nbt1240">A probability-based approach for high-throughput protein phosphorylation analysis and site localization</a> (Nat. Biotechnol., 2006, PMID: 16964243).


  <!-- <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_PhosphoScoring.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_PhosphoScoring.html -->
*/


// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPPhosphoScoring :
  public TOPPBase
{
public:
  TOPPPhosphoScoring() :
    TOPPBase("PhosphoScoring", "Scores potential phosphorylation sites in order to localize the most probable sites.")
  {
  }

protected:

  // spectrum must not contain 0 intensity peaks and must be sorted by m/z
  template <typename SpectrumType>
  static void deisotopeAndSingleChargeMSSpectrum_(SpectrumType& in, Int min_charge, Int max_charge, double fragment_tolerance, bool fragment_unit_ppm, bool keep_only_deisotoped = false, Size min_isopeaks = 3, Size max_isopeaks = 10, bool make_single_charged = true)
  {
    if (in.empty())
    {
      return;
    }

    SpectrumType old_spectrum = in;

    // determine charge seeds and extend them
    vector<Size> mono_isotopic_peak(old_spectrum.size(), 0);
    vector<Int> features(old_spectrum.size(), -1);
    Int feature_number = 0;

    for (Size current_peak = 0; current_peak != old_spectrum.size(); ++current_peak)
    {
      double current_mz = old_spectrum[current_peak].getPosition()[0];

      for (Int q = max_charge; q >= min_charge; --q)   // important: test charge hypothesis from high to low
      {
        // try to extend isotopes from mono-isotopic peak
        // if extension larger then min_isopeaks possible:
        //   - save charge q in mono_isotopic_peak[]
        //   - annotate all isotopic peaks with feature number
        if (features[current_peak] == -1)   // only process peaks which have no assigned feature number
        {
          bool has_min_isopeaks = true;
          vector<Size> extensions;
          for (Size i = 0; i < max_isopeaks; ++i)
          {
            double expected_mz = current_mz + i * Constants::C13C12_MASSDIFF_U / q;
            Size p = old_spectrum.findNearest(expected_mz);
            double tolerance_dalton = fragment_unit_ppm ? fragment_tolerance * old_spectrum[p].getPosition()[0] * 1e-6 : fragment_tolerance;
            if (fabs(old_spectrum[p].getPosition()[0] - expected_mz) > tolerance_dalton)   // test for missing peak
            {
              if (i < min_isopeaks)
              {
                has_min_isopeaks = false;
              }
              break;
            }
            else
            {
              // TODO: include proper averagine model filtering. for now start at the second peak to test hypothesis
              Size n_extensions = extensions.size();
              if (n_extensions != 0)
              {
                if (old_spectrum[p].getIntensity() > old_spectrum[extensions[n_extensions - 1]].getIntensity())
                {
                  if (i < min_isopeaks)
                  {
                    has_min_isopeaks = false;
                  }
                  break;
                }
              }

              // averagine check passed
              extensions.push_back(p);
            }
          }

          if (has_min_isopeaks)
          {
            //cout << "min peaks at " << current_mz << " " << " extensions: " << extensions.size() << endl;
            mono_isotopic_peak[current_peak] = q;
            for (Size i = 0; i != extensions.size(); ++i)
            {
              features[extensions[i]] = feature_number;
            }
            feature_number++;
          }
        }
      }
    }

    in.clear(false);
    for (Size i = 0; i != old_spectrum.size(); ++i)
    {
      Int z = mono_isotopic_peak[i];
      if (keep_only_deisotoped)
      {
        if (z == 0)
        {
          continue;
        }

        // if already single charged or no decharging selected keep peak as it is
        if (!make_single_charged)
        {
          in.push_back(old_spectrum[i]);
        }
        else
        {
          RichPeak1D p = old_spectrum[i];
          p.setMZ(p.getMZ() * z - (z - 1) * Constants::PROTON_MASS_U);
          in.push_back(p);
        }
      }
      else
      {
        // keep all unassigned peaks
        if (features[i] < 0)
        {
          in.push_back(old_spectrum[i]);
          continue;
        }

        // convert mono-isotopic peak with charge assigned by deisotoping
        if (z != 0)
        {
          if (!make_single_charged)
          {
            in.push_back(old_spectrum[i]);
          }
          else
          {
            RichPeak1D p = old_spectrum[i];
            p.setMZ(p.getMZ() * z - (z - 1) * Constants::PROTON_MASS_U);
            in.push_back(p);
          }
        }
      }
    }

    in.sortByPosition();
  }

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "Input file with MS/MS spectra", false);
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerInputFile_("id", "<file>", "", "Identification input file which contains a search against a concatenated sequence database", false);
    setValidFormats_("id", ListUtils::create<String>("idXML"));
    registerOutputFile_("out", "<file>", "", "Identification output with annotated phosphorylation scores");
    registerDoubleOption_("fragment_mass_tolerance", "<tolerance>", 0.05, "Fragment mass error", false);

    StringList fragment_mass_tolerance_unit_valid_strings;
    fragment_mass_tolerance_unit_valid_strings.push_back("Da");
    fragment_mass_tolerance_unit_valid_strings.push_back("ppm");
    registerStringOption_("fragment_mass_unit", "<unit>", "Da", "Unit of fragment mass error", false, false);
    setValidStrings_("fragment_mass_unit", fragment_mass_tolerance_unit_valid_strings);

    addEmptyLine_();
  }

  ExitCodes main_(int, const char**)
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    String in(getStringOption_("in"));
    String id(getStringOption_("id"));
    String out(getStringOption_("out"));
    double fragment_mass_tolerance(getDoubleOption_("fragment_mass_tolerance"));
    bool fragment_mass_unit_ppm = getStringOption_("fragment_mass_unit") == "Da" ? false : true;
    AScore ascore;

    RichPeakMap exp;
    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------
    vector<PeptideIdentification> pep_ids;
    vector<ProteinIdentification> prot_ids;
    vector<PeptideIdentification> pep_out;
    IdXMLFile().load(id, prot_ids, pep_ids);

    MzMLFile f;
    f.setLogType(log_type_);

    PeakFileOptions options;
    options.clearMSLevels();
    options.addMSLevel(2);
    f.getOptions() = options;
    f.load(in, exp);
    exp.sortSpectra(true);

    for (Size i = 0; i != exp.getNrSpectra(); ++i)
    {
      deisotopeAndSingleChargeMSSpectrum_(exp.getSpectrum(i), 2, 5, fragment_mass_tolerance, fragment_mass_unit_ppm);
    }

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    // map the ids to the spectra
    IDMapper id_mapper;
    id_mapper.annotate(exp, pep_ids, prot_ids);
    for (RichPeakMap::iterator it = exp.begin(); it != exp.end(); ++it)
    {
      if (it->getPeptideIdentifications().empty())
      {
        continue;
      }

      for (vector<PeptideIdentification>::iterator hits = it->getPeptideIdentifications().begin(); hits < it->getPeptideIdentifications().end(); ++hits)
      {
        vector<PeptideHit> scored_peptides;
        for (vector<PeptideHit>::const_iterator hit = hits->getHits().begin(); hit < hits->getHits().end(); ++hit)
        {
          PeptideHit scored_hit = *hit;
          RichPeakSpectrum& temp  = *it;

          PeptideHit phospho_sites;
          phospho_sites = ascore.compute(scored_hit, temp, fragment_mass_tolerance, fragment_mass_unit_ppm);
          scored_peptides.push_back(phospho_sites);
        }

        PeptideIdentification new_hits(*hits);
        new_hits.setScoreType("PhosphoScore");
        new_hits.setHits(scored_peptides);
        pep_out.push_back(new_hits);
      }
    }
    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    IdXMLFile().store(out, prot_ids, pep_out);
    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPPhosphoScoring tool;
  return tool.main(argc, argv);
}

/// @endcond
