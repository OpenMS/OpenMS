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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignment.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/MATH/STATISTICS/Histogram.h>
#include <OpenMS/MATH/STATISTICS/GaussFitter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

using namespace OpenMS;
using namespace std;
using namespace Math;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_IDMassAccuracy IDMassAccuracy

    @brief Calculates a distribution of the mass error from given mass spectra and IDs.

    @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_IDMassAccuracy.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_IDMassAccuracy.html

    Given a number of peak maps and for each of the maps an idXML file which contains
    peptide identifications the theoretical masses of the identifications and the peaks
    of the spectra are compared. This can be done for precursor information stored in
    the spectra as well as for fragment information.

    The result is a distribution of errors of experimental vs. theoretical masses.
    Having such distributions given
    the search parameters of the sequence database search can be adjusted to speed-up
    the identification process and to get a higher performance.
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

// simple struct which can hold the
// measured and expected masses
struct MassDifference
{
  double exp_mz = 0.0;
  Int charge = 0;
  double theo_mz = 0.0;
  double intensity = 0.0;
};

class TOPPIDMassAccuracy :
  public TOPPBase
{
public:
  TOPPIDMassAccuracy() :
    TOPPBase("IDMassAccuracy", "Calculates a distribution of the mass error from given mass spectra and IDs.", false)
  {

  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<file list>", StringList(), "Input mzML file list, containing the spectra.");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerInputFileList_("id_in", "<file list>", StringList(), "Input idXML file list, containing the identifications.");
    setValidFormats_("id_in", ListUtils::create<String>("idXML"));

    registerOutputFile_("precursor_out", "<file>", "", "Output file which contains the deviations from the precursors", false, false);
    setValidFormats_("precursor_out", ListUtils::create<String>("csv"));

    registerStringList_("precursor_columns", "<columns>", ListUtils::create<String>("MassDifference"), "Columns which will be written to the output file", false);
    setValidStrings_("precursor_columns", ListUtils::create<String>("MassDifference"));
    registerFlag_("precursor_error_ppm", "If this flag is used, the precursor mass tolerances are estimated in ppm instead of Da.");

    registerOutputFile_("fragment_out", "<file>", "", "Output file which contains the fragment ion m/z deviations", false, false);
    setValidFormats_("fragment_out", ListUtils::create<String>("csv"));
    registerStringList_("fragment_columns", "<columns>", ListUtils::create<String>("MassDifference"), "Columns which will be written to the output file", false);
    setValidStrings_("fragment_columns", ListUtils::create<String>("MassDifference"));
    registerFlag_("fragment_error_ppm", "If this flag is used, the fragment mass tolerances are estimated in ppm instead of Da.");

    registerDoubleOption_("fragment_mass_tolerance", "<tolerance>", 0.5, "Maximal fragment mass tolerance which is allowed for MS/MS spectra, used for the calculation of matching ions.", false, false);
    registerStringOption_("separator", "<character>", " ", "character which should be used to separate the columns in the output files", false);

    registerIntOption_("number_of_bins", "<#bins>", 100, "Number of bins that should be used to calculate the histograms for the fitting.", false, true);
    setMinInt_("number_of_bins", 10);

    registerStringOption_("generate_gnuplot_scripts", "<false>", "false", "If this option is set to true, the distributions and the fits are used to generate a gnuplot script, that can be used to generate plots. The options 'precursor_out' and 'fragment_out' must be set to take this effect.", false, true);
    setValidStrings_("generate_gnuplot_scripts", ListUtils::create<String>("true,false"));
  }

  double getMassDifference(double theo_mz, double exp_mz, bool use_ppm)
  {
    double error(exp_mz - theo_mz);
    if (use_ppm)
    {
      error = error / theo_mz * (double)1e6;
    }
    return error;
  }

  ExitCodes main_(int, const char **) override
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------

    StringList id_in(getStringList_("id_in"));
    StringList in_raw(getStringList_("in"));
    Size number_of_bins((UInt)getIntOption_("number_of_bins"));
    bool precursor_error_ppm(getFlag_("precursor_error_ppm"));
    bool fragment_error_ppm(getFlag_("fragment_error_ppm"));
    bool generate_gnuplot_scripts(DataValue(getStringOption_("generate_gnuplot_scripts")).toBool());

    if (in_raw.size() != id_in.size())
    {
      writeLog_("Number of spectrum files and identification files differs...");
      return ILLEGAL_PARAMETERS;
    }

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------

    vector<vector<PeptideIdentification> > pep_ids;
    vector<vector<ProteinIdentification> > prot_ids;
    pep_ids.resize(id_in.size());
    prot_ids.resize(id_in.size());

    IdXMLFile idxmlfile;
    for (Size i = 0; i != id_in.size(); ++i)
    {
      String doc_id;
      idxmlfile.load(id_in[i], prot_ids[i], pep_ids[i], doc_id);
    }

    // read mzML files
    vector<PeakMap> maps_raw;
    maps_raw.resize(in_raw.size());

    MzMLFile mzml_file;
    for (Size i = 0; i != in_raw.size(); ++i)
    {
      mzml_file.load(in_raw[i], maps_raw[i]);
    }

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    // mapping ids
    IDMapper mapper;
    for (Size i = 0; i != maps_raw.size(); ++i)
    {
      mapper.annotate(maps_raw[i], pep_ids[i], prot_ids[i]);
    }

    // normalize the spectra
    Normalizer normalizer;
    for (vector<PeakMap>::iterator it1 = maps_raw.begin(); it1 != maps_raw.end(); ++it1)
    {
      for (PeakMap::Iterator it2 = it1->begin(); it2 != it1->end(); ++it2)
      {
        normalizer.filterSpectrum(*it2);
      }
    }

    // generate precursor statistics
    vector<MassDifference> precursor_diffs;
    if (getStringOption_("precursor_out") != "")
    {
      for (Size i = 0; i != maps_raw.size(); ++i)
      {
        for (Size j = 0; j != maps_raw[i].size(); ++j)
        {
          if (maps_raw[i][j].getPeptideIdentifications().empty())
          {
            continue;
          }
          for (vector<PeptideIdentification>::const_iterator it = maps_raw[i][j].getPeptideIdentifications().begin(); it != maps_raw[i][j].getPeptideIdentifications().end(); ++it)
          {
            if (it->getHits().size() > 0)
            {
              PeptideHit hit = *it->getHits().begin();
              MassDifference md;
              Int charge = hit.getCharge();
              if (charge == 0)
              {
                charge = 1;
              }
              md.exp_mz = it->getMZ();
              md.theo_mz = (hit.getSequence().getMonoWeight() + (double)charge * Constants::PROTON_MASS_U) / (double)charge;
              md.charge = charge;
              precursor_diffs.push_back(md);
            }
          }
        }
      }
    }

    // generate fragment ions statistics
    vector<MassDifference> fragment_diffs;
    TheoreticalSpectrumGenerator tsg;
    SpectrumAlignment sa;
    double fragment_mass_tolerance(getDoubleOption_("fragment_mass_tolerance"));
    Param sa_param(sa.getParameters());
    sa_param.setValue("tolerance", fragment_mass_tolerance);
    sa.setParameters(sa_param);

    if (getStringOption_("fragment_out") != "")
    {
      for (Size i = 0; i != maps_raw.size(); ++i)
      {
        for (Size j = 0; j != maps_raw[i].size(); ++j)
        {
          if (maps_raw[i][j].getPeptideIdentifications().empty())
          {
            continue;
          }
          for (vector<PeptideIdentification>::const_iterator it = maps_raw[i][j].getPeptideIdentifications().begin(); it != maps_raw[i][j].getPeptideIdentifications().end(); ++it)
          {
            if (it->getHits().size() > 0)
            {
              PeptideHit hit = *it->getHits().begin();

              PeakSpectrum theo_spec;
              tsg.getSpectrum(theo_spec, hit.getSequence(), 1, 1);

              vector<pair<Size, Size> > pairs;
              sa.getSpectrumAlignment(pairs, theo_spec, maps_raw[i][j]);
              //cerr << hit.getSequence() << " " << hit.getSequence().getSuffix(1).getFormula() << " " << hit.getSequence().getSuffix(1).getFormula().getMonoWeight() << endl;
              for (vector<pair<Size, Size> >::const_iterator pit = pairs.begin(); pit != pairs.end(); ++pit)
              {
                MassDifference md;
                md.exp_mz = maps_raw[i][j][pit->second].getMZ();
                md.theo_mz = theo_spec[pit->first].getMZ();
                //cerr.precision(15);
                //cerr << md.exp_mz << " " << md.theo_mz << " " << md.exp_mz - md.theo_mz << endl;
                md.intensity = maps_raw[i][j][pit->second].getIntensity();
                md.charge = hit.getCharge();
                fragment_diffs.push_back(md);
              }
            }
          }
        }
      }
    }

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    String precursor_out_file(getStringOption_("precursor_out"));
    if (precursor_out_file != "")
    {
      vector<double> errors;
      ofstream precursor_out(precursor_out_file.c_str());
      double min_diff(numeric_limits<double>::max()), max_diff(numeric_limits<double>::min());
      for (Size i = 0; i != precursor_diffs.size(); ++i)
      {
        double diff = getMassDifference(precursor_diffs[i].theo_mz, precursor_diffs[i].exp_mz, precursor_error_ppm);
        precursor_out << diff << "\n";
        errors.push_back(diff);

        if (diff > max_diff)
        {
          max_diff = diff;
        }
        if (diff < min_diff)
        {
          min_diff = diff;
        }
      }
      precursor_out.close();

      // fill histogram with the collected values
      double bin_size = (max_diff - min_diff) / (double)number_of_bins;
      Histogram<double, double> hist(min_diff, max_diff, bin_size);
      for (Size i = 0; i != errors.size(); ++i)
      {
        hist.inc(errors[i], 1.0);
      }

      writeDebug_("min_diff=" + String(min_diff) + ", max_diff=" + String(max_diff) + ", number_of_bins=" + String(number_of_bins), 1);

      // transform the histogram into a vector<DPosition<2> > for the fitting
      vector<DPosition<2> > values;
      for (Size i = 0; i != hist.size(); ++i)
      {
        DPosition<2> p;
        p.setX((double)i / (double)number_of_bins * (max_diff - min_diff) + min_diff);
        p.setY(hist[i]);
        values.push_back(p);
      }

      double mean = Math::mean(errors.begin(), errors.end());
      double abs_dev = Math::absdev(errors.begin(), errors.end(), mean);
      double sdv = Math::sd(errors.begin(), errors.end(), mean);
      sort(errors.begin(), errors.end());
      double median = errors[(Size)(errors.size() / 2.0)];

      writeDebug_("Precursor mean error: " + String(mean), 1);
      writeDebug_("Precursor abs. dev.:  " + String(abs_dev), 1);
      writeDebug_("Precursor std. dev.:  " + String(sdv), 1);
      writeDebug_("Precursor median error:  " + String(median), 1);


      // calculate histogram for gauss fitting
      GaussFitter gf;
      GaussFitter::GaussFitResult init_param (hist.maxValue(), median, sdv/500.0);
      gf.setInitialParameters(init_param);

      try
      {
        gf.fit(values);

        // write gnuplot scripts
        if (generate_gnuplot_scripts)
        {
          ofstream out(String(precursor_out_file + "_gnuplot.dat").c_str());
          for (vector<DPosition<2> >::const_iterator it = values.begin(); it != values.end(); ++it)
          {
            out << it->getX() << " " << it->getY() << endl;
          }
          out.close();

          ofstream gpl_out(String(precursor_out_file + "_gnuplot.gpl").c_str());
          gpl_out << "set terminal png" << endl;
          gpl_out << "set output \"" << precursor_out_file  << "_gnuplot.png\"" << endl;
          if (precursor_error_ppm)
          {
            gpl_out << "set xlabel \"error in ppm\"" << endl;
          }
          else
          {
            gpl_out << "set xlabel \"error in Da\"" << endl;
          }
          gpl_out << "set ylabel \"frequency\"" << endl;
          gpl_out << "plot '" << precursor_out_file << "_gnuplot.dat' title 'Precursor mass error distribution' w boxes, f(x) w lp title 'Gaussian fit of the error distribution'" << endl;
          gpl_out.close();
        }

      }
      catch (Exception::UnableToFit)
      {
        writeLog_("Unable to fit a Gaussian distribution to the precursor mass errors");
      }
    }

    String fragment_out_file(getStringOption_("fragment_out"));
    if (fragment_out_file != "")
    {
      vector<double> errors;
      ofstream fragment_out(fragment_out_file.c_str());
      double min_diff(numeric_limits<double>::max()), max_diff(numeric_limits<double>::min());
      for (Size i = 0; i != fragment_diffs.size(); ++i)
      {
        double diff = getMassDifference(fragment_diffs[i].theo_mz, fragment_diffs[i].exp_mz, fragment_error_ppm);
        fragment_out << diff << endl;
        errors.push_back(diff);

        if (diff > max_diff)
        {
          max_diff = diff;
        }
        if (diff < min_diff)
        {
          min_diff = diff;
        }
      }
      fragment_out.close();

      // fill histogram with the collected values
      // here we use the intensities to scale the error
      // low intensity peaks are likely to be random matches
      double bin_size = (max_diff - min_diff) / (double)number_of_bins;
      Histogram<double, double> hist(min_diff, max_diff, bin_size);
      for (Size i = 0; i != fragment_diffs.size(); ++i)
      {
        double diff = getMassDifference(fragment_diffs[i].theo_mz, fragment_diffs[i].exp_mz, fragment_error_ppm);
        hist.inc(diff, fragment_diffs[i].intensity);
      }

      writeDebug_("min_diff=" + String(min_diff) + ", max_diff=" + String(max_diff) + ", number_of_bins=" + String(number_of_bins), 1);

      // transform the histogram into a vector<DPosition<2> > for the fitting
      vector<DPosition<2> > values;
      for (Size i = 0; i != hist.size(); ++i)
      {
        DPosition<2> p;
        p.setX((double)i / (double)number_of_bins * (max_diff - min_diff) + min_diff);
        p.setY(hist[i]);
        values.push_back(p);
      }

      double mean = Math::mean(errors.begin(), errors.end());
      double abs_dev = Math::absdev(errors.begin(), errors.end(), mean);
      double sdv = Math::sd(errors.begin(), errors.end(), mean);
      sort(errors.begin(), errors.end());
      double median = errors[(Size)(errors.size() / 2.0)];

      writeDebug_("Fragment mean error:  " + String(mean), 1);
      writeDebug_("Fragment abs. dev.:   " + String(abs_dev), 1);
      writeDebug_("Fragment std. dev.:   " + String(sdv), 1);
      writeDebug_("Fragment median error:   " + String(median), 1);

      // calculate histogram for gauss fitting
      GaussFitter gf;
      GaussFitter::GaussFitResult init_param (hist.maxValue(), median, sdv / 100.0);
      gf.setInitialParameters(init_param);

      try
      {
        gf.fit(values);


        // write gnuplot script
        if (generate_gnuplot_scripts)
        {
          ofstream out(String(fragment_out_file + "_gnuplot.dat").c_str());
          for (vector<DPosition<2> >::const_iterator it = values.begin(); it != values.end(); ++it)
          {
            out << it->getX() << " " << it->getY() << endl;
          }
          out.close();

          ofstream gpl_out(String(fragment_out_file + "_gnuplot.gpl").c_str());
          gpl_out << "set terminal png" << endl;
          gpl_out << "set output \"" << fragment_out_file  << "_gnuplot.png\"" << endl;
          if (fragment_error_ppm)
          {
            gpl_out << "set xlabel \"error in ppm\"" << endl;
          }
          else
          {
            gpl_out << "set xlabel \"error in Da\"" << endl;
          }
          gpl_out << "set ylabel \"frequency\"" << endl;
          gpl_out << "plot '" << fragment_out_file << "_gnuplot.dat' title 'Fragment mass error distribution' w boxes, f(x) w lp title 'Gaussian fit of the error distribution'" << endl;
          gpl_out.close();
        }
      }
      catch (Exception::UnableToFit)
      {
        writeLog_("Unable to fit a Gaussian distribution to the fragment mass errors");
      }
    }

    return EXECUTION_OK;
  }

};


int main(int argc, const char ** argv)
{
  TOPPIDMassAccuracy tool;
  return tool.main(argc, argv);
}

/// @endcond
