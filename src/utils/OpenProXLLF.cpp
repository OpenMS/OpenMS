
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
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/XLMS/OpenProXLUtils.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/CHEMISTRY/EnzymesDB.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/ANALYSIS/RNPXL/ModifiedPeptideGenerator.h>
#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/FORMAT/Base64.h>
#include <OpenMS/ANALYSIS/ID/PeptideIndexing.h>

// TESTING SCORES
#include <OpenMS/ANALYSIS/RNPXL/HyperScore.h>
#include <OpenMS/ANALYSIS/RNPXL/PScore.h>

// preprocessing and filtering
#include <OpenMS/FILTERING/TRANSFORMERS/ThresholdMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/NLargest.h>
#include <OpenMS/FILTERING/TRANSFORMERS/WindowMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignment.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGeneratorXLinks.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGeneratorXLMS.h>

// results
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <iostream>
#include <cmath>
#include <numeric>

#include <boost/unordered_map.hpp>
#include <boost/math/special_functions/binomial.hpp>

#ifdef _OPENMP
#include <omp.h>
#define NUMBER_OF_THREADS (omp_get_num_threads())
#else
#define NUMBER_OF_THREADS (1)
#endif

using namespace std;
using namespace OpenMS;

/**
    @page UTILS_OpenProXLLF OpenProXLLF

    @brief Perform protein-protein cross-linking experiment search.

    <CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ OpenProXLLF \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> - </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> - </td>
        </tr>
    </table>
</CENTER>

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_OpenProXLLF.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_OpenProXLLF.html
*/

class TOPPOpenProXLLF :
  public TOPPBase
{
public:
  TOPPOpenProXLLF() :
    TOPPBase("OpenProXLLF", "Tool for protein-protein cross linking with label-free linkers.", false)
  {
  }

protected:
  static void wrap_(const String& input, Size width, String & output)
  {
    Size start = 0;

    while (start + width < input.size())
    {
      output += input.substr(start, width) + "\n";
      start += width;
    }

    if (start < input.size())
    {
      output += input.substr(start, input.size() - start) + "\n";
    }
  }

protected:
  static String getxQuestBase64EncodedSpectrum_(const PeakSpectrum& spec, String header)
  {
    vector<String> in_strings;
    StringList sl;

    double precursor_mz = spec.getPrecursors()[0].getMZ();
    double precursor_z = spec.getPrecursors()[0].getCharge();

    // header lines
    if (!header.empty()) // common or xlinker spectrum will be reported
    {
      sl.push_back(header + "\n"); // e.g. GUA1372-S14-A-LRRK2_DSS_1A3.03873.03873.3.dta,GUA1372-S14-A-LRRK2_DSS_1A3.03863.03863.3.dta
      sl.push_back(String(precursor_mz) + "\n");
      sl.push_back(String(precursor_z) + "\n");
    }
    else // light or heavy spectrum will be reported
    {
      sl.push_back(String(precursor_mz) + "\t" + String(precursor_z) + "\n");
    }

    // write peaks
    for (Size i = 0; i != spec.size(); ++i)
    {
      String s;
      s += String(spec[i].getMZ()) + "\t";
      s += String(spec[i].getIntensity()) + "\t";

      // add fragment charge if meta value exists (must be present for 'common' and 'xlinker'.
//      if (spec[i].metaValueExists("z"))
//      {
//        s += String(spec[i].getMetaValue("z"));
//      }
      s += "0";

      s += "\n";

      sl.push_back(s);
    }

    String out;
    out.concatenate(sl.begin(), sl.end(), "");
    in_strings.push_back(out);

    String out_encoded;
    Base64().encodeStrings(in_strings, out_encoded, false, false);
    String out_wrapped;
    wrap_(out_encoded, 76, out_wrapped);
    return out_wrapped;
  }

protected:
  void registerOptionsAndFlags_()
  {
    // input files
    registerInputFile_("in", "<file>", "", "Input file containing the spectra.");
    setValidFormats_("in", ListUtils::create<String>("mzML"));

    //registerInputFile_("consensus", "<file>", "", "Input file containing the linked mass peaks.");
    //setValidFormats_("consensus", ListUtils::create<String>("consensusXML"));

    registerInputFile_("database", "<file>", "", "Input file containing the protein database.");
    setValidFormats_("database", ListUtils::create<String>("fasta"));

    registerInputFile_("decoy_database", "<file>", "", "Input file containing the decoy protein database.", false);
    setValidFormats_("decoy_database", ListUtils::create<String>("fasta"));

    registerStringOption_("decoy_string", "<string>", "decoy", "String that was appended (or prefixed - see 'prefix' flag below) to the accessions in the protein database to indicate decoy proteins.", false);
    registerFlag_("decoy_prefix", "Set flag, if the decoy_string is a prefix of accessions in the protein database. Otherwise it is a suffix.");

    registerTOPPSubsection_("precursor", "Precursor (Parent Ion) Options");
    registerDoubleOption_("precursor:mass_tolerance", "<tolerance>", 10.0, "Width of precursor mass tolerance window", false);

    StringList precursor_mass_tolerance_unit_valid_strings;
    precursor_mass_tolerance_unit_valid_strings.push_back("ppm");
    precursor_mass_tolerance_unit_valid_strings.push_back("Da");

    registerStringOption_("precursor:mass_tolerance_unit", "<unit>", "ppm", "Unit of precursor mass tolerance.", false, false);
    setValidStrings_("precursor:mass_tolerance_unit", precursor_mass_tolerance_unit_valid_strings);

    registerIntOption_("precursor:min_charge", "<num>", 3, "Minimum precursor charge to be considered.", false, true);
    registerIntOption_("precursor:max_charge", "<num>", 7, "Maximum precursor charge to be considered.", false, true);

    registerTOPPSubsection_("fragment", "Fragments (Product Ion) Options");
    registerDoubleOption_("fragment:mass_tolerance", "<tolerance>", 3, "Fragment mass tolerance", false);
    registerDoubleOption_("fragment:mass_tolerance_xlinks", "<tolerance>", 3, "Fragment mass tolerance for cross-link ions", false);

    StringList fragment_mass_tolerance_unit_valid_strings;
    fragment_mass_tolerance_unit_valid_strings.push_back("ppm");
    fragment_mass_tolerance_unit_valid_strings.push_back("Da");

    registerStringOption_("fragment:mass_tolerance_unit", "<unit>", "ppm", "Unit of fragment m", false, false);
    setValidStrings_("fragment:mass_tolerance_unit", fragment_mass_tolerance_unit_valid_strings);

    registerTOPPSubsection_("modifications", "Modifications Options");
    vector<String> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
    registerStringList_("modifications:fixed", "<mods>", ListUtils::create<String>(""), "Fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)'", false);
    setValidStrings_("modifications:fixed", all_mods);
    registerStringList_("modifications:variable", "<mods>", ListUtils::create<String>(""), "Variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Oxidation (M)'", false);
    setValidStrings_("modifications:variable", all_mods);
    registerIntOption_("modifications:variable_max_per_peptide", "<num>", 2, "Maximum number of residues carrying a variable modification per candidate peptide", false, false);

    registerTOPPSubsection_("peptide", "Peptide Options");
    registerIntOption_("peptide:min_size", "<num>", 5, "Minimum size a peptide must have after digestion to be considered in the search.", false, true);
    registerIntOption_("peptide:missed_cleavages", "<num>", 2, "Number of missed cleavages.", false, false);
    vector<String> all_enzymes;
    EnzymesDB::getInstance()->getAllNames(all_enzymes);
    registerStringOption_("peptide:enzyme", "<cleavage site>", "Trypsin", "The enzyme used for peptide digestion.", false);
    setValidStrings_("peptide:enzyme", all_enzymes);


    registerTOPPSubsection_("cross_linker", "Cross Linker Options");
    registerStringList_("cross_linker:residue1", "<one letter code>", ListUtils::create<String>("K"), "Comma separated residues, that the first side of a bifunctional cross-linker can attach to", false);
    registerStringList_("cross_linker:residue2", "<one letter code>", ListUtils::create<String>("K"), "Comma separated residues, that the second side of a bifunctional cross-linker can attach to", false);
    registerDoubleOption_("cross_linker:mass", "<mass>", 138.0680796, "Mass of the light cross-linker, linking two residues on one or two peptides", false);
    //registerDoubleOption_("cross_linker:mass_iso_shift", "<mass>", 12.075321, "Mass of the isotopic shift between the light and heavy linkers", false);
    registerDoubleList_("cross_linker:mass_mono_link", "<mass>", ListUtils::create<double>("156.07864431, 155.094628715"), "Possible masses of the linker, when attached to only one peptide", false);
    registerStringOption_("cross_linker:name", "<string>", "DSS" ,  "Name of the searched cross-link, used to resolve ambiguity of equal masses (e.g. DSS or BS3)", false);

    registerTOPPSubsection_("algorithm", "Algorithm Options");
    registerStringOption_("algorithm:candidate_search", "<param>", "enumeration", "Mode used to generate candidate peptides.", false, false);
    StringList candidate_search_modes_strings;
    candidate_search_modes_strings.push_back("index");
    candidate_search_modes_strings.push_back("enumeration");
    setValidStrings_("algorithm:candidate_search", candidate_search_modes_strings);

    registerIntOption_("algorithm:number_top_hits", "<num>", 5, "Number of top hits reported for each spectrum pair", false, true);

    // output file
    registerOutputFile_("out_xquestxml", "<file>", "", "Results in the original xquest.xml format", false);
    setValidFormats_("out_xquestxml", ListUtils::create<String>("xml"));

    registerOutputFile_("out_idXML", "<file>", "", "Results in idXML format", false);
    setValidFormats_("out_idXML", ListUtils::create<String>("idXML"));

    registerOutputFile_("out_mzIdentML", "<file>","", "Results in mzIdentML (.mzid) format", false);
    setValidFormats_("out_mzIdentML", ListUtils::create<String>("mzid"));
  }

  vector<ResidueModification> getModifications_(StringList modNames)
  {
    vector<ResidueModification> modifications;

    // iterate over modification names and add to vector
    for (StringList::iterator mod_it = modNames.begin(); mod_it != modNames.end(); ++mod_it)
    {
      String modification(*mod_it);
      modifications.push_back(ModificationsDB::getInstance()->getModification(modification));
    }

    return modifications;
  }

  // check if for minimum size
  class HasInvalidPeptideLengthPredicate
  {
      public:
        explicit HasInvalidPeptideLengthPredicate(Size min_size)
          :min_size_(min_size)
        {
        }

        bool operator()(const AASequence& aas)
        {
          return (aas.size() < min_size_);
        }
    private:
        Size min_size_;
  };

protected:
  static double map_add (double value, const std::map<Size, double>::value_type& p)
  {
    return value + p.second;
  }

//// Smallest range with percent% number of peaks
//protected:
//  static pair<Size, Size> find_cover_region_peaknumber(PeakSpectrum spectrum, double percent)
//  {
//    Size n = spectrum.size();

//    Size k = floor((1-percent) * n);
//    vector<int> a_vec = boost::copy_range<std::vector<int>>(boost::irange(n-k, n));
//    vector<int> b_vec = boost::copy_range<std::vector<int>>(boost::irange(0, k+1));
//    //Size i = which.min(x[seq.int(n-k, n)] - x[seq_len(k+1L)])
//    vector<Size> possible_pos;

//    for (Size i = 0; i < a_vec.size(); ++i)
//    {
//      possible_pos.push_back(a_vec[i] - b_vec[i]);
//    }
//    vector::iterator min_pos_it= min_element(possible_pos.begin(), possible_pos.end());
//    Size min_pos = min_pos_it - possible_pos.begin();

//    pair<Size, Size> interval = make_pair(x[min_pos], x[n-k+min_pos-1]);
//    return interval;
//  }

// Smallest range with percent% total peak intensity
protected:
  double find_cover_region(PeakSpectrum spectrum, double percent, double total_intensity)
  {
    double threshold = percent * total_intensity;

    Size pos1 = 0, pos2 = 0;
    double mzrange = INFINITY;

    for (Size k = 0; k < spectrum.size(); ++k)
    {
      for (Size n = 0; n < k; ++n)
      {
        double int_sum = 0;
        for (Size i = n; i <= k; ++i)
        {
          int_sum += spectrum[i].getIntensity();
        }
        if (int_sum >= threshold)
        {
          //pos1 = n;
          //pos2 = k;
          double current_range = spectrum[k].getMZ() - spectrum[n].getMZ();
          if (current_range < mzrange)
          {
            mzrange = current_range;
            pos1 = n;
            pos2 = k;
          }
        }
      }
    }

    return mzrange;
  }

protected:
  PeakMap preprocessSpectra_(PeakMap& exp)
  {

    Size peptide_min_size = getIntOption_("peptide:min_size") * 2; //  x2 because cross-links
    Int min_precursor_charge = getIntOption_("precursor:min_charge");
    Int max_precursor_charge = getIntOption_("precursor:max_charge");

    //String precursor_mass_tolerance_unit = getStringOption_("precursor:mass_tolerance_unit");
    String fragment_mass_tolerance_unit = getStringOption_("fragment:mass_tolerance_unit");
    bool fragment_mass_tolerance_ppm = fragment_mass_tolerance_unit == "ppm";
    //double fragment_mass_tolerance = getDoubleOption_("fragment:mass_tolerance");
    double fragment_mass_tolerance_xlinks = getDoubleOption_("fragment:mass_tolerance_xlinks");

    // filter MS2 map
    // remove 0 intensities
    ThresholdMower threshold_mower_filter;
    threshold_mower_filter.filterPeakMap(exp);
    // TODO perl code filters by dynamic range (1000), meaning everything below max_intensity / 1000 is filtered out additionally to 0 int, before scaling / normalizing

    Normalizer normalizer;
    normalizer.filterPeakMap(exp);
    // TODO perl code scales to 0-100: int / max_int * 100

    // sort by rt
    exp.sortSpectra(false);

    // filter settings
//    WindowMower window_mower_filter;
//    Param filter_param = window_mower_filter.getParameters();
//    filter_param.setValue("windowsize", 100.0, "The size of the sliding window along the m/z axis.");
//    filter_param.setValue("peakcount", 20, "The number of peaks that should be kept.");
//    filter_param.setValue("movetype", "jump", "Whether sliding window (one peak steps or jumping window window size steps) should be used.");
//    window_mower_filter.setParameters(filter_param);
    //NLargest nlargest_filter = NLargest(500);   // De-noising in xQuest: Dynamic range = 1000, 250 most intense peaks?
    LOG_DEBUG << "Deisotoping and filtering spectra." << endl;


    // EXPERIMENTAL SPECTRUM QUALITY FEATURES (turned out to be useless, at least in the implemenetd form)
//    map<Size, double> peak_number_map;
//    map<Size, double> std_var_map;
//    map<Size, double> avg_neighbor_map;
//    map<Size, double> current_per_mz_map;
//    map<Size, double> mz_diff_stdvar_map;
//    map<Size, double> peak_number_with_charge_map;
//    map<Size, double> cover_region_95_map;
//    map<Size, double> cover_region_50_map;
//    map<Size, double> spectrum_quality_score_map;

    // Calculation of weighted score
    // weights, in the order of the vectors above:
//    double a = 0.3, b = 0.05, c = -0.005, d = 0.05, e = -0.03, f = 0.3; // This scoring with threshold > 3 = about 7k spectra, all pLink analyzed included

//    double a = -0.063, b = -0.024, c = 0.348, d = 0.432, e = -0.382, f = 0.959; // Scoring from LDA, with 57 / 60 class1 spectra (3 bad pLink examples), about 5500 spectra

//    double threshold_multiplier = 1.5;

//#ifdef _OPENMP
//#pragma omp parallel for
//#endif
//    for (SignedSize exp_index = 0; exp_index < static_cast<SignedSize>(exp.size()); ++exp_index)
//    {
//#ifdef _OPENMP
//#pragma omp critical
//#endif
//      cout << "Preprocessing spectrum " << exp_index << " of " << exp.size() << " ; Spectrum ID: " << exp[exp_index].getNativeID() << " MS Level: " << exp[exp_index].getMSLevel() << endl;

      // there should only one precursor and MS2 should contain at least a few peaks to be considered (e.g. at least for every AA in the peptide)
//      vector<Precursor> precursor = exp[exp_index].getPrecursors();
//      bool process_this_spectrum = false;
//      if (precursor.size() == 1 && exp[exp_index].size() >= peptide_min_size*2)
//      {
//        int precursor_charge = precursor[0].getCharge();
//        if (precursor_charge >= min_precursor_charge && precursor_charge <= max_precursor_charge)
//        {
//          process_this_spectrum = true;
//        }
//      }



//#ifdef _OPENMP
//#pragma omp critical
//#endif
//      peak_number_with_charge_map.insert(make_pair(exp_index, 0));
//      int counter = 0;

//      if (!process_this_spectrum)
//      {
//#ifdef _OPENMP
//#pragma omp critical
//#endif
//        {
//          peak_number_map.insert(make_pair(exp_index, 0));
//          std_var_map.insert(make_pair(exp_index, 0));
//          avg_neighbor_map.insert(make_pair(exp_index, 0));
//          current_per_mz_map.insert(make_pair(exp_index, 0));
//          mz_diff_stdvar_map.insert(make_pair(exp_index, 0));
//          cover_region_95_map.insert(make_pair(exp_index, 0));
//          cover_region_50_map.insert(make_pair(exp_index, 0));
//          spectrum_quality_score_map.insert(make_pair(exp_index, -999));
//        }
//      }
//      else
//      {

//        counter++;
//        // sort by mz
//        exp[exp_index].sortByPosition();

        // TODO more sophisticated spectra filter, maybe similar to pLink?
        // General Spectrum Quality features from: Nesvizhskii, A. I.  et al . Mol Cell Proteomics 5 , 652 -­‐ 670 (2006)

        //Number of peaks, square root
//        double peak_number = sqrt(exp[exp_index].size());

        // Arithmetic mean of peak intensities, log-transformed
//        double total_current = 0;

        // Standard deviation of the peak intensities, log-transformed
//        double std_var = 0;

        // Average number of neighbor peaks within a 2-Da interval around any peak
//        double avg_neighbors = 0;

//        double mean = 0;
//        for (Size i = 0; i < exp[exp_index].size(); ++i)
//        {
//          total_current += exp[exp_index][i].getIntensity();
//          PeakSpectrum neighbors = getToleranceWindowPeaks(exp[exp_index], exp[exp_index][i].getMZ(), 2.0, false);
//          avg_neighbors += neighbors.size();
//        }
//        mean = log(total_current / exp[exp_index].size());
//        avg_neighbors = avg_neighbors / exp[exp_index].size();

//        for (Size i = 0; i < exp[exp_index].size(); ++i)
//        {
//          std_var += pow(mean - exp[exp_index][i].getIntensity(), 2);
//        }
//        std_var = log(sqrt(std_var / exp[exp_index].size()));

        // Total ion current per m/z (total ion current divided by feature d), log-transformed
//        double feature_d = exp[exp_index][exp[exp_index].size()-1].getMZ() - exp[exp_index][0].getMZ();
//        double current_per_mz = -log(total_current / feature_d);

        // Standard deviation of the consecutive m/z gaps between all peaks, log-transformed
//        double mz_diff_stdvar = 0;
//        double mz_diff_mean = 0;

//        for (Size i = 0; i < exp[exp_index].size()-1; ++i)
//        {
//          mz_diff_mean += exp[exp_index][i+1].getMZ() - exp[exp_index][i].getMZ();
//        }
//        mz_diff_mean = mz_diff_mean / (exp[exp_index].size()-1);

//        for (Size i = 0; i < exp[exp_index].size()-1; ++i)
//        {
//          mz_diff_stdvar += pow(mz_diff_mean - (exp[exp_index][i+1].getMZ() - exp[exp_index][i].getMZ()), 2);
//        }

//        mz_diff_stdvar = log(sqrt(mz_diff_stdvar / (exp[exp_index].size()-1)));

        // COMPUTATIONALLY VERY EXPENSIVE, BUT SEEM TO BE USELESS, maybe try another combo, e.g. distance between them, number of peaks inside?

//        // Smallest m/z range containing 95% of the total peak number             // TODO % of total peak intensity
//        double cover_region_95 = find_cover_region(exp[exp_index], 0.95, total_current);


//        // Smallest m/z range containing 50% of the total peak intensity
//        double cover_region_50 =  find_cover_region(exp[exp_index], 0.5, total_current);

        // Calculation of weighted score
//        double spectrum_quality_score = a * peak_number + b * std_var + c * avg_neighbors + d * current_per_mz + e * mz_diff_stdvar;

//#ifdef _OPENMP
//#pragma omp critical
//#endif
//        {
//          peak_number_map.insert(make_pair(exp_index, peak_number));
//          std_var_map.insert(make_pair(exp_index, std_var));
//          avg_neighbor_map.insert(make_pair(exp_index, avg_neighbors));
//          current_per_mz_map.insert(make_pair(exp_index, current_per_mz));
//          mz_diff_stdvar_map.insert(make_pair(exp_index, mz_diff_stdvar));
////          cover_region_95_map.insert(make_pair(exp_index, cover_region_95));
////          cover_region_50_map.insert(make_pair(exp_index, cover_region_50));
//          cover_region_95_map.insert(make_pair(exp_index, 0));
//          cover_region_50_map.insert(make_pair(exp_index, 0));
//          spectrum_quality_score_map.insert(make_pair(exp_index, spectrum_quality_score));
//          //cout << "Preprocessing spectrum " << exp_index << " of " << exp.size() << " ; Spectrum ID: " << exp[exp_index].getNativeID() << " MS Level: " << exp[exp_index].getMSLevel() << endl;
//          //cout << "Spectrum quality metrics: \nPeak number: " << peak_number << "\tstd_var: " << std_var << "\tavg_neighbors: " << avg_neighbors << "\tcurrent_per_mz: " << current_per_mz << "\tmz_diff_stdvar: " << mz_diff_stdvar << "\tSQS: " << spectrum_quality_score << endl;
//        }
//      }
//    }


//    double peak_number_mean = accumulate(peak_number_map.begin(), peak_number_map.end(), 0.0, map_add) / peak_number_map.size();
//    double std_var_mean = accumulate(std_var_map.begin(), std_var_map.end(), 0.0, map_add) / std_var_map.size();
//    double avg_neighbors_mean = accumulate(avg_neighbor_map.begin(), avg_neighbor_map.end(), 0.0, map_add) / avg_neighbor_map.size();
//    double current_per_mz_mean = accumulate(current_per_mz_map.begin(), current_per_mz_map.end(), 0.0, map_add) / current_per_mz_map.size();
//    double mz_diff_stdvar_mean = accumulate(mz_diff_stdvar_map.begin(), mz_diff_stdvar_map.end(), 0.0, map_add) / mz_diff_stdvar_map.size();
//    double spectrum_quality_score_mean = a * peak_number_mean + b * std_var_mean + c * avg_neighbors_mean + d * current_per_mz_mean + e * mz_diff_stdvar_mean;
////    cout << "Spectrum quality metrics means: \n" << "peak_number_mean" << "\t" << "std_var_mean" << "\t" << "avg_neighbor_mean" << "\t" << "current_per_mz_mean" << "\t" << "mz_diff_stdvar_mean" << "\t" << "average QS" << "\t" << "threshold\n"
////                                                                              << peak_number_mean << "\t\t\t" << std_var_mean << "\t\t" << avg_neighbors_mean << "\t\t\t" << current_per_mz_mean << "\t\t\t" << mz_diff_stdvar_mean << "\t\t\t" << spectrum_quality_score_mean << "\t\t\t" << spectrum_quality_score_mean * threshold_multiplier  << endl;

//    double passed_threshold = 0;

    PeakMap deisotoped_spectra;

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (SignedSize exp_index = 0; exp_index < static_cast<SignedSize>(exp.size()); ++exp_index)
    {
      vector<Precursor> precursor = exp[exp_index].getPrecursors();
      bool process_this_spectrum = false;
      if (precursor.size() == 1 && exp[exp_index].size() >= peptide_min_size*2)
      {
        int precursor_charge = precursor[0].getCharge();
        if (precursor_charge >= min_precursor_charge && precursor_charge <= max_precursor_charge)
        {
          process_this_spectrum = true;
        }
      }

      if (!process_this_spectrum)
      {
        continue;
      }
      exp[exp_index].sortByPosition();

      // TODO think of a better boundary
      // if ( spectrum_quality_score_map.find(exp_index)->second > (spectrum_quality_score_mean * threshold_multiplier) )
//      if(true)
//      {

//#ifdef _OPENMP
//#pragma omp atomic
//#endif
//      passed_threshold++;


//#ifdef _OPENMP
//#pragma omp critical
//#endif
//      cout << "Spectrum passed threshold: " << exp[exp_index].getNativeID() << "\twith SQS: " << spectrum_quality_score_map.find(exp_index)->second << " > " << (spectrum_quality_score_mean * threshold_multiplier)  << endl;

        // sort by mz (done above when calculating quality scores already)
        // exp[exp_index].sortByPosition();
        // params:                                                                                             (PeakSpectrum,  min_charge, max_charge,  fragment_tol,                                 tol_unit_ppm,                              only_keep_deiso, min_iso_peaks, max_iso_peaks, make_single_charged);
        PeakSpectrum deisotoped = deisotopeAndSingleChargeMSSpectrum(exp[exp_index], 1,                 7,                    10,                                                 fragment_mass_tolerance_ppm, false,                    3,                      10,                     false);

        //check, whether at least one charge is > 0, e.g. sum over the charge array > 0
//        int charge_sum = 0;
//        cout << "Deisotoped Spectrum size: " << deisotoped.size() << endl;
//        if (deisotoped.size() > peptide_min_size * 2)
//        {
//          if (deisotoped.getIntegerDataArrays().size() > 0)
//          {
//            PeakSpectrum::IntegerDataArray charge_array = deisotoped.getIntegerDataArrays()[0];
//            charge_sum = accumulate(charge_array.begin(), charge_array.end(), 0);
//          }
//        }


//        if (deisotoped.size() > peptide_min_size)
//        {
//#ifdef _OPENMP
//#pragma omp critical
//#endif
//          double charge_peaks = 0.0;
//          for (Size i = 0; i < deisotoped.size(); ++i)
//          {
//            if( double(deisotoped[i].getMetaValue("z")) != 0)
//            {
//              charge_peaks += 1.0;
//            }
//          }
//          peak_number_with_charge_map[exp_index] = sqrt(charge_peaks);
//          double new_score = spectrum_quality_score_map.find(exp_index)->second;
//          new_score += f * sqrt(charge_peaks);
//          spectrum_quality_score_map[exp_index] = new_score;

//          passed_threshold++;
          //cout << "Deisotoped a spectrum, original size: " << exp[exp_index].size() << ", new size: " << deisotoped.size() << endl;
//          deisotoped_spectra[exp_index] = deisotoped;

          if (deisotoped.size() > peptide_min_size * 2)
          {
            // TODO filterSpectrum() has to remove DataArray Entries too
//            nlargest_filter.filterSpectrum(deisotoped);
            OpenProXLUtils::nLargestSpectrumFilter(deisotoped, 500);
            deisotoped.sortByPosition();

#ifdef _OPENMP
#pragma omp critical
#endif
            deisotoped_spectra.addSpectrum(deisotoped);
          }

//        }
//        else
//        {
////          deisotoped_spectra[exp_index] = PeakSpectrum();

//          PeakSpectrum dummy = PeakSpectrum();
//          dummy.setPrecursors(exp[exp_index].getPrecursors());
//          dummy.setRT(exp[exp_index].getRT());
//          dummy.setNativeID(exp[exp_index].getNativeID());
//          dummy.setInstrumentSettings(exp[exp_index].getInstrumentSettings());
//          dummy.setAcquisitionInfo(exp[exp_index].getAcquisitionInfo());
//          dummy.setSourceFile(exp[exp_index].getSourceFile());
//          dummy.setDataProcessing(exp[exp_index].getDataProcessing());
//          dummy.setType(exp[exp_index].getType());
//          dummy.setMSLevel(exp[exp_index].getMSLevel());
//          dummy.setName(exp[exp_index].getName());

//#ifdef _OPENMP
//#pragma omp critical
//#endif
//          deisotoped_spectra.addSpectrum(dummy);
//          cout << "Spectrum did not pass: " << exp[exp_index].getNativeID() << "\twith SQS: " << spectrum_quality_score_map.find(exp_index)->second << " > " << (spectrum_quality_score_mean * threshold_multiplier) << " , kicked by deisotoping."  << endl;
//        }

//      }
//      else
//      {

//        PeakSpectrum dummy = PeakSpectrum();
//        dummy.setPrecursors(exp[exp_index].getPrecursors());
//        dummy.setRT(exp[exp_index].getRT());
//        dummy.setNativeID(exp[exp_index].getNativeID());
//        dummy.setInstrumentSettings(exp[exp_index].getInstrumentSettings());
//        dummy.setAcquisitionInfo(exp[exp_index].getAcquisitionInfo());
//        dummy.setSourceFile(exp[exp_index].getSourceFile());
//        dummy.setDataProcessing(exp[exp_index].getDataProcessing());
//        dummy.setType(exp[exp_index].getType());
//        dummy.setMSLevel(exp[exp_index].getMSLevel());
//        dummy.setName(exp[exp_index].getName());

//#ifdef _OPENMP
//#pragma omp critical
//#endif
//        deisotoped_spectra.addSpectrum(dummy);
//        cout << "Spectrum did not pass: " << exp[exp_index].getNativeID() << "\twith SQS: " << spectrum_quality_score_map.find(exp_index)->second << " < " << (spectrum_quality_score_mean * threshold_multiplier)  << endl;
//      }

    }

    //cout << "Number of Spectra passed_threshold (and deisotoping): " << passed_threshold << endl;


    //############### TEST CODE (filtering scores) ######################
//    for (Size exp_index = 0; exp_index < peak_number_map.size(); ++exp_index)
//    {

//      String sqs_filename = "SQ_Scores.txt";
//      ofstream sqs_file;
//      sqs_file.open(sqs_filename.c_str(), ios::app);

//      sqs_file << exp[exp_index].getNativeID() << " " << peak_number_map.find(exp_index)->second << " " << std_var_map.find(exp_index)->second << " " << avg_neighbor_map.find(exp_index)->second << " "
//                                                                                  << current_per_mz_map.find(exp_index)->second << " " << mz_diff_stdvar_map.find(exp_index)->second << " "
//                                                                                   << peak_number_with_charge_map.find(exp_index)->second << " "  << cover_region_95_map.find(exp_index)->second  << " "
//                                                                                   << cover_region_50_map.find(exp_index)->second << " " << spectrum_quality_score_map.find(exp_index)->second << endl;


//      }
      //############### TEST CODE END ####################

    // TODO another loop to determine which spectra to throw out based on final SQS, right now all spectra are considered
//    PeakMap sqs_filtered_spectra;
//    for (Size exp_index = 0; exp_index < spectrum_quality_score_map.size(); ++exp_index)
//    {
//      // TODO VERIFY HARDCODED THRESHOLD AND SCORING SCHEME (maybe threshold s advanced parameter?)
//      // Threshold based on GUA0001_2 :  1.2
//      if (spectrum_quality_score_map.find(exp_index)->second > 0)
//      {
//        //window_mower_filter.filterPeakSpectrum(deisotoped_spectra[exp_index]);
//        nlargest_filter.filterSpectrum(deisotoped_spectra[exp_index]);
//        // sort (nlargest changes order)
//        deisotoped_spectra[exp_index].sortByPosition();
//        sqs_filtered_spectra.addSpectrum(deisotoped_spectra[exp_index]);
//      }
//    }

    //cout << "Spectra left after filtering: " << sqs_filtered_spectra.size() << " / " << deisotoped_spectra.size() << endl;

//    return sqs_filtered_spectra;
    return deisotoped_spectra;
  }

  // Transform a PeakSpectrum into a RichPeakSpectrum
//  RichPeakSpectrum makeRichPeakSpectrum(PeakSpectrum spectrum, bool is_common_or_xlink_spectrum)
//  {
//    RichPeakSpectrum rich_spectrum;

//    rich_spectrum.setPrecursors(spectrum.getPrecursors());
//    rich_spectrum.setRT(spectrum.getRT());

//    for(Size i = 0; i < spectrum.size(); ++i)
//    {
//      Peak1D peak = spectrum[i];
//      RichPeak1D p;
//      p.setMZ(peak.getMZ());
//      p.setIntensity(peak.getIntensity());
//      if (is_common_or_xlink_spectrum)
//      {
//        p.setMetaValue("z", 0); // TODO Where to get the actual charge?
//      }
//      rich_spectrum.push_back(p);
//    }
//    return rich_spectrum;
//  }

    PeakSpectrum deisotopeAndSingleChargeMSSpectrum(PeakSpectrum& old_spectrum, Int min_charge, Int max_charge, double fragment_tolerance, bool fragment_tolerance_unit_ppm, bool keep_only_deisotoped = false, Size min_isopeaks = 3, Size max_isopeaks = 10, bool make_single_charged = false)
    {

      // Input Spectrum originally called "in"
      //PeakSpectrum old_spectrum = in;
      PeakSpectrum out;
      PeakSpectrum::IntegerDataArray charge_array;

      vector<Size> mono_isotopic_peak(old_spectrum.size(), 0);
      if (old_spectrum.empty())
      {
        return out;
      }

      // determine charge seeds and extend them
      vector<Int> features(old_spectrum.size(), -1);
      Int feature_number = 0;

      for (Size current_peak = 0; current_peak != old_spectrum.size(); ++current_peak)
      {
        double current_mz = old_spectrum[current_peak].getMZ();

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
              double tolerance_dalton = fragment_tolerance_unit_ppm ? fragment_tolerance * old_spectrum[p].getMZ() * 1e-6 : fragment_tolerance;
              if (fabs(old_spectrum[p].getMZ() - expected_mz) > tolerance_dalton)   // test for missing peak
              {
                if (i < min_isopeaks)
                {
                  has_min_isopeaks = false;
                }
                break;
              }
              else
              {
                // TODO: include proper averagine model filtering. assuming the intensity gets lower for heavier peaks does not work for the high masses of cross-linked peptides
//                Size n_extensions = extensions.size();
//                if (n_extensions != 0)
//                {
//                  if (old_spectrum[p].getIntensity() > old_spectrum[extensions[n_extensions - 1]].getIntensity())
//                  {
//                    if (i < min_isopeaks)
//                    {
//                      has_min_isopeaks = false;
//                    }
//                    break;
//                  }
//                }

                // averagine check passed
                extensions.push_back(p);
              }
            }

            if (has_min_isopeaks)
            {
              //LOG_DEBUG << "min peaks at " << current_mz << " " << " extensions: " << extensions.size() << endl;
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


      // creating PeakSpectrum containing charges
      //out.clear(false);

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
            RichPeak1D p;
            p.setMZ(old_spectrum[i].getMZ());
            p.setIntensity(old_spectrum[i].getIntensity());
//            p.setMetaValue("z", z);
            charge_array.push_back(z);
            out.push_back(p);
          }
          else
          {
            RichPeak1D p;
            p.setIntensity(old_spectrum[i].getIntensity());
            p.setMZ(old_spectrum[i].getMZ() * z - (z - 1) * Constants::PROTON_MASS_U);
//            p.setMetaValue("z", 1);
            charge_array.push_back(1);
            out.push_back(p);
          }
        }
        else
        {
          // keep all unassigned peaks
          if (features[i] < 0)
          {
            RichPeak1D p;
            p.setMZ(old_spectrum[i].getMZ());
            p.setIntensity(old_spectrum[i].getIntensity());
//            p.setMetaValue("z", 0);
            charge_array.push_back(0);
            out.push_back(p);
            continue;
          }

          // convert mono-isotopic peak with charge assigned by deisotoping
          if (z != 0)
          {
            if (!make_single_charged)
            {
              RichPeak1D p;
              p.setMZ(old_spectrum[i].getMZ());
              p.setIntensity(old_spectrum[i].getIntensity());
//              p.setMetaValue("z", z);
              charge_array.push_back(z);
              out.push_back(p);
            }
            else
            {
              RichPeak1D p;
              p.setMZ(old_spectrum[i].getMZ() * z - (z - 1) * Constants::PROTON_MASS_U);
              p.setIntensity(old_spectrum[i].getIntensity());
//              p.setMetaValue("z", z);
              charge_array.push_back(z);
              out.push_back(p);
            }
          }
        }
      }
      out.setPrecursors(old_spectrum.getPrecursors());
      out.setRT(old_spectrum.getRT());

      out.setNativeID(old_spectrum.getNativeID());
      out.setInstrumentSettings(old_spectrum.getInstrumentSettings());
      out.setAcquisitionInfo(old_spectrum.getAcquisitionInfo());
      out.setSourceFile(old_spectrum.getSourceFile());
      out.setDataProcessing(old_spectrum.getDataProcessing());
      out.setType(old_spectrum.getType());
      out.setMSLevel(old_spectrum.getMSLevel());
      out.setName(old_spectrum.getName());

      out.getIntegerDataArrays().push_back(charge_array);

//      out.sortByPosition();
      return out;
    }

    // Spectrum Alignment function adapted from SpectrumAlignment.h, intensity_cutoff: 0 for not considering, 0.3 = lower intensity has to be at least 30% of higher intensity (< 70% difference)
    //template <typename SpectrumType1, typename SpectrumType2>
    void getSpectrumAlignment(std::vector<std::pair<Size, Size> > & alignment, const PeakSpectrum & s1, PeakSpectrum s2, double tolerance, bool relative_tolerance, double intensity_cutoff = 0.0) const
    {
      if (!s1.isSorted() || !s2.isSorted())
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Input to SpectrumAlignment is not sorted!");
        }

      // clear result
      alignment.clear();
      //double tolerance = (double)param_.getValue("tolerance");

      if (!(relative_tolerance))
        {
          std::map<Size, std::map<Size, std::pair<Size, Size> > > traceback;
          std::map<Size, std::map<Size, double> > matrix;

          // init the matrix with "gap costs" tolerance
          matrix[0][0] = 0;
          for (Size i = 1; i <= s1.size(); ++i)
            {
              matrix[i][0] = i * tolerance;
              traceback[i][0]  = std::make_pair(i - 1, 0);
            }
          for (Size j = 1; j <= s2.size(); ++j)
            {
              matrix[0][j] = j * tolerance;
              traceback[0][j] = std::make_pair(0, j - 1);
            }

          // fill in the matrix
          Size left_ptr(1);
          Size last_i(0), last_j(0);

          //Size off_band_counter(0);
          for (Size i = 1; i <= s1.size(); ++i)
            {
              double pos1(s1[i - 1].getMZ());

              for (Size j = left_ptr; j <= s2.size(); ++j)
                {
                  bool off_band(false);
                  // find min of the three possible directions
                  double pos2(s2[j - 1].getMZ());
                  double diff_align = fabs(pos1 - pos2);

                  // running off the right border of the band?
                  if (pos2 > pos1 && diff_align >= tolerance)
                    {
                      if (i < s1.size() && j < s2.size() && s1[i].getMZ() < pos2)
                        {
                          off_band = true;
                        }
                    }

                  // can we tighten the left border of the band?
                  if (pos1 > pos2 && diff_align >= tolerance && j > left_ptr + 1)
                    {
                      ++left_ptr;
                    }

                  double score_align = diff_align;

                  if (matrix.find(i - 1) != matrix.end() && matrix[i - 1].find(j - 1) != matrix[i - 1].end())
                    {
                      score_align += matrix[i - 1][j - 1];
                    }
                  else
                    {
                      score_align += (i - 1 + j - 1) * tolerance;
                    }

                  double score_up = tolerance;
                  if (matrix.find(i) != matrix.end() && matrix[i].find(j - 1) != matrix[i].end())
                    {
                      score_up += matrix[i][j - 1];
                    }
                  else
                    {
                      score_up += (i + j - 1) * tolerance;
                    }

                  double score_left = tolerance;
                  if (matrix.find(i - 1) != matrix.end() && matrix[i - 1].find(j) != matrix[i - 1].end())
                    {
                      score_left += matrix[i - 1][j];
                    }
                  else
                    {
                      score_left += (i - 1 + j) * tolerance;
                    }

                  // check for similar intensity values
                  double intensity1(s1[i - 1].getIntensity());
                  double intensity2(s2[j - 1].getIntensity());
                  bool diff_int_clear = (min(intensity1, intensity2) / max(intensity1, intensity2)) > intensity_cutoff;

                  // check for same charge
//                  PeakSpectrum::IntegerDataArray s1_charges = s1.getIntegerDataArrays()[0];
//                  PeakSpectrum::IntegerDataArray s2_charges = s2.getIntegerDataArrays()[0];
//                  bool charge_fits = s1_charges[i - 1] == s2_charges[j - 1] || s2_charges[j - 1] == 0;
//                  LOG_DEBUG << "s1 charge: " << s1_charges[i - 1] << " | s2 charge: " << s2_charges[j - 1] << endl;

                  // TODO SET A CHARGE FOR EXPERIMENTAL SPECTRA, then use this again
                  bool charge_fits = true;

//                  LOG_DEBUG << "CHARGE FIT TEST 1### charge1 = " << charge1 << " ### charge2 = " << charge2 << " ### charge fits = " << charge_fits << endl;


                  if (score_align <= score_up && score_align <= score_left && diff_align < tolerance && diff_int_clear && charge_fits)
                    {
                      matrix[i][j] = score_align;
                      traceback[i][j] = std::make_pair(i - 1, j - 1);
                      last_i = i;
                      last_j = j;

//#ifdef _OPENMP
//#pragma omp critical
//#endif
//                      cout << "TEST aligned peaks: " << pos1 << " = " << pos2 << "\t| tolerance: " << tolerance << "\t| diff_align: " << diff_align << "\t| ion name: " << s1[i-1].getMetaValue("IonName") << endl;
                    }
                  else
                    {
                      if (score_up <= score_left)
                        {
                          matrix[i][j] = score_up;
                          traceback[i][j] = std::make_pair(i, j - 1);
                        }
                      else
                        {
                          matrix[i][j] = score_left;
                          traceback[i][j] = std::make_pair(i - 1, j);
                        }
                    }

                  if (off_band)
                    {
                      break;
                    }
                }
            }

          // do traceback
          Size i = last_i;
          Size j = last_j;

          while (i >= 1 && j >= 1)
            {
              if (traceback[i][j].first == i - 1 && traceback[i][j].second == j - 1)
                {
                  alignment.push_back(std::make_pair(i - 1, j - 1));
                }
              Size new_i = traceback[i][j].first;
              Size new_j = traceback[i][j].second;

              i = new_i;
              j = new_j;
            }

          std::reverse(alignment.begin(), alignment.end());

        }
      else  // relative alignment (ppm tolerance)
        {
          for (Size i = 0; i != s1.size(); ++i)
            {
              const double& theo_mz = s1[i].getMZ();
              double max_dist_dalton = theo_mz * tolerance * 1e-6;

              // iterate over peaks in experimental spectrum in given fragment tolerance around theoretical peak
              Size j = s2.findNearest(theo_mz);
              double exp_mz = s2[j].getMZ();

              // check for similar intensity values
              double intensity1(s1[i].getIntensity());
              double intensity2(s2[j].getIntensity());
              bool diff_int_clear = (min(intensity1, intensity2) / max(intensity1, intensity2)) > intensity_cutoff;

              // check for same charge
              PeakSpectrum::IntegerDataArray s1_charges = s1.getIntegerDataArrays()[0];
              PeakSpectrum::IntegerDataArray s2_charges = s2.getIntegerDataArrays()[0];
              bool charge_fits = s1_charges[i - 1] == s2_charges[j - 1] || s2_charges[j - 1] == 0;

//              LOG_DEBUG << "CHARGE FIT TEST 2### charge1 = " << charge1 << " ### charge2 = " << charge2 << " ### charge fits = " << charge_fits << endl;

              // found peak match
              if (std::abs(theo_mz - exp_mz) <= max_dist_dalton && diff_int_clear && charge_fits)
                {
                  alignment.push_back(std::make_pair(i, j));
//#ifdef _OPENMP
//#pragma omp critical
//#endif
//                  cout << "TEST aligned peaks: " << theo_mz << " = " << exp_mz << "\t| Da tolerance: " << max_dist_dalton << "\t| diff_align: " << std::abs(theo_mz - exp_mz) << "\t| ion name: " << s1[i].getMetaValue("IonName") << endl;
                } else if (std::abs(theo_mz - exp_mz) <= max_dist_dalton && !diff_int_clear)
                {
                  s2.erase(s2.begin() + s2.findNearest(theo_mz));
//#ifdef _OPENMP
//#pragma omp critical
//#endif
//                  cout << "TEST peak was jumped, because of intensity check" << endl;
                  if (i > 0)
                  {
                    i--;
                   }
                }
            }
        }
    }

    PeakSpectrum getToleranceWindowPeaks(PeakSpectrum spec, double mz, double tolerance, bool relative_tolerance) const
    {
      PeakSpectrum peaks;
      double max_dist;

      if (relative_tolerance)
      {
        max_dist = mz + (mz * tolerance * 1e-6);
      } else
      {
        max_dist = tolerance;
      }

      bool inside = true;
      while (inside)
      {
        Size index = spec.findNearest(mz);
        double peak_mz = spec[index].getMZ();
        double dist = abs(mz - peak_mz);

        if (dist <= max_dist)
        {
          peaks.push_back(spec[index]);
          spec.erase(spec.begin() + index);
        } else
        {
          inside = false;
        }
      }
      return peaks;
    }

    // Write xQuest.xml output
//    void writeXQuestXML(const String& out_file, String base_name, const vector< PeptideIdentification >& peptide_ids, const vector< vector< CrossLinkSpectrumMatch > >& all_top_csms, const RichPeakMap& spectra)
//    {
//      String spec_xml_name = base_name + "_matched";

//      cout << "Writing xquest.xml to " << out_file << endl;
//      ofstream xml_file;
//      xml_file.open(out_file.c_str(), ios::trunc);
//      // XML Header
//      xml_file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
//      xml_file << "<?xml-stylesheet type=\"text/xsl\" href=\"\"?>" << endl;

//      // TODO!!! write actual experiment data
//      // original date/time: Fri Dec 18 12:28:23 2015
//      DateTime time= DateTime::now();
//      String timestring = time.getDate() + " " + time.getTime();

//      String precursor_mass_tolerance_unit = getStringOption_("precursor:mass_tolerance_unit");
//      String fragment_mass_tolerance_unit = getStringOption_("fragment:mass_tolerance_unit");
//      double precursor_mass_tolerance = getDoubleOption_("precursor:mass_tolerance");
//      double fragment_mass_tolerance = getDoubleOption_("fragment:mass_tolerance");
//      double fragment_mass_tolerance_xlinks = getDoubleOption_("fragment:mass_tolerance_xlinks");

//      String cross_link_name = getStringOption_("cross_linker:name");
//      double cross_link_mass = getDoubleOption_("cross_linker:mass");
//      DoubleList cross_link_mass_mono_link = getDoubleList_("cross_linker:mass_mono_link");
//      String mono_masses;
//      for (Size k = 0; k < cross_link_mass_mono_link.size()-1; ++k)
//      {
//        mono_masses += String(cross_link_mass_mono_link[k]) + ", ";
//      }
//      mono_masses += cross_link_mass_mono_link[cross_link_mass_mono_link.size()-1];

//      const string in_fasta(getStringOption_("database"));
//      const string in_decoy_fasta(getStringOption_("decoy_database"));
//      StringList cross_link_residue1 = getStringList_("cross_linker:residue1");
//      StringList cross_link_residue2 = getStringList_("cross_linker:residue2");
//      String aarequired1, aarequired2;
//      for (Size k= 0; k < cross_link_residue1.size()-1; ++k)
//      {
//        aarequired1 += cross_link_residue1[k] + ",";
//      }
//      aarequired1 += cross_link_residue1[cross_link_residue1.size()-1];
//      for (Size k= 0; k < cross_link_residue2.size()-1; ++k)
//      {
//        aarequired2 += cross_link_residue2[k] + ",";
//      }
//      aarequired2 += cross_link_residue2[cross_link_residue2.size()-1];

//      //double cross_link_mass_iso_shift = getDoubleOption_("cross_linker:mass_iso_shift");
//      String enzyme_name = getStringOption_("peptide:enzyme");
//      Size missed_cleavages = getIntOption_("peptide:missed_cleavages");

//      xml_file << "<xquest_results xquest_version=\"OpenProXL 1.0\" date=\"" << timestring <<
//             "\" author=\"Eugen Netz, Timo Sachsenberg\" tolerancemeasure_ms1=\"" << precursor_mass_tolerance_unit  <<
//             "\" tolerancemeasure_ms2=\"" << fragment_mass_tolerance_unit << "\" ms1tolerance=\"" << precursor_mass_tolerance <<
//             "\" ms2tolerance=\"" << fragment_mass_tolerance << "\" xlink_ms2tolerance=\"" << fragment_mass_tolerance_xlinks <<
//             "\" crosslinkername=\"" << cross_link_name << "\" xlinkermw=\"" << cross_link_mass <<
//             "\" monolinkmw=\"" << mono_masses << "\" database=\"" << in_fasta << "\" database_dc=\"" << in_decoy_fasta <<
//             "\" xlinktypes=\"1111\" AArequired1=\"" << aarequired1 << "\" AArequired2=\"" << aarequired2 <<  "\" cp_isotopediff=\"" << 0 <<
//             "\" enzyme_name=\"" << enzyme_name << "\" outputpath=\"" << spec_xml_name <<
//             "\" Iontag_charges_for_index=\"1\" missed_cleavages=\"" << missed_cleavages <<
//             "\" ntermxlinkable=\"0\" CID_match2ndisotope=\"1" <<
//             "\" variable_mod=\"TODO\" nocutatxlink=\"1\" xcorrdelay=\"5\" >" << endl;



//      for (vector< vector< CrossLinkSpectrumMatch > >::const_iterator top_csms_spectrum = all_top_csms.begin(); top_csms_spectrum != all_top_csms.end(); ++top_csms_spectrum)
//      {
//        vector< CrossLinkSpectrumMatch > top_vector = (*top_csms_spectrum);

//        if (!top_vector.empty())
//        {
//          // Spectrum Data, for each spectrum
//          Size scan_index = top_vector[0].scan_index_light;
//          //Size scan_index_heavy = top_vector[0].scan_index_heavy;
//          const PeakSpectrum& spectrum = spectra[scan_index];
//          double precursor_charge = spectrum.getPrecursors()[0].getCharge();

//          double precursor_mz = spectrum.getPrecursors()[0].getMZ();
//          double precursor_rt = spectrum.getRT();
//          double precursor_mass = precursor_mz * static_cast<double>(precursor_charge) - static_cast<double>(precursor_charge) * Constants::PROTON_MASS_U;

//          // print information about new peak to file (starts with <spectrum_search..., ends with </spectrum_search>
//          String spectrum_light_name = base_name + ".light." + scan_index;
//          String spectrum_heavy_name = base_name + ".heavy." + scan_index;

//          String spectrum_name = spectrum_light_name + "_" + spectrum_heavy_name;
//          String rt_scans = String(precursor_rt) + ":" + String(precursor_rt);
//          String mz_scans = String(precursor_mz) + ":" + String(precursor_mz);

//          // Mean ion intensity (light spectrum, TODO add heavy spectrum?)
//          double mean_intensity= 0;
//          for (SignedSize j = 0; j < static_cast<SignedSize>(spectrum.size()); ++j) mean_intensity += spectrum[j].getIntensity();
//          mean_intensity = mean_intensity / spectrum.size();

//          xml_file << "<spectrum_search spectrum=\"" << spectrum_name << "\" mean_ionintensity=\"" << mean_intensity << "\" ionintensity_stdev=\"" << "TODO" << "\" addedMass=\"" << "TODO" << "\" iontag_ncandidates=\"" << "TODO"
//            << "\"  apriori_pmatch_common=\"" << "TODO" << "\" apriori_pmatch_xlink=\"" << "TODO" << "\" ncommonions=\"" << "TODO" << "\" nxlinkions=\"" << "TODO" << "\" mz_precursor=\"" << precursor_mz
//            << "\" scantype=\"" << "light" << "\" charge_precursor=\"" << precursor_charge << "\" Mr_precursor=\"" << precursor_mass <<  "\" rtsecscans=\"" << rt_scans << "\" mzscans=\"" << mz_scans << "\" >" << endl;


//          for (vector< CrossLinkSpectrumMatch>::const_iterator top_csm = top_csms_spectrum->begin(); top_csm != top_csms_spectrum->end(); ++top_csm)
//          {
//            String xltype = "monolink";
//            String structure = top_csm->cross_link.alpha.toUnmodifiedString();
//            String letter_first = structure.substr(top_csm->cross_link.cross_link_position.first, 1);


//             // TODO track or otherwise find out, which kind of mono-link it was (if there are several possibilities for the weigths)
//            double weight = top_csm->cross_link.alpha.getMonoWeight() + top_csm->cross_link.cross_linker_mass;
////            bool is_monolink = (top_csm->cross_link.cross_link_position.second == -1);
//            int alpha_pos = top_csm->cross_link.cross_link_position.first + 1;
//            int beta_pos = top_csm->cross_link.cross_link_position.second + 1;

//            String topology = String("a") + alpha_pos;
//            String id = structure + String("-") + letter_first + alpha_pos + String("-") + static_cast<int>(top_csm->cross_link.cross_linker_mass);

//            if (top_csm->cross_link.getType() == ProteinProteinCrossLink::CROSS)
//            {
//              xltype = "xlink";
//              structure += "-" + top_csm->cross_link.beta.toUnmodifiedString();
//              topology += String("-b") + beta_pos;
//              weight += top_csm->cross_link.beta.getMonoWeight();
//              id = structure + "-" + topology;
//            }
//            else if (top_csm->cross_link.getType() == ProteinProteinCrossLink::LOOP)
//            {
//              xltype = "intralink";
//              topology += String("-b") + beta_pos;
//              String letter_second = structure.substr(top_csm->cross_link.cross_link_position.second, 1);
//              id = structure + String("-") + letter_first + alpha_pos + String("-") + letter_second + beta_pos;
//            }

//            // Error calculation
//            double cl_mz = (weight + (static_cast<double>(precursor_charge) * Constants::PROTON_MASS_U)) / static_cast<double>(precursor_charge);

//            double error = precursor_mz - cl_mz;
//            //double rel_error = (error / precursor_mz) / 1e-6;
//            double rel_error = (error / cl_mz) / 1e-6;
//            //cout << "rel_error[ppm]: " << rel_error << "\terror: " << error << "\tcl_mz: " << cl_mz << "\tweight: " << weight << endl;

//            PeptideIdentification pep_id = peptide_ids[top_csm->peptide_id_index];
//            vector< PeptideHit > pep_hits = pep_id.getHits();

//            String prot_alpha = pep_hits[0].getPeptideEvidences()[0].getProteinAccession();
//            if (pep_hits[0].getPeptideEvidences().size() > 1)
//            {
//              for (Size i = 1; i < pep_hits[0].getPeptideEvidences().size(); ++i)
//              {
//                prot_alpha = prot_alpha + "," + pep_hits[0].getPeptideEvidences()[i].getProteinAccession();
//              }
//            }

//            String prot_beta = "";

//            if (pep_hits.size() > 1)
//            {
//              prot_beta= pep_hits[1].getPeptideEvidences()[0].getProteinAccession();
//              if (pep_hits[1].getPeptideEvidences().size() > 1)
//              {
//                for (Size i = 1; i < pep_hits[1].getPeptideEvidences().size(); ++i)
//                {
//                  prot_alpha = prot_alpha + "," + pep_hits[1].getPeptideEvidences()[i].getProteinAccession();
//                }
//              }
//            }
//            // Hit Data, for each cross-link to Spectrum Hit (e.g. top 5 per spectrum)
//            xml_file << "<search_hit search_hit_rank=\"" <<top_csm->rank << "\" id=\"" << id << "\" type=\"" << xltype << "\" structure=\"" << structure << "\" seq1=\"" << top_csm->cross_link.alpha.toUnmodifiedString() << "\" seq2=\"" << top_csm->cross_link.beta.toUnmodifiedString()
//                << "\" prot1=\"" << prot_alpha << "\" prot2=\"" << prot_beta << "\" topology=\"" << topology << "\" xlinkposition=\"" << (top_csm->cross_link.cross_link_position.first+1) << "," << (top_csm->cross_link.cross_link_position.second+1)
//                << "\" Mr=\"" << weight << "\" mz=\"" << cl_mz << "\" charge=\"" << precursor_charge << "\" xlinkermass=\"" << top_csm->cross_link.cross_linker_mass << "\" measured_mass=\"" << precursor_mass << "\" error=\"" << error
//                << "\" error_rel=\"" << rel_error << "\" xlinkions_matched=\"" << (top_csm->matched_xlink_alpha + top_csm->matched_xlink_beta) << "\" backboneions_matched=\"" << (top_csm->matched_common_alpha + top_csm->matched_common_beta)
//                << "\" weighted_matchodds_mean=\"" << "TODO" << "\" weighted_matchodds_sum=\"" << "TODO" << "\" match_error_mean=\"" << "TODO" << "\" match_error_stdev=\"" << "TODO" << "\" xcorrx=\"" << top_csm->xcorrx_max << "\" xcorrb=\"" << top_csm->xcorrc_max << "\" match_odds=\"" <<top_csm->match_odds << "\" prescore=\"" << top_csm->pre_score
//                << "\" prescore_alpha=\"" << "TODO" << "\" prescore_beta=\"" << "TODO" << "\" match_odds_alphacommon=\"" << "TODO" << "\" match_odds_betacommon=\"" << "TODO" << "\" match_odds_alphaxlink=\"" << "TODO"
//                << "\" match_odds_betaxlink=\"" << "TODO" << "\" num_of_matched_ions_alpha=\"" << (top_csm->matched_common_alpha + top_csm->matched_xlink_alpha) << "\" num_of_matched_ions_beta=\"" << (top_csm->matched_common_beta + top_csm->matched_xlink_beta) << "\" num_of_matched_common_ions_alpha=\"" << top_csm->matched_common_alpha
//                << "\" num_of_matched_common_ions_beta=\"" << top_csm->matched_common_beta << "\" num_of_matched_xlink_ions_alpha=\"" << top_csm->matched_xlink_alpha << "\" num_of_matched_xlink_ions_beta=\"" << top_csm->matched_xlink_beta << "\" xcorrall=\"" << "TODO" << "\" TIC=\"" << top_csm->percTIC
//                << "\" TIC_alpha=\"" << "TODO" << "\" TIC_beta=\"" << "TODO" << "\" wTIC=\"" << top_csm->wTIC << "\" intsum=\"" << top_csm->int_sum * 100 << "\" apriori_match_probs=\"" << "TODO" << "\" apriori_match_probs_log=\"" << "TODO"
//                << "\" series_score_mean=\"" << "TODO" << "\" annotated_spec=\"" << "" << "\" score=\"" << top_csm->score << "\" >" << endl;
//            xml_file << "</search_hit>" << endl;
//          }
//          // Closing tag for Spectrum
//          xml_file << "</spectrum_search>" << endl;
//        }
//      }

//      // Closing tag for results (end of file)
//      xml_file << "</xquest_results>" << endl;
//      xml_file.close();

//      return;
//    }

    void writeXQuestXMLSpec(String out_file, String base_name, const vector< vector< CrossLinkSpectrumMatch > >& all_top_csms, const PeakMap& spectra)
    {
      // String spec_xml_filename = base_name + "_matched.spec.xml";
      // XML Header
      ofstream spec_xml_file;
      cout << "Writing spec.xml to " << out_file << endl;
      spec_xml_file.open(out_file.c_str(), ios::trunc); // ios::app = append to file, ios::trunc = overwrites file
      // TODO write actual data
      spec_xml_file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?><xquest_spectra compare_peaks_version=\"3.4\" date=\"Tue Nov 24 12:41:18 2015\" author=\"Thomas Walzthoeni,Oliver Rinner\" homepage=\"http://proteomics.ethz.ch\" resultdir=\"aleitner_M1012_004_matched\" deffile=\"xquest.def\" >" << endl;

      for (Size i = 0; i < spectra.size(); ++i)
      {
        if (!all_top_csms[i].empty())
        {
          String spectrum_light_name = base_name + ".light." + i;
          String spectrum_heavy_name = base_name + ".heavy." + i;

          String spectrum_name = spectrum_light_name + String("_") + spectrum_heavy_name;

          // 4 Spectra resulting from a light/heavy spectra pair.  Write for each spectrum, that is written to xquest.xml (should be all considered pairs, or better only those with at least one sensible Hit, meaning a score was computed)
          spec_xml_file << "<spectrum filename=\"" << spectrum_light_name << ".dta" << "\" type=\"light\">" << endl;
          spec_xml_file << TOPPOpenProXLLF::getxQuestBase64EncodedSpectrum_(spectra[i], String(""));
          spec_xml_file << "</spectrum>" << endl;

          spec_xml_file << "<spectrum filename=\"" << spectrum_heavy_name << ".dta" << "\" type=\"heavy\">" << endl;
          spec_xml_file << TOPPOpenProXLLF::getxQuestBase64EncodedSpectrum_(spectra[i], String(""));
          spec_xml_file << "</spectrum>" << endl;

          String spectrum_common_name = spectrum_name + String("_common.txt");
          spec_xml_file << "<spectrum filename=\"" << spectrum_common_name << "\" type=\"common\">" << endl;
          spec_xml_file << TOPPOpenProXLLF::getxQuestBase64EncodedSpectrum_(spectra[i], spectrum_light_name + ".dta," + spectrum_heavy_name + ".dta");
          spec_xml_file << "</spectrum>" << endl;

          String spectrum_xlink_name = spectrum_name + String("_xlinker.txt");
          spec_xml_file << "<spectrum filename=\"" << spectrum_xlink_name << "\" type=\"xlinker\">" << endl;
          spec_xml_file <<TOPPOpenProXLLF::getxQuestBase64EncodedSpectrum_(spectra[i], spectrum_light_name + ".dta," + spectrum_heavy_name + ".dta");
          spec_xml_file << "</spectrum>" << endl;
        }
      }

      spec_xml_file << "</xquest_spectra>" << endl;
      spec_xml_file.close();

      return;
    }

  ExitCodes main_(int, const char**)
  {
    ProgressLogger progresslogger;
    progresslogger.setLogType(log_type_);

    const string in_mzml(getStringOption_("in"));
    const string in_fasta(getStringOption_("database"));
    const string in_decoy_fasta(getStringOption_("decoy_database"));
    //const string in_consensus(getStringOption_("consensus"));
    const string out_idXML(getStringOption_("out_idXML"));
    const string out_xquest = getStringOption_("out_xquestxml");
    const string out_mzIdentML = getStringOption_("out_mzIdentML");

    const bool decoy_prefix(getFlag_("decoy_prefix"));
    const string decoy_string(getStringOption_("decoy_string"));

    Int min_precursor_charge = getIntOption_("precursor:min_charge");
    Int max_precursor_charge = getIntOption_("precursor:max_charge");
    double precursor_mass_tolerance = getDoubleOption_("precursor:mass_tolerance");
    bool precursor_mass_tolerance_unit_ppm = (getStringOption_("precursor:mass_tolerance_unit") == "ppm");

    double fragment_mass_tolerance = getDoubleOption_("fragment:mass_tolerance");
    double fragment_mass_tolerance_xlinks = getDoubleOption_("fragment:mass_tolerance_xlinks");
    bool fragment_mass_tolerance_unit_ppm = (getStringOption_("fragment:mass_tolerance_unit") == "ppm");

    SpectrumAlignment ms2_aligner;
    Param ms2_alinger_param = ms2_aligner.getParameters();
    String ms2_relative_tolerance = fragment_mass_tolerance_unit_ppm ? "true" : "false";
    ms2_alinger_param.setValue("is_relative_tolerance", ms2_relative_tolerance);
    ms2_alinger_param.setValue("tolerance", fragment_mass_tolerance);
    ms2_aligner.setParameters(ms2_alinger_param);

    SpectrumAlignment ms2_aligner_xlinks;
    Param ms2_aligner_xlinks_param = ms2_aligner_xlinks.getParameters();
    ms2_aligner_xlinks_param.setValue("is_relative_tolerance", ms2_relative_tolerance);
    ms2_aligner_xlinks_param.setValue("tolerance", fragment_mass_tolerance_xlinks);
    ms2_aligner_xlinks.setParameters(ms2_aligner_xlinks_param);

    StringList cross_link_residue1 = getStringList_("cross_linker:residue1");
    StringList cross_link_residue2 = getStringList_("cross_linker:residue2");
    double cross_link_mass = getDoubleOption_("cross_linker:mass");
    //double cross_link_mass_iso_shift = getDoubleOption_("cross_linker:mass_iso_shift");
    DoubleList cross_link_mass_mono_link = getDoubleList_("cross_linker:mass_mono_link");
    String cross_link_name = getStringOption_("cross_linker:name");

    StringList fixedModNames = getStringList_("modifications:fixed");
    set<String> fixed_unique(fixedModNames.begin(), fixedModNames.end());

    Size peptide_min_size = getIntOption_("peptide:min_size");

    bool ion_index_mode = (getStringOption_("algorithm:candidate_search") == "index");
    Int number_top_hits = getIntOption_("algorithm:number_top_hits");

    if (fixed_unique.size() != fixedModNames.size())
    {
      LOG_DEBUG << "duplicate fixed modification provided." << endl;
      return ILLEGAL_PARAMETERS;
    }

    StringList varModNames = getStringList_("modifications:variable");
    set<String> var_unique(varModNames.begin(), varModNames.end());
    if (var_unique.size() != varModNames.size())
    {
      LOG_DEBUG << "duplicate variable modification provided." << endl;
      return ILLEGAL_PARAMETERS;
    }
    //TODO add crosslinker

    vector<ResidueModification> fixed_modifications = getModifications_(fixedModNames);
    vector<ResidueModification> variable_modifications = getModifications_(varModNames);
    Size max_variable_mods_per_peptide = getIntOption_("modifications:variable_max_per_peptide");

    // load MS2 map
    PeakMap unprocessed_spectra;
    MzMLFile f;
    f.setLogType(log_type_);

    PeakFileOptions options;
    options.clearMSLevels();
    options.addMSLevel(2);
    f.getOptions() = options;
    f.load(in_mzml, unprocessed_spectra);
    unprocessed_spectra.sortSpectra(true);

    // preprocess spectra (filter out 0 values, sort by position)
    progresslogger.startProgress(0, 1, "Filtering spectra...");
    PeakMap spectra = preprocessSpectra_(unprocessed_spectra);
    progresslogger.endProgress();

    // load linked features
    ConsensusMap cfeatures;
    ConsensusXMLFile cf;
    //cf.load(in_consensus, cfeatures);

    // load fasta database
    progresslogger.startProgress(0, 1, "Load database from FASTA file...");
    FASTAFile fastaFile;
    vector<FASTAFile::FASTAEntry> fasta_db;
    fastaFile.load(in_fasta, fasta_db);


    if (!in_decoy_fasta.empty())
    {
      vector<FASTAFile::FASTAEntry> fasta_decoys;
      fastaFile.load(in_decoy_fasta, fasta_decoys);
      fasta_db.reserve(fasta_db.size() + fasta_decoys.size());
      fasta_db.insert(fasta_db.end(), fasta_decoys.begin(), fasta_decoys.end());
    }

    progresslogger.endProgress();

    const Size missed_cleavages = getIntOption_("peptide:missed_cleavages");
    EnzymaticDigestion digestor;
    String enzyme_name = getStringOption_("peptide:enzyme");
    digestor.setEnzyme(enzyme_name);
    digestor.setMissedCleavages(missed_cleavages);

    // lookup for processed peptides. must be defined outside of omp section and synchronized
    multimap<StringView, AASequence> processed_peptides;
    vector<OpenProXLUtils::PeptideMass> peptide_masses;

    // set minimum size of peptide after digestion
    Size min_peptide_length = getIntOption_("peptide:min_size");

    // one identification run
    vector<ProteinIdentification> protein_ids(1);
    protein_ids[0].setDateTime(DateTime::now());
    protein_ids[0].setSearchEngine("OpenXQuest");
    protein_ids[0].setSearchEngineVersion(VersionInfo::getVersion());
    protein_ids[0].setPrimaryMSRunPath(unprocessed_spectra.getPrimaryMSRunPath());
    protein_ids[0].setMetaValue("SpectrumIdentificationProtocol", DataValue("MS:1002494")); // cross-linking search = MS:1002494

    ProteinIdentification::SearchParameters search_params;
    search_params.charges = "2,3,4,5,6";
    search_params.db = in_fasta;
    search_params.digestion_enzyme = (*EnzymesDB::getInstance()->getEnzyme(enzyme_name));
    search_params.fixed_modifications = fixedModNames;
    search_params.variable_modifications = varModNames;
    search_params.mass_type = ProteinIdentification::MONOISOTOPIC;
    search_params.missed_cleavages = missed_cleavages;
    search_params.fragment_mass_tolerance = fragment_mass_tolerance;
    search_params.fragment_mass_tolerance_ppm =  fragment_mass_tolerance_unit_ppm;
    search_params.precursor_mass_tolerance = precursor_mass_tolerance;
    search_params.precursor_mass_tolerance_ppm = precursor_mass_tolerance_unit_ppm;

    // As MetaValues
    //search_params.setMetaValue("input_consensusXML", in_consensus);
    protein_ids[0].setMetaValue("input_mzML", in_mzml);
    //protein_ids[0].setMetaValue("input_database", in_fasta);

    search_params.setMetaValue("input_decoys", in_decoy_fasta);
    search_params.setMetaValue("decoy_prefix", decoy_prefix);
    search_params.setMetaValue("decoy_string", decoy_string);

    search_params.setMetaValue("precursor:min_charge", min_precursor_charge);
    search_params.setMetaValue("precursor:max_charge", max_precursor_charge);
    //protein_ids[0].setMetaValue("precursor:mass_tolerance", precursor_mass_tolerance);
    //protein_ids[0].setMetaValue("precursor:mass_tolerance_unit", precursor_mass_tolerance_unit_ppm ? "ppm" : "Da");

    search_params.setMetaValue("fragment:mass_tolerance_xlinks", fragment_mass_tolerance_xlinks);
    //protein_ids[0].setMetaValue("fragment:mass_tolerance", fragment_mass_tolerance);
    //protein_ids[0].setMetaValue("fragment:mass_tolerance_unit", fragment_mass_tolerance_unit_ppm ? "ppm" : "Da");

    search_params.setMetaValue("peptide:min_size", peptide_min_size);
    //protein_ids[0].setMetaValue("peptide:missed_cleavages", missed_cleavages);
    //protein_ids[0].setMetaValue("peptide:enzyme", enzyme_name);

    search_params.setMetaValue("cross_link:residue1", cross_link_residue1);
    search_params.setMetaValue("cross_link:residue2", cross_link_residue2);
    search_params.setMetaValue("cross_link:mass", cross_link_mass);
    //search_params.setMetaValue("cross_link:mass_isoshift", cross_link_mass_iso_shift);
    search_params.setMetaValue("cross_link:mass_monolink", cross_link_mass_mono_link);

    //protein_ids[0].setMetaValue("modifications:fixed", fixedModNames);
    //protein_ids[0].setMetaValue("modifications:variable", varModNames);
    search_params.setMetaValue("modifications:variable_max_per_peptide", max_variable_mods_per_peptide);

    search_params.setMetaValue("algorithm:candidate_search", ion_index_mode ? "ion-tag" : "enumeration");

    protein_ids[0].setSearchParameters(search_params);

    vector<PeptideIdentification> peptide_ids;


    progresslogger.startProgress(0, 1, "Digesting peptides...");

//#ifdef _OPENMP
//#pragma omp parallel for
//#endif
    // digest and filter database
    for (SignedSize fasta_index = 0; fasta_index < static_cast<SignedSize>(fasta_db.size()); ++fasta_index)
    {

      IF_MASTERTHREAD
      {
        progresslogger.setProgress(static_cast<SignedSize>(fasta_index) * NUMBER_OF_THREADS);
      }

      // store vector of substrings pointing in fasta database (bounded by pairs of begin, end iterators)
      vector<StringView> current_digest;
      digestor.digestUnmodifiedString(fasta_db[fasta_index].sequence, current_digest, min_peptide_length);

      for (vector<StringView>::iterator cit = current_digest.begin(); cit != current_digest.end(); ++cit)
      {
        // skip peptides with invalid AAs // TODO is this necessary? How can we adress such AAs?
        if (cit->getString().has('B') || cit->getString().has('O') || cit->getString().has('U') || cit->getString().has('X') || cit->getString().has('Z')) continue;

        // skip if no cross-linked residue
        bool skip = true;
        for (Size k = 0; k < cross_link_residue1.size(); k++)
        {
          if (cit->getString().find(cross_link_residue1[k]) < cit->getString().size()-1) skip = false;
        }
        for (Size k = 0; k < cross_link_residue2.size(); k++)
        {
            if (cit->getString().find(cross_link_residue2[k]) < cit->getString().size()-1) skip = false;
        }
        if (skip) continue;

        bool already_processed = false;
//#ifdef _OPENMP
//#pragma omp critical (processed_peptides_access)
//#endif
        {
          if (processed_peptides.find(*cit) != processed_peptides.end())
          {
            // peptide (and all modified variants) already processed so skip it
            already_processed = true;
          }
        }

        if (already_processed)
        {
          continue;
        }

        OpenProXLUtils::PeptidePosition position = OpenProXLUtils::INTERNAL;
        if (fasta_db[fasta_index].sequence.hasPrefix(cit->getString()))
        {
          position = OpenProXLUtils::C_TERM;
        } else if (fasta_db[fasta_index].sequence.hasSuffix(cit->getString()))
        {
          position = OpenProXLUtils::N_TERM;
        }

        vector<AASequence> all_modified_peptides;

        // generate all modified variants of a peptide
        // Note: no critial section is needed despite ResidueDB not beeing thread sage.
        //       It is only written to on introduction of novel modified residues. These residues have been already added above (single thread context).
        {
          AASequence aas = AASequence::fromString(cit->getString());
          ModifiedPeptideGenerator::applyFixedModifications(fixed_modifications.begin(), fixed_modifications.end(), aas);
          ModifiedPeptideGenerator::applyVariableModifications(variable_modifications.begin(), variable_modifications.end(), aas, max_variable_mods_per_peptide, all_modified_peptides);
        }

        for (SignedSize mod_pep_idx = 0; mod_pep_idx < static_cast<SignedSize>(all_modified_peptides.size()); ++mod_pep_idx)
        {
          const AASequence& candidate = all_modified_peptides[mod_pep_idx];
          OpenProXLUtils::PeptideMass pep_mass;
          pep_mass.peptide_mass = candidate.getMonoWeight();
          pep_mass.peptide_seq = candidate;
          pep_mass.position = position;

//#ifdef _OPENMP
//#pragma omp critical (processed_peptides_access)
//#endif
          {
            processed_peptides.insert(pair<StringView, AASequence>(*cit, candidate));
            peptide_masses.push_back(pep_mass);
          }
        }
      }
    }
    processed_peptides.clear();
    progresslogger.endProgress();

    // create spectrum generator
    TheoreticalSpectrumGenerator spectrum_generator;

    TheoreticalSpectrumGeneratorXLinks specGen_old;
    TheoreticalSpectrumGeneratorXLMS specGen;

    // Setting parameters for cross-link fragmentation
    Param specGenParams = specGen.getParameters();
    specGenParams.setValue("add_isotopes", "true", "If set to 1 isotope peaks of the product ion peaks are added");
    specGenParams.setValue("max_isotope", 2, "Defines the maximal isotopic peak which is added, add_isotopes must be set to 1");
    specGenParams.setValue("add_losses", "false", "Adds common losses to those ion expect to have them, only water and ammonia loss is considered");
    specGenParams.setValue("add_precursor_peaks", "false", "Adds peaks of the precursor to the spectrum, which happen to occur sometimes");
    specGenParams.setValue("add_abundant_immonium_ions", "false", "Add most abundant immonium ions");
    specGenParams.setValue("add_first_prefix_ion", "true", "If set to true e.g. b1 ions are added");
    specGenParams.setValue("add_y_ions", "true", "Add peaks of y-ions to the spectrum");
    specGenParams.setValue("add_b_ions", "true", "Add peaks of b-ions to the spectrum");
    specGenParams.setValue("add_a_ions", "false", "Add peaks of a-ions to the spectrum");
    specGenParams.setValue("add_c_ions", "false", "Add peaks of c-ions to the spectrum");
    specGenParams.setValue("add_x_ions", "false", "Add peaks of  x-ions to the spectrum");
    specGenParams.setValue("add_z_ions", "false", "Add peaks of z-ions to the spectrum");
    // TODO does nothing yet
    specGenParams.setValue("multiple_fragmentation_mode" , "false", "If set to true, multiple fragmentation events on the same cross-linked peptide pair are considered (HCD fragmentation)");
    specGen.setParameters(specGenParams);

    // TODO constant binsize for HashGrid computation
    double tolerance_binsize = 0.2;

    LOG_DEBUG << "Peptide candidates: " << processed_peptides.size() << endl;
    search_params = protein_ids[0].getSearchParameters();
    search_params.setMetaValue("MS:1001029", processed_peptides.size()); // number of sequences searched = MS:1001029
    protein_ids[0].setSearchParameters(search_params);

    //TODO refactor, so that only the used mode is initialized and the pre-scoring code only appears once
    // Initialize enumeration mode
    //multimap<double, pair<const AASequence*, const AASequence*> > enumerated_cross_link_masses;
    vector<OpenProXLUtils::XLPrecursor> enumerated_cross_link_masses;

    //TODO remove, adapt to ppm
//    HashGrid1D hg(0.0, 20000.0, tolerance_binsize);

//    map<AASequence*, MSSpectrum<RichPeak1D> > peptide_spectra;

    // Collect precursor MZs for filtering enumerated peptide pairs
    vector< double > spectrum_precursors;
    for (Size i = 0; i < spectra.size(); i++)
    {
      double current_precursor_mz = spectra[i].getPrecursors()[0].getMZ();
      double current_precursor_charge = spectra[i].getPrecursors()[0].getCharge();
      double current_precursor_mass = (current_precursor_mz * current_precursor_charge) - (current_precursor_charge * Constants::PROTON_MASS_U);
      spectrum_precursors.push_back(current_precursor_mass);
    }
    sort(spectrum_precursors.begin(), spectrum_precursors.end());
    cout << "Number of precursor masses in the spectra: " << spectrum_precursors.size() << endl;

    sort(peptide_masses.begin(), peptide_masses.end());
    // The largest peptides given a fixed maximal precursor mass are possible with loop links
    // Filter peptides using maximal loop link mass first
    double max_precursor_mass = spectrum_precursors[spectrum_precursors.size()-1];

    // compute absolute tolerance from relative, if necessary
    double allowed_error = 0;
    if (precursor_mass_tolerance_unit_ppm) // ppm
    {
      allowed_error = max_precursor_mass * precursor_mass_tolerance * 1e-6;
    }
    else // Dalton
    {
      allowed_error = precursor_mass_tolerance;
    }

    double max_peptide_mass = max_precursor_mass - cross_link_mass + allowed_error;

    cout << "Filtering peptides with precursors" << endl;

    // TODO peptides are sorted!
    // search for the closest one, cut off everything lager

    // iterating over the vector, while changing its size
    vector<OpenProXLUtils::PeptideMass>::iterator a = peptide_masses.begin();
    while(true)
    {
      if ( a->peptide_mass > max_peptide_mass )
      {
        peptide_masses.erase(a);
      }
      else
      {
        ++a;
      }
      if (a == peptide_masses.end())
      {
        break;
      }
    }



    if (!ion_index_mode)
    {
      progresslogger.startProgress(0, 1, "Enumerating cross-links...");
      enumerated_cross_link_masses = OpenProXLUtils::enumerateCrossLinksAndMasses_(peptide_masses, cross_link_mass, cross_link_mass_mono_link, cross_link_residue1, cross_link_residue2,
                                                                                                                                                    spectrum_precursors, precursor_mass_tolerance, precursor_mass_tolerance_unit_ppm);
      progresslogger.endProgress();
      cout << "Enumerated cross-links: " << enumerated_cross_link_masses.size() << endl;
      sort(enumerated_cross_link_masses.begin(), enumerated_cross_link_masses.end());
      cout << "Sorting of enumerated precursors finished" << endl;
    }
    else
    {
      // TODO refactor to remove ion_index mode completely
    }

    // for PScore, precompute ranks
    vector<vector<Size> > rankMap = PScore::calculateRankMap(spectra);

    // TODO test variables, can be removed, or set to be used in debug mode?
    double pScoreMax = 0;
    double TICMax = 0;
    double wTICMax = 0;
    double intsumMax = 0;
    double matchOddsMax = 0;
    double xcorrxMax = 0;
    double xcorrcMax = 0;
    double maxMatchCount = 0;
    double sumMatchCount = 0;

    // iterate over all spectra
    progresslogger.startProgress(0, 1, "Matching to theoretical spectra and scoring...");
    vector< vector< CrossLinkSpectrumMatch > > all_top_csms;

    Size spectrum_counter = 0;

    cout << "Spectra left after preprocessing and filtering: " << spectra.size() << " of " << unprocessed_spectra.size() << endl;

// Multithreading options: schedule: static, dynamic, guided with chunk size
#ifdef _OPENMP
#pragma omp parallel for schedule(guided)
#endif
    for (SignedSize scan_index = 0; scan_index < static_cast<SignedSize>(spectra.size()); ++scan_index)
    {

      const PeakSpectrum& spectrum = spectra[scan_index];

      // TODO probably not necessary, if this is appropriately done in preprocessing
      if ( spectrum.size() > peptide_min_size )
      {
        LOG_DEBUG << "################## NEW SPECTRUM ##############################" << endl;
        //LOG_DEBUG << "Scan index: " << scan_index << "\tSpectrum native ID: " << spectrum.getNativeID()  << endl;

#ifdef _OPENMP
#pragma omp critical
#endif
        {
          spectrum_counter++;
          cout << "Processing spectrum " << spectrum_counter << " / " << spectra.size() << "\tSpectrum ID: " << spectrum.getNativeID()  << endl;
        }

        const double precursor_charge = spectrum.getPrecursors()[0].getCharge();
        const double precursor_mz = spectrum.getPrecursors()[0].getMZ();
        const double precursor_mass = (precursor_mz * static_cast<double>(precursor_charge)) - (static_cast<double>(precursor_charge) * Constants::PROTON_MASS_U);

        // Mean ion intensity (light spectrum, TODO add heavy spectrum?)
        double mean_intensity= 0;
        for (SignedSize j = 0; j < static_cast<SignedSize>(spectrum.size()); ++j)
        {
          mean_intensity += spectrum[j].getIntensity();
        }
        mean_intensity = mean_intensity / spectrum.size();

        const PeakSpectrum& common_peaks = spectra[scan_index];

        vector< CrossLinkSpectrumMatch > top_csms_spectrum;


        // determine candidates
        //vector<pair<const AASequence*, const AASequence*> > candidates;
        vector< OpenProXLUtils::XLPrecursor > candidates;
        double allowed_error = 0;

        if (ion_index_mode)
        {
          // TODO refactor to remove ion_index mode completely
        } else // enumeration mode
        {
          //LOG_DEBUG << "Number of common peaks, xlink peaks: " << preprocessed_pair_spectra.spectra_common_peaks[scan_index].size() << "\t" << preprocessed_pair_spectra.spectra_xlink_peaks[scan_index].size();

          // determine MS2 precursors that match to the current peptide mass
          vector< OpenProXLUtils::XLPrecursor >::const_iterator low_it;
          vector< OpenProXLUtils::XLPrecursor >::const_iterator up_it;

          if (precursor_mass_tolerance_unit_ppm) // ppm
          {
            allowed_error = precursor_mass * precursor_mass_tolerance * 1e-6;
          }
          else // Dalton
          {
            allowed_error = precursor_mass_tolerance;
          }
#ifdef _OPENMP
#pragma omp critical (enumerated_cross_link_masses_access)
#endif
          {
            low_it = lower_bound(enumerated_cross_link_masses.begin(), enumerated_cross_link_masses.end(), precursor_mass - allowed_error);
            up_it =  upper_bound(enumerated_cross_link_masses.begin(), enumerated_cross_link_masses.end(), precursor_mass + allowed_error);
          }

          if (low_it != up_it) // no matching precursor in data
          {
            for (; low_it != up_it; ++low_it)
            {
              candidates.push_back(*low_it);
            }
          }
        }

        // Find all positions of lysine (K) in the peptides (possible scross-linking sites), create cross_link_candidates with all combinations
        vector <ProteinProteinCrossLink> cross_link_candidates;
        for (Size i = 0; i != candidates.size(); ++i)
        {
          //pair<const AASequence*, const AASequence*> candidate = candidates[i];
          OpenProXLUtils::XLPrecursor candidate = candidates[i];
          vector <SignedSize> link_pos_first;
          vector <SignedSize> link_pos_second;
          AASequence peptide_first = peptide_masses[candidate.alpha_index].peptide_seq;
          AASequence peptide_second;
          if (candidate.beta_index)
          {
            peptide_second = peptide_masses[candidate.beta_index].peptide_seq;
          }
          String seq_first = peptide_first.toUnmodifiedString();
          String seq_second =  peptide_second.toUnmodifiedString();


          // TODO mono-links and loop-links with different masses can be generated for the same precursor mass, but only one of them can be valid each time.
          // Find out which is the case. But it should not happen often enough to slow down the tool significantly.
          bool is_loop = abs(precursor_mass - (peptide_first.getMonoWeight() + cross_link_mass)) <= allowed_error;

          for (Size k = 0; k < seq_first.size()-1; ++k)
          {
            for (Size x = 0; x < cross_link_residue1.size(); ++x)
            {
              if (seq_first.substr(k, 1) == cross_link_residue1[x]) link_pos_first.push_back(k);
            }
          }
          if (candidate.beta_index)
          {
            for (Size k = 0; k < seq_second.size()-1; ++k)
            {
              for (Size x = 0; x < cross_link_residue2.size(); ++x)
              {
                if (seq_second.substr(k, 1) == cross_link_residue2[x]) link_pos_second.push_back(k);
              }
            }
          } else
          {
            // Second position defining a mono-link and the second positions on the same peptide for loop links (only one of these two is valid for any specific precursor)
            if (!is_loop)
            {
              link_pos_second.push_back(-1);
            }
            else
            {
              for (Size k = 0; k < seq_first.size()-1; ++k)
              {
                for (Size x = 0; x < cross_link_residue2.size(); ++x)
                {
                  if (seq_first.substr(k, 1) == cross_link_residue2[x]) link_pos_second.push_back(k);
                }
              }
            }
          }

          // Determine larger peptide (alpha) by sequence length, use mass as tie breaker
          bool alpha_first = true;

          if (seq_second.size() > seq_first.size())
          {
            alpha_first = false;
          } else if (seq_second.size() == seq_first.size() && peptide_second.getMonoWeight() > peptide_first.getMonoWeight())
          {
            alpha_first = false;
          }

          // generate cross_links for all valid combinations
          for (Size x = 0; x < link_pos_first.size(); ++x)
          {
            for (Size y = 0; y < link_pos_second.size(); ++y)
            {
              ProteinProteinCrossLink cross_link_candidate;
              if (seq_second.size() == 0)
              {
                // if loop link, and the positions are the same, then it is linking the same residue with itself,  skip this combination, also pos1 > pos2 would be the same link as pos1 < pos2
                if ( (link_pos_first[x] >= link_pos_second[y]) && (link_pos_second[y] != -1) )
                {
                  continue;
                }
              }
              if (alpha_first)
              {
                cross_link_candidate.alpha = peptide_first;
                cross_link_candidate.beta = peptide_second;
                cross_link_candidate.cross_link_position.first = link_pos_first[x];
                cross_link_candidate.cross_link_position.second = link_pos_second[y];
              }
              else
              {
                cross_link_candidate.alpha = peptide_second;
                cross_link_candidate.beta = peptide_first;
                cross_link_candidate.cross_link_position.first = link_pos_second[y];
                cross_link_candidate.cross_link_position.second = link_pos_first[x];
              }
              // Cross-linker mass is only one of the mono-link masses, if there is no second position (second == -1), otherwise the normal linker mass
              if (link_pos_second[y] != -1)
              {
                cross_link_candidate.cross_linker_mass = cross_link_mass;
                cross_link_candidate.cross_linker_name = cross_link_name;
                cross_link_candidates.push_back(cross_link_candidate);
              }
              else
              {
                for (Size k = 0; k < cross_link_mass_mono_link.size(); ++k)
                {
                  // only use the correct mono-links (at this point we know it is a mono-link, but not which one, so loop over all and compare precursors)
                  if (abs(precursor_mass - (peptide_first.getMonoWeight() + cross_link_mass_mono_link[k])) <= allowed_error)
                  {
                    cross_link_candidate.cross_linker_mass = cross_link_mass_mono_link[k];
                    cross_link_candidate.cross_linker_name = cross_link_name;
                    cross_link_candidates.push_back(cross_link_candidate);
                  }
                }
              }
            }
          }
        }

#ifdef _OPENMP
#pragma omp critical
#endif
        cout << "Number of candidates for this spectrum: " << candidates.size() << endl;

        // lists for one spectrum, to determine best match to the spectrum
        vector< CrossLinkSpectrumMatch > all_csms_spectrum;

        // TODO variables for benchmarking and testing purposes
#ifdef _OPENMP
#pragma omp critical (max_subscore_variable_access)
#endif
        {
          if (cross_link_candidates.size() > maxMatchCount) maxMatchCount = cross_link_candidates.size();
          sumMatchCount += cross_link_candidates.size();
        }

        for (Size i = 0; i != cross_link_candidates.size(); ++i)
        {
          ProteinProteinCrossLink cross_link_candidate = cross_link_candidates[i];
          double candidate_mz = (cross_link_candidate.alpha.getMonoWeight() + cross_link_candidate.beta.getMonoWeight() +  cross_link_candidate.cross_linker_mass+ (static_cast<double>(precursor_charge) * Constants::PROTON_MASS_U)) / precursor_charge;

          LOG_DEBUG << "Pair: " << cross_link_candidate.alpha.toString() << "-" << cross_link_candidate.beta.toString() << " matched to light spectrum " << scan_index << "\t and heavy spectrum " << scan_index
              << " with m/z: " << precursor_mz << "\t" << "and candidate m/z: " << candidate_mz << "\tK Positions: " << cross_link_candidate.cross_link_position.first << "\t" << cross_link_candidate.cross_link_position.second << endl;
//          LOG_DEBUG << a->second.getMonoWeight() << ", " << b->second.getMonoWeight() << " cross_link_mass: " <<  cross_link_mass <<  endl;

	  CrossLinkSpectrumMatch csm;
	  csm.cross_link = cross_link_candidate;

	  PeakSpectrum theoretical_spec_common_alpha;
	  PeakSpectrum theoretical_spec_common_beta;
	  PeakSpectrum theoretical_spec_xlinks_alpha;
	  PeakSpectrum theoretical_spec_xlinks_beta;

	  bool type_is_cross_link = cross_link_candidate.getType() == ProteinProteinCrossLink::CROSS;
	  bool type_is_loop = cross_link_candidate.getType() == ProteinProteinCrossLink::LOOP;
	  Size link_pos_B = 0;
	  if (type_is_loop)
	  {
	    link_pos_B = cross_link_candidate.cross_link_position.second;
	  }

//          specGen.getCommonIonSpectrum(theoretical_spec_common_alpha, cross_link_candidate, 2, true);
          specGen.getCommonIonSpectrum(theoretical_spec_common_alpha, cross_link_candidate.alpha, cross_link_candidate.cross_link_position.first, true, 2, link_pos_B);
          if (type_is_cross_link)
          {
//            specGen.getCommonIonSpectrum(theoretical_spec_common_beta, cross_link_candidate, 2, false);
//            specGen.getXLinkIonSpectrum(theoretical_spec_xlinks_alpha, theoretical_spec_xlinks_beta, cross_link_candidate, 2, precursor_charge);
            specGen.getCommonIonSpectrum(theoretical_spec_common_beta, cross_link_candidate.beta, cross_link_candidate.cross_link_position.second, false, 2);
            specGen.getXLinkIonSpectrum(theoretical_spec_xlinks_alpha, cross_link_candidate.alpha, cross_link_candidate.cross_link_position.first, precursor_mass, true, 1, precursor_charge);
            specGen.getXLinkIonSpectrum(theoretical_spec_xlinks_beta, cross_link_candidate.beta, cross_link_candidate.cross_link_position.second, precursor_mass, false, 1, precursor_charge);
          } else
          {
            // Function for mono-links or loop-links
//            specGen.getXLinkIonSpectrum(theoretical_spec_xlinks_alpha, cross_link_candidate, 2, precursor_charge);
            specGen.getXLinkIonSpectrum(theoretical_spec_xlinks_alpha, cross_link_candidate.alpha, cross_link_candidate.cross_link_position.first, precursor_mass, true, 2, precursor_charge, link_pos_B);
          }

          // Something like this can happen, e.g. with a loop link connecting the first and last residue of a peptide
          if ( (theoretical_spec_common_alpha.size() < 1) || (theoretical_spec_xlinks_alpha.size() < 1) )
          {
            continue;
          }

          vector< pair< Size, Size > > matched_spec_common_alpha;
          vector< pair< Size, Size > > matched_spec_common_beta;
          vector< pair< Size, Size > > matched_spec_xlinks_alpha;
          vector< pair< Size, Size > > matched_spec_xlinks_beta;

          getSpectrumAlignment(matched_spec_common_alpha, theoretical_spec_common_alpha, spectrum, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm);
          getSpectrumAlignment(matched_spec_common_beta, theoretical_spec_common_beta, spectrum, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm);
          getSpectrumAlignment(matched_spec_xlinks_alpha, theoretical_spec_xlinks_alpha, spectrum, fragment_mass_tolerance_xlinks, fragment_mass_tolerance_unit_ppm);
          getSpectrumAlignment(matched_spec_xlinks_beta, theoretical_spec_xlinks_beta, spectrum, fragment_mass_tolerance_xlinks, fragment_mass_tolerance_unit_ppm);

          LOG_DEBUG << "Spectrum sizes: " << spectrum.size() << " || " << theoretical_spec_common_alpha.size() <<  " | " << theoretical_spec_common_beta.size()
                                <<  " | " << theoretical_spec_xlinks_alpha.size() <<  " | " << theoretical_spec_xlinks_beta.size() << endl;
          LOG_DEBUG << "Matched peaks: " << matched_spec_common_alpha.size() << " | " << matched_spec_common_beta.size()
                                <<  " | " << matched_spec_xlinks_alpha.size() <<  " | " << matched_spec_xlinks_beta.size() << endl;

          // Pre-Score calculations
          Size matched_alpha_count = matched_spec_common_alpha.size() + matched_spec_xlinks_alpha.size();
          Size theor_alpha_count = theoretical_spec_common_alpha.size() + theoretical_spec_xlinks_alpha.size();
          Size matched_beta_count = matched_spec_common_beta.size() + matched_spec_xlinks_beta.size();
          Size theor_beta_count = theoretical_spec_common_beta.size() + theoretical_spec_xlinks_beta.size();

          if (matched_alpha_count + matched_beta_count > 0)
          {
            // Simplified pre-Score
            //float pre_score = preScore(matched_fragments_theor_spec.size(), theoretical_spec.size());
            double pre_score = 0;
            if (type_is_cross_link)
            {
               pre_score = OpenProXLUtils::preScore(matched_alpha_count, theor_alpha_count, matched_beta_count, theor_beta_count);
             }
             else
             {
              pre_score = OpenProXLUtils::preScore(matched_alpha_count, theor_alpha_count);
             }
             //LOG_DEBUG << "Number of matched peaks to theor. spectrum: " << matched_alpha_count << "\t" << matched_beta_count << endl;
             //LOG_DEBUG << "Number of theoretical ions: " << theor_alpha_count << "\t" << theor_beta_count << endl;
             //LOG_DEBUG << "Pre Score: " << pre_score << endl;
             //LOG_DEBUG << "Peptide size: " << a->second.size() << "\t" << b->second.size() << "\t" << "K Pos:" << K_pos_a[i] << "\t" << K_pos_b[i] << endl;

#ifdef _OPENMP
#pragma omp critical (max_subscore_variable_access)
#endif
             if (pre_score > pScoreMax) pScoreMax = pre_score;

             // compute intsum score
             double intsum = OpenProXLUtils::total_matched_current(matched_spec_common_alpha, matched_spec_common_beta, matched_spec_xlinks_alpha, matched_spec_xlinks_beta, spectrum, spectrum);


              // Total ion intensity of light spectrum
              // sum over common and xlink ion spectra instead of unfiltered
              double total_current = 0;
              for (SignedSize j = 0; j < static_cast<SignedSize>(spectrum.size()); ++j)
              {
                total_current += spectrum[j].getIntensity();
              }
              double TIC = intsum / total_current;

#ifdef _OPENMP
#pragma omp critical (max_subscore_variable_access)
#endif
              if (TIC > TICMax) TICMax = TIC;

              // TIC_alpha and _beta
              double intsum_alpha = OpenProXLUtils::matched_current_chain(matched_spec_common_alpha, matched_spec_xlinks_alpha, spectrum, spectrum);
              double intsum_beta = 0;
              if (type_is_cross_link)
              {
                intsum_beta = OpenProXLUtils::matched_current_chain(matched_spec_common_beta, matched_spec_xlinks_beta, spectrum, spectrum);
              }

              // normalize TIC_alpha and  _beta
              if ((intsum_alpha + intsum_beta) > 0.0)
              {
                intsum_alpha = intsum_alpha * intsum / (intsum_alpha + intsum_beta);
                intsum_beta = intsum_beta *  intsum / (intsum_alpha + intsum_beta);
              }

              // compute wTIC
              double wTIC = OpenProXLUtils::weighted_TIC_score(cross_link_candidate.alpha.size(), cross_link_candidate.beta.size(), intsum_alpha, intsum_beta, intsum, total_current, type_is_cross_link);

#ifdef _OPENMP
#pragma omp critical (max_subscore_variable_access)
#endif
              {
                if (wTIC > wTICMax) wTICMax = wTIC;
                if (intsum > intsumMax) intsumMax = intsum;
              }

              // maximal xlink ion charge = (Precursor charge - 1), minimal xlink ion charge: 2
              Size n_xlink_charges = (precursor_charge - 1) - 2;
              if (n_xlink_charges < 1) n_xlink_charges = 1;

              // compute match odds (unweighted), the 3 is the number of charge states in the theoretical spectra
              double match_odds_c_alpha = OpenProXLUtils::match_odds_score(theoretical_spec_common_alpha, matched_spec_common_alpha, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, false);
              double match_odds_x_alpha = OpenProXLUtils::match_odds_score(theoretical_spec_xlinks_alpha, matched_spec_xlinks_alpha, fragment_mass_tolerance_xlinks , fragment_mass_tolerance_unit_ppm, true, n_xlink_charges);
              double match_odds = 0;
              //cout << "TEST Match_odds_c_alpha: " << match_odds_c_alpha << "\t x_alpha: " << match_odds_x_alpha << endl;
              if (type_is_cross_link)
              {
                double match_odds_c_beta = OpenProXLUtils::match_odds_score(theoretical_spec_common_beta, matched_spec_common_beta, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, false);
                double match_odds_x_beta = OpenProXLUtils::match_odds_score(theoretical_spec_xlinks_beta, matched_spec_xlinks_beta, fragment_mass_tolerance_xlinks, fragment_mass_tolerance_unit_ppm, true, n_xlink_charges);
                match_odds = (match_odds_c_alpha + match_odds_x_alpha + match_odds_c_beta + match_odds_x_beta) / 4;
                //cout << "TEST Match_odds_c_beta: " << match_odds_c_beta << "\t x_beta: " << match_odds_x_beta<< endl;
              }
              else
              {
                match_odds = (match_odds_c_alpha + match_odds_x_alpha) / 2;
              }

#ifdef _OPENMP
#pragma omp critical (max_subscore_variable_access)
#endif
              if (match_odds > matchOddsMax) matchOddsMax = match_odds;

              //Cross-correlation
              PeakSpectrum theoretical_spec_common;
              PeakSpectrum theoretical_spec_xlinks;

              if (type_is_cross_link)
              {
//                theoretical_spec_common.reserve(theoretical_spec_common_alpha.size() + theoretical_spec_common_beta.size());
//                theoretical_spec_xlinks.reserve( theoretical_spec_xlinks_alpha.size() + theoretical_spec_xlinks_beta.size());
//                theoretical_spec_common.insert(theoretical_spec_common.end(), theoretical_spec_common_alpha.begin(), theoretical_spec_common_alpha.end());
//                theoretical_spec_common.insert(theoretical_spec_common.end(), theoretical_spec_common_beta.begin(), theoretical_spec_common_beta.end());
//                theoretical_spec_xlinks.insert(theoretical_spec_xlinks.end(), theoretical_spec_xlinks_alpha.begin(), theoretical_spec_xlinks_alpha.end());
//                theoretical_spec_xlinks.insert(theoretical_spec_xlinks.end(), theoretical_spec_xlinks_beta.begin(), theoretical_spec_xlinks_beta.end());
//                theoretical_spec_common.sortByPosition();
//                theoretical_spec_xlinks.sortByPosition();

                theoretical_spec_common = OpenProXLUtils::mergeAnnotatedSpectra(theoretical_spec_common_alpha, theoretical_spec_common_beta);
                theoretical_spec_xlinks = OpenProXLUtils::mergeAnnotatedSpectra(theoretical_spec_xlinks_alpha, theoretical_spec_xlinks_beta);
              }
              else
              {
                theoretical_spec_common = theoretical_spec_common_alpha;
                theoretical_spec_xlinks = theoretical_spec_xlinks_alpha;
              }

              PeakSpectrum theoretical_spec = OpenProXLUtils::mergeAnnotatedSpectra(theoretical_spec_common, theoretical_spec_xlinks);
//              theoretical_spec.reserve(theoretical_spec_common.size() + theoretical_spec_xlinks.size());
//              theoretical_spec.insert(theoretical_spec.end(), theoretical_spec_common.begin(), theoretical_spec_common.end());
//              theoretical_spec.insert(theoretical_spec.end(), theoretical_spec_xlinks.begin(), theoretical_spec_xlinks.end());
//              theoretical_spec.sortByPosition();

              PeakSpectrum theoretical_spec_alpha = OpenProXLUtils::mergeAnnotatedSpectra(theoretical_spec_common_alpha, theoretical_spec_xlinks_alpha);
//              theoretical_spec_alpha.reserve(theoretical_spec_common_alpha.size() + theoretical_spec_xlinks_alpha.size());
//              theoretical_spec_alpha.insert(theoretical_spec_alpha.end(), theoretical_spec_common_alpha.begin(), theoretical_spec_common_alpha.end());
//              theoretical_spec_alpha.insert(theoretical_spec_alpha.end(), theoretical_spec_xlinks_alpha.begin(), theoretical_spec_xlinks_alpha.end());

              PeakSpectrum theoretical_spec_beta;
              if (type_is_cross_link)
              {
                theoretical_spec_beta = OpenProXLUtils::mergeAnnotatedSpectra(theoretical_spec_common_beta, theoretical_spec_xlinks_beta);
              }
//              theoretical_spec_beta.reserve(theoretical_spec_common_beta.size() + theoretical_spec_xlinks_beta.size());
//              theoretical_spec_beta.insert(theoretical_spec_beta.end(), theoretical_spec_common_beta.begin(), theoretical_spec_common_beta.end());
//              theoretical_spec_beta.insert(theoretical_spec_beta.end(), theoretical_spec_xlinks_beta.begin(), theoretical_spec_xlinks_beta.end());

              vector< double > xcorrx = OpenProXLUtils::xCorrelation(spectrum, theoretical_spec_xlinks, 5, 0.3);
              vector< double > xcorrc = OpenProXLUtils::xCorrelation(spectrum, theoretical_spec_common, 5, 0.2);

                // TODO save time: only needs to be done once per light spectrum, here it is repeated for cross-link each candidate
              vector< double > aucorrx = OpenProXLUtils::xCorrelation(spectrum, spectrum, 5, 0.3);
              vector< double > aucorrc = OpenProXLUtils::xCorrelation(spectrum, spectrum, 5, 0.2);
              double aucorr_sumx = accumulate(aucorrx.begin(), aucorrx.end(), 0.0);
              double aucorr_sumc = accumulate(aucorrc.begin(), aucorrc.end(), 0.0);
              double xcorrx_max = accumulate(xcorrx.begin(), xcorrx.end(), 0.0) / aucorr_sumx;
              double xcorrc_max = accumulate(xcorrc.begin(), xcorrc.end(), 0.0) / aucorr_sumc;
//                LOG_DEBUG << "xCorrelation X: " << xcorrx << endl;
//                LOG_DEBUG << "xCorrelation C: " << xcorrc << endl;
//                LOG_DEBUG << "Autocorr: " << aucorr << "\t Autocorr_sum: " << aucorr_sum << "\t xcorrx_max: " << xcorrx_max << "\t xcorrc_max: " << xcorrc_max << endl;

//############################# TESTING SCORES ##############################################

                map<Size, PeakSpectrum> peak_level_spectra = PScore::calculatePeakLevelSpectra(spectrum, rankMap[scan_index]);
                csm.PScoreCommon = PScore::computePScore(fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, peak_level_spectra, theoretical_spec_common);
                csm.PScoreXlink = PScore::computePScore(fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, peak_level_spectra, theoretical_spec_xlinks);
                csm.PScoreBoth = PScore::computePScore(fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, peak_level_spectra, theoretical_spec);
                csm.PScoreAlpha = PScore::computePScore(fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, peak_level_spectra, theoretical_spec_alpha);
                csm.PScoreBeta = PScore::computePScore(fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, peak_level_spectra, theoretical_spec_beta);

                csm.HyperCommon = HyperScore::compute(fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, spectrum, theoretical_spec_common);
                csm.HyperAlpha = HyperScore::compute(fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, spectrum, theoretical_spec_alpha);
                csm.HyperBeta = HyperScore::compute(fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, spectrum, theoretical_spec_beta);
                csm.HyperXlink = HyperScore::compute(fragment_mass_tolerance_xlinks, fragment_mass_tolerance_unit_ppm, spectrum, theoretical_spec_xlinks);
                csm.HyperBoth = HyperScore::compute(fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, spectrum, theoretical_spec);




//############################# END TESTING SCORES ###########################################



#ifdef _OPENMP
#pragma omp critical (max_subscore_variable_access)
#endif
                {
                  if (xcorrx_max > xcorrxMax) xcorrxMax = xcorrx_max;
                  if (xcorrc_max > xcorrcMax) xcorrcMax = xcorrc_max;
                }

                // Compute score from the 4 scores and 4 weights
                double xcorrx_weight = 2.488;
                double xcorrc_weight = 21.279;
                //double match_odds_weight = 1.973;
                // TODO match-.odds score does not work for high res data, very low tolerances lead to low probabilities for random matches and the maximal score is reached very often
                double match_odds_weight = 0.1;
                double wTIC_weight = 12.829;
                double intsum_weight = 1.8;

                double score = xcorrx_weight * xcorrx_max + xcorrc_weight * xcorrc_max + match_odds_weight * match_odds + wTIC_weight * wTIC + intsum_weight * intsum;

                csm.score = score;
                csm.pre_score = pre_score;
                csm.percTIC = TIC;
                csm.wTIC = wTIC;
                csm.int_sum = intsum;
                csm.match_odds = match_odds;
//                csm.xcorrx = xcorrx;
                csm.xcorrx_max = xcorrx_max;
//                csm.xcorrc = xcorrc;
                csm.xcorrc_max = xcorrc_max;
                csm.matched_common_alpha = matched_spec_common_alpha.size();
                csm.matched_common_beta = matched_spec_common_beta.size();
                csm.matched_xlink_alpha = matched_spec_xlinks_alpha.size();
                csm.matched_xlink_beta = matched_spec_xlinks_beta.size();
                csm.scan_index_light = scan_index;
                csm.scan_index_heavy = -1;


                // write fragment annotations
                LOG_DEBUG << "Start writing annotations" << endl;
                vector<PeptideHit::FragmentAnnotation> frag_annotations;

                PeakSpectrum::IntegerDataArray common_alpha_charges = theoretical_spec_common_alpha.getIntegerDataArrays()[0];
                PeakSpectrum::StringDataArray common_alpha_names = theoretical_spec_common_alpha.getStringDataArrays()[0];
                for (Size k = 0; k < matched_spec_common_alpha.size(); ++k)
                {
                  PeptideHit::FragmentAnnotation frag_anno;
//                  frag_anno.charge = static_cast<int>(theoretical_spec_common_alpha[matched_spec_common_alpha[k].first].getMetaValue("z"));
                  frag_anno.mz = spectrum[matched_spec_common_alpha[k].second].getMZ();
                  frag_anno.intensity = spectrum[matched_spec_common_alpha[k].second].getIntensity();
//                  frag_anno.annotation = theoretical_spec_common_alpha[matched_spec_common_alpha[k].first].getMetaValue("IonName");

                  frag_anno.charge = common_alpha_charges[matched_spec_common_alpha[k].first];
                  frag_anno.annotation = common_alpha_names[matched_spec_common_alpha[k].first];
                  frag_annotations.push_back(frag_anno);
                }

                if (theoretical_spec_common_beta.getIntegerDataArrays().size() > 0)
                {
                  PeakSpectrum::IntegerDataArray common_beta_charges = theoretical_spec_common_beta.getIntegerDataArrays()[0];
                  PeakSpectrum::StringDataArray common_beta_names = theoretical_spec_common_beta.getStringDataArrays()[0];
                  for (Size k = 0; k < matched_spec_common_beta.size(); ++k)
                  {
                    PeptideHit::FragmentAnnotation frag_anno;
//                    frag_anno.charge = static_cast<int>(theoretical_spec_common_beta[matched_spec_common_beta[k].first].getMetaValue("z"));
                    frag_anno.mz = spectrum[matched_spec_common_beta[k].second].getMZ();
                    frag_anno.intensity = spectrum[matched_spec_common_beta[k].second].getIntensity();
//                    frag_anno.annotation = theoretical_spec_common_beta[matched_spec_common_beta[k].first].getMetaValue("IonName");

                    frag_anno.charge = common_beta_charges[matched_spec_common_beta[k].first];
                    frag_anno.annotation = common_beta_names[matched_spec_common_beta[k].first];
                    frag_annotations.push_back(frag_anno);
                  }
                }

                PeakSpectrum::IntegerDataArray xlinks_alpha_charges = theoretical_spec_xlinks_alpha.getIntegerDataArrays()[0];
                PeakSpectrum::StringDataArray xlinks_alpha_names = theoretical_spec_xlinks_alpha.getStringDataArrays()[0];
                for (Size k = 0; k < matched_spec_xlinks_alpha.size(); ++k)
                {
                  PeptideHit::FragmentAnnotation frag_anno;
//                  frag_anno.charge = static_cast<int>(theoretical_spec_xlinks_alpha[matched_spec_xlinks_alpha[k].first].getMetaValue("z"));
                  frag_anno.mz = spectrum[matched_spec_xlinks_alpha[k].second].getMZ();
                  frag_anno.intensity = spectrum[matched_spec_xlinks_alpha[k].second].getIntensity();
//                  frag_anno.annotation = theoretical_spec_xlinks_alpha[matched_spec_xlinks_alpha[k].first].getMetaValue("IonName");

                  frag_anno.charge = xlinks_alpha_charges[matched_spec_xlinks_alpha[k].first];
                  frag_anno.annotation = xlinks_alpha_names[matched_spec_xlinks_alpha[k].first];
                  frag_annotations.push_back(frag_anno);
                }

                if (theoretical_spec_xlinks_beta.getIntegerDataArrays().size() > 0)
                {
                  PeakSpectrum::IntegerDataArray xlinks_beta_charges = theoretical_spec_xlinks_beta.getIntegerDataArrays()[0];
                  PeakSpectrum::StringDataArray xlinks_beta_names = theoretical_spec_xlinks_beta.getStringDataArrays()[0];
                  for (Size k = 0; k < matched_spec_xlinks_beta.size(); ++k)
                  {
                    PeptideHit::FragmentAnnotation frag_anno;
//                    frag_anno.charge = static_cast<int>(theoretical_spec_xlinks_beta[matched_spec_xlinks_beta[k].first].getMetaValue("z"));
                    frag_anno.mz = spectrum[matched_spec_xlinks_beta[k].second].getMZ();
                    frag_anno.intensity = spectrum[matched_spec_xlinks_beta[k].second].getIntensity();
//                    frag_anno.annotation = theoretical_spec_xlinks_beta[matched_spec_xlinks_beta[k].first].getMetaValue("IonName");

                    frag_anno.charge = xlinks_beta_charges[matched_spec_xlinks_beta[k].first];
                    frag_anno.annotation = xlinks_beta_names[matched_spec_xlinks_beta[k].first];
                    frag_annotations.push_back(frag_anno);
                  }
                }
                LOG_DEBUG << "End writing annotations, size: " << frag_annotations.size() << endl;

                // make annotations unique
                sort(frag_annotations.begin(), frag_annotations.end());
                vector<PeptideHit::FragmentAnnotation>::iterator last_unique_anno = unique(frag_annotations.begin(), frag_annotations.end());
                if (last_unique_anno != frag_annotations.end())
                {
                  LOG_DEBUG << "uniqiifying: " << endl;
                  for (vector<PeptideHit::FragmentAnnotation>::iterator double_frag = last_unique_anno; double_frag != frag_annotations.end(); ++double_frag)
                  {
                    LOG_DEBUG << "anno: " << double_frag->annotation << "\tcharge: " << double_frag->charge << "\tmz: " << double_frag->mz << "\tint: " << double_frag->intensity << endl;
                  }
                  frag_annotations.erase(last_unique_anno, frag_annotations.end());
                  LOG_DEBUG << "Fragment annotations were uniquified, new size: " << frag_annotations.size() << endl;
                }

                csm.frag_annotations = frag_annotations;

                all_csms_spectrum.push_back(csm);
              }
            } // candidates for peak finished, determine best matching candidate

            Int top = 0;

            // collect top n matches to spectrum
            while(!all_csms_spectrum.empty() && top < number_top_hits)
            {
              top++;

              //double max_score = *max_element(candidate_score.begin(), candidate_score.end());
              Int max_position = distance(all_csms_spectrum.begin(), max_element(all_csms_spectrum.begin(), all_csms_spectrum.end()));
              all_csms_spectrum[max_position].rank = top;
              top_csms_spectrum.push_back(all_csms_spectrum[max_position]);
              all_csms_spectrum.erase(all_csms_spectrum.begin() + max_position);

              LOG_DEBUG << "Score: " << all_csms_spectrum[max_position].score << "\t wTIC: " << all_csms_spectrum[max_position].wTIC << "\t xcorrx: " << all_csms_spectrum[max_position].xcorrx_max
                      << "\t xcorrc: " << all_csms_spectrum[max_position].xcorrc_max << "\t match-odds: " << all_csms_spectrum[max_position].match_odds << "\t Intsum: " << all_csms_spectrum[max_position].int_sum << endl;

              if (all_csms_spectrum[max_position].cross_link.getType() == ProteinProteinCrossLink::CROSS)
              {
                LOG_DEBUG << "Matched ions calpha , cbeta , xalpha , xbeta" << "\t" << all_csms_spectrum[max_position].matched_common_alpha << "\t" << all_csms_spectrum[max_position].matched_common_beta
                        << "\t" << all_csms_spectrum[max_position].matched_xlink_alpha <<  "\t" << all_csms_spectrum[max_position].matched_xlink_beta << endl;
              }
              else
              {
                LOG_DEBUG << "Matched ions common, cross-links " << all_csms_spectrum[max_position].matched_common_alpha << "\t" << all_csms_spectrum[max_position].matched_xlink_alpha << endl;
              }
            }

            Size all_top_csms_current_index = 0;
#ifdef _OPENMP
#pragma omp critical (all_top_csms_access)
#endif
            {
              all_top_csms.push_back(top_csms_spectrum);
              all_top_csms_current_index = all_top_csms.size()-1;
            }

            // Write PeptideIdentifications and PeptideHits for n top hits
            for (Size i = 0; i < top_csms_spectrum.size(); ++i)
            {
              PeptideIdentification peptide_id;

              String xltype = "cross-link";
              SignedSize alpha_pos = top_csms_spectrum[i].cross_link.cross_link_position.first;
              SignedSize beta_pos = top_csms_spectrum[i].cross_link.cross_link_position.second;

              if (top_csms_spectrum[i].cross_link.getType() == ProteinProteinCrossLink::MONO)
              {
                xltype = "mono-link";
              }
              else if (top_csms_spectrum[i].cross_link.getType() == ProteinProteinCrossLink::LOOP)
              {
                xltype = "loop-link";
              }

              PeptideHit ph_alpha, ph_beta;
              // Set monolink as a modification or add MetaValue for cross-link identity and mass
              AASequence seq_alpha = top_csms_spectrum[i].cross_link.alpha;
              if (top_csms_spectrum[i].cross_link.getType() == ProteinProteinCrossLink::MONO)
              {
                //AASequence seq_alpha = top_csms_spectrum[i].cross_link.alpha;
                vector< String > mods;
                const String residue = seq_alpha[alpha_pos].getOneLetterCode();
                ModificationsDB::getInstance()->getModificationsByDiffMonoMass(mods, residue, top_csms_spectrum[i].cross_link.cross_linker_mass, 0.001);

                LOG_DEBUG << "number of modifications fitting the diff mass: " << mods.size() << "\t" << mods << endl;
                bool mod_set = false;
                if (mods.size() > 0) // If several mods have the same diff mass, try to resolve ambiguity by cross-linker name (e.g. DSS and BS3 are different reagents, but have the same result after the reaction)
                {
                  for (Size s = 0; s < mods.size(); ++s)
                  {
                    if (mods[s].hasSubstring(cross_link_name))
                    {
                      LOG_DEBUG << "applied modification: " << mods[s] << endl;
                      seq_alpha.setModification(alpha_pos, mods[s]);
                      mod_set = true;
                      break;
                    }
                  }
                }
                if ( (mods.size() > 0) && (!mod_set) ) // If resolving by name did not work, use any with matching diff mass
                {
                  seq_alpha.setModification(alpha_pos, mods[0]);
                  mod_set = true;
                }
                if (!mod_set) // If no equivalent mono-link exists in the UNIMOD or XLMOD databases, use the given name to construct a placeholder
                {
                  String mod_name = String("unknown mono-link " + cross_link_name + " mass " + String(top_csms_spectrum[i].cross_link.cross_linker_mass));
                  //seq_alpha.setModification(alpha_pos, mod_name);
                  LOG_DEBUG << "unknown mono-link" << endl;
                  ph_alpha.setMetaValue("xl_mod", mod_name);
                }
              }
              else
              {
                ph_alpha.setMetaValue("xl_mod", top_csms_spectrum[i].cross_link.cross_linker_name);
                ph_alpha.setMetaValue("xl_mass", DataValue(top_csms_spectrum[i].cross_link.cross_linker_mass));
              }


              if (top_csms_spectrum[i].cross_link.getType() == ProteinProteinCrossLink::LOOP)
              {
                ph_alpha.setMetaValue("xl_pos2", DataValue(beta_pos));
              }

              vector<PeptideHit> phs;

              ph_alpha.setSequence(seq_alpha);
              ph_alpha.setCharge(precursor_charge);
              ph_alpha.setScore(top_csms_spectrum[i].score);
              ph_alpha.setRank(DataValue(i+1));
              ph_alpha.setMetaValue("xl_chain", "MS:1002509");  // donor (longer, heavier, alphabetically earlier)
              ph_alpha.setMetaValue("xl_pos", DataValue(alpha_pos));
              //ph_alpha.setMetaValue("spec_heavy_RT", spectra[scan_index_heavy].getRT());
              //ph_alpha.setMetaValue("spec_heavy_MZ", spectra[scan_index_heavy].getPrecursors()[0].getMZ());
              ph_alpha.setMetaValue("spectrum_reference", spectra[scan_index].getNativeID());
              //ph_alpha.setMetaValue("spectrum_reference_heavy", spectra[scan_index_heavy].getNativeID());
              ph_alpha.setMetaValue("xl_type", xltype);
              ph_alpha.setMetaValue("xl_rank", DataValue(i + 1));

              ph_alpha.setMetaValue("OpenXQuest:xcorr xlink", top_csms_spectrum[i].xcorrx_max);
              ph_alpha.setMetaValue("OpenXQuest:xcorr common", top_csms_spectrum[i].xcorrc_max);
              ph_alpha.setMetaValue("OpenXQuest:match-odds", top_csms_spectrum[i].match_odds);
              ph_alpha.setMetaValue("OpenXQuest:intsum", top_csms_spectrum[i].int_sum);
              ph_alpha.setMetaValue("OpenXQuest:wTIC", top_csms_spectrum[i].wTIC);

              ph_alpha.setMetaValue("OpenProXL:HyperCommon",top_csms_spectrum[i].HyperCommon);
              ph_alpha.setMetaValue("OpenProXL:HyperXlink",top_csms_spectrum[i].HyperXlink);
              ph_alpha.setMetaValue("OpenProXL:HyperAlpha", top_csms_spectrum[i].HyperAlpha);
              ph_alpha.setMetaValue("OpenProXL:HyperBeta", top_csms_spectrum[i].HyperBeta);
              ph_alpha.setMetaValue("OpenProXL:HyperBoth",top_csms_spectrum[i].HyperBoth);

              ph_alpha.setMetaValue("OpenProXL:AndroCommon",top_csms_spectrum[i].PScoreCommon);
              ph_alpha.setMetaValue("OpenProXL:AndroXlink",top_csms_spectrum[i].PScoreXlink);
              ph_alpha.setMetaValue("OpenProXL:AndroAlpha",top_csms_spectrum[i].PScoreAlpha);
              ph_alpha.setMetaValue("OpenProXL:AndroBeta",top_csms_spectrum[i].PScoreBeta);
              ph_alpha.setMetaValue("OpenProXL:AndroBoth",top_csms_spectrum[i].PScoreBoth);

              ph_alpha.setMetaValue("OpenProXL:HyperABMulti",top_csms_spectrum[i].HyperAlpha * top_csms_spectrum[i].HyperBeta);
              ph_alpha.setMetaValue("OpenProXL:HyperABAddi",top_csms_spectrum[i].HyperAlpha + top_csms_spectrum[i].HyperBeta);
              ph_alpha.setMetaValue("OpenProXL:HyperCXMulti",top_csms_spectrum[i].HyperCommon * top_csms_spectrum[i].HyperXlink);
              ph_alpha.setMetaValue("OpenProXL:HyperCXAddi",top_csms_spectrum[i].HyperCommon + top_csms_spectrum[i].HyperXlink);

              ph_alpha.setMetaValue("OpenProXL:AndroABMulti",top_csms_spectrum[i].PScoreAlpha * top_csms_spectrum[i].PScoreBeta);
              ph_alpha.setMetaValue("OpenProXL:AndroABAddi",top_csms_spectrum[i].PScoreAlpha + top_csms_spectrum[i].PScoreBeta);
              ph_alpha.setMetaValue("OpenProXL:AndroCXMulti",top_csms_spectrum[i].PScoreCommon * top_csms_spectrum[i].PScoreXlink);
              ph_alpha.setMetaValue("OpenProXL:AndroCXAddi",top_csms_spectrum[i].PScoreCommon + top_csms_spectrum[i].PScoreXlink);

              ph_alpha.setFragmentAnnotations(top_csms_spectrum[i].frag_annotations);
              LOG_DEBUG << "Annotations of size " << ph_alpha.getFragmentAnnotations().size() << endl;
              phs.push_back(ph_alpha);

              if (top_csms_spectrum[i].cross_link.getType() == ProteinProteinCrossLink::CROSS)
              {
                ph_beta.setSequence(top_csms_spectrum[i].cross_link.beta);
                ph_beta.setCharge(precursor_charge);
                ph_beta.setScore(top_csms_spectrum[i].score);
                ph_beta.setRank(DataValue(i+1));
                ph_beta.setMetaValue("xl_chain", "MS:1002510"); // receiver
                ph_beta.setMetaValue("xl_pos", DataValue(beta_pos));
                //ph_beta.setMetaValue("spec_heavy_RT", spectra[scan_index_heavy].getRT());
                //ph_beta.setMetaValue("spec_heavy_MZ", spectra[scan_index_heavy].getPrecursors()[0].getMZ());
                ph_beta.setMetaValue("spectrum_reference", spectra[scan_index].getNativeID());
                //ph_beta.setMetaValue("spectrum_reference_heavy", spectra[scan_index_heavy].getNativeID());

                ph_beta.setMetaValue("OpenXQuest:xcorr xlink", top_csms_spectrum[i].xcorrx_max);
                ph_beta.setMetaValue("OpenXQuest:xcorr common", top_csms_spectrum[i].xcorrc_max);
                ph_beta.setMetaValue("OpenXQuest:match-odds", top_csms_spectrum[i].match_odds);
                ph_beta.setMetaValue("OpenXQuest:intsum", top_csms_spectrum[i].int_sum);
                ph_beta.setMetaValue("OpenXQuest:wTIC", top_csms_spectrum[i].wTIC);

                ph_beta.setMetaValue("OpenProXL:HyperCommon",top_csms_spectrum[i].HyperCommon);
                ph_beta.setMetaValue("OpenProXL:HyperXlink",top_csms_spectrum[i].HyperXlink);
                ph_beta.setMetaValue("OpenProXL:HyperAlpha",top_csms_spectrum[i].HyperAlpha);
                ph_beta.setMetaValue("OpenProXL:HyperBeta",top_csms_spectrum[i].HyperBeta);
                ph_beta.setMetaValue("OpenProXL:HyperBoth",top_csms_spectrum[i].HyperBoth);

                ph_beta.setMetaValue("OpenProXL:AndroCommon",top_csms_spectrum[i].PScoreCommon);
                ph_beta.setMetaValue("OpenProXL:AndroXlink",top_csms_spectrum[i].PScoreXlink);
                ph_beta.setMetaValue("OpenProXL:AndroAlpha",top_csms_spectrum[i].PScoreAlpha);
                ph_beta.setMetaValue("OpenProXL:AndroBeta",top_csms_spectrum[i].PScoreBeta);
                ph_beta.setMetaValue("OpenProXL:AndroBoth",top_csms_spectrum[i].PScoreBoth);

                ph_beta.setMetaValue("OpenProXL:HyperABMulti",top_csms_spectrum[i].HyperAlpha * top_csms_spectrum[i].HyperBeta);
                ph_beta.setMetaValue("OpenProXL:HyperABAddi",top_csms_spectrum[i].HyperAlpha + top_csms_spectrum[i].HyperBeta);
                ph_beta.setMetaValue("OpenProXL:HyperCXMulti",top_csms_spectrum[i].HyperCommon * top_csms_spectrum[i].HyperXlink);
                ph_beta.setMetaValue("OpenProXL:HyperCXAddi",top_csms_spectrum[i].HyperCommon + top_csms_spectrum[i].HyperXlink);

                ph_beta.setMetaValue("OpenProXL:AndroABMulti",top_csms_spectrum[i].PScoreAlpha * top_csms_spectrum[i].PScoreBeta);
                ph_beta.setMetaValue("OpenProXL:AndroABAddi",top_csms_spectrum[i].PScoreAlpha + top_csms_spectrum[i].PScoreBeta);
                ph_beta.setMetaValue("OpenProXL:AndroCXMulti",top_csms_spectrum[i].PScoreCommon * top_csms_spectrum[i].PScoreXlink);
                ph_beta.setMetaValue("OpenProXL:AndroCXAddi",top_csms_spectrum[i].PScoreCommon + top_csms_spectrum[i].PScoreXlink);

                phs.push_back(ph_beta);
              }

              peptide_id.setRT(spectrum.getRT());
              peptide_id.setMZ(precursor_mz);
              String specIDs = spectra[scan_index].getNativeID();
              peptide_id.setMetaValue("spectrum_reference", specIDs);
//              peptide_id.setMetaValue("spectrum_reference", spectra[scan_index].getNativeID());
//              peptide_id.setMetaValue("xl_type", xltype); // TODO: needs CV term
//              peptide_id.setMetaValue("xl_rank", DataValue(i + 1));

              peptide_id.setHits(phs);
              peptide_id.setScoreType("OpenXQuest:combined score");

#ifdef _OPENMP
#pragma omp critical (peptides_ids_access)
#endif
              peptide_ids.push_back(peptide_id);

#ifdef _OPENMP
#pragma omp critical (all_top_csms_access)
#endif
              all_top_csms[all_top_csms_current_index][i].peptide_id_index = peptide_ids.size()-1;
            }


            LOG_DEBUG << "###################### Spectrum was processed ################################## \n";
          }
    }
    // end of matching / scoring
    progresslogger.endProgress();

    LOG_DEBUG << "Pre Score maximum: " << pScoreMax << "\t TIC maximum: " << TICMax << "\t wTIC maximum: " << wTICMax << "\t Match-Odds maximum: " << matchOddsMax << endl;
    LOG_DEBUG << "XLink Cross-correlation maximum: " << xcorrxMax << "\t Common Cross-correlation maximum: " << xcorrcMax << "\t Intsum maximum: " << intsumMax << endl;
    LOG_DEBUG << "Total number of matched candidates: " << sumMatchCount << "\t Maximum number of matched candidates to one spectrum pair: " << maxMatchCount << "\t Average: " << sumMatchCount / spectra.size() << endl;

    // Add protein identifications
    PeptideIndexing pep_indexing;
    Param indexing_param = pep_indexing.getParameters();

    // TODO update additional parameters of PeptideIndexing (enzyme etc.)
    String d_prefix = decoy_prefix ? "true" : "false";
    indexing_param.setValue("prefix", d_prefix, "If set, protein accessions in the database contain 'decoy_string' as prefix.");
    indexing_param.setValue("decoy_string", decoy_string, "String that was appended (or prefixed - see 'prefix' flag below) to the accessions in the protein database to indicate decoy proteins.");
    indexing_param.setValue("missing_decoy_action", "warn");
    indexing_param.setValue("enzyme:name", enzyme_name);
    pep_indexing.setParameters(indexing_param);

    pep_indexing.run(fasta_db, protein_ids, peptide_ids);

    // write output
    progresslogger.startProgress(0, 1, "Writing output...");
    if (out_idXML.size() > 0)
    {
      IdXMLFile().store(out_idXML, protein_ids, peptide_ids);
    }
    if (out_mzIdentML.size() > 0)
    {
      MzIdentMLFile().store(out_mzIdentML, protein_ids, peptide_ids);
    }
    if (out_xquest.size() > 0)
    {
      vector<String> input_split_dir;
      vector<String> input_split;
      getStringOption_("in").split(String("/"), input_split_dir);
      input_split_dir[input_split_dir.size()-1].split(String("."), input_split);
      String base_name = input_split[0];

      vector<String> output_split_dir;
      vector<String> output_split;
      Size found;
      found = out_xquest.find_last_of("/\\");
      // TODO "/" is Unix specific
      String matched_spec_xml_name = out_xquest.substr(0, found) + "/" + base_name + "_matched.spec.xml";

      // TODO old version here
//      writeXQuestXML(out_xquest, base_name, peptide_ids, all_top_csms, spectra);

      // TODO use this after RichPeaks are no longer used
      String precursor_mass_tolerance_unit_string = precursor_mass_tolerance_unit_ppm ? "ppm" : "Da";
      String fragment_mass_tolerance_unit_string = fragment_mass_tolerance_unit_ppm ? "ppm" : "Da";
      double cross_link_mass_iso_shift = 0;
      double cross_link_mass_light = cross_link_mass;
      OpenProXLUtils::writeXQuestXML(out_xquest, base_name, peptide_ids, all_top_csms, spectra,
                                                            precursor_mass_tolerance_unit_string, fragment_mass_tolerance_unit_string, precursor_mass_tolerance, fragment_mass_tolerance, fragment_mass_tolerance_xlinks, cross_link_name,
                                                            cross_link_mass_light, cross_link_mass_mono_link, in_fasta, in_decoy_fasta, cross_link_residue1, cross_link_residue2, cross_link_mass_iso_shift, enzyme_name, missed_cleavages);

      writeXQuestXMLSpec(matched_spec_xml_name, base_name, all_top_csms, spectra);
    }
    progresslogger.endProgress();

    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{

  TOPPOpenProXLLF tool;

  return tool.main(argc, argv);
}

