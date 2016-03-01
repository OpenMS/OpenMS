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
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/CHEMISTRY/EnzymesDB.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/ANALYSIS/RNPXL/ModifiedPeptideGenerator.h>
#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/FORMAT/Base64.h>
#include <OpenMS/ANALYSIS/ID/PeptideIndexing.h>

// preprocessing and filtering
#include <OpenMS/FILTERING/TRANSFORMERS/ThresholdMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/NLargest.h>
#include <OpenMS/FILTERING/TRANSFORMERS/WindowMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignment.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGeneratorXLinks.h>

// results
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <iostream>
#include <cmath>
#include <numeric>

#include <boost/unordered_map.hpp>
#include <boost/math/special_functions/binomial.hpp>

using namespace std;
using namespace OpenMS;


#ifdef _OPENMP
#include <omp.h>
#define NUMBER_OF_THREADS (omp_get_num_threads())
#else
#define NUMBER_OF_THREADS (1)
#endif


/**
    @page UTILS_xQuest xQuest

    @brief Perform protein-protein cross-linking experiments.

    <CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ xQuest \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> - </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> - </td>
        </tr>
    </table>
</CENTER>

    @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_xQuest.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_xQuest.html
*/

class TOPPxQuest :
  public TOPPBase
{
public:
  TOPPxQuest() :
    TOPPBase("xQuest", "Tool for protein-protein cross linking using the xQuest algorithm.", false)
  {
  }

  struct CrossLinkSpectrumMatch
  {
    /// structure of the cross-link
    TheoreticalSpectrumGeneratorXLinks::ProteinProteinCrossLink cross_link;

    /// reference to pair of spectra
    Size scan_index_light;
    Size scan_index_heavy;

    /// final score
    double score;

    /// rank among the matches to the same spectrum
    Size rank;

    /// counts, scores and other data for xQuest-like output
    double pre_score;
    double percTIC;
    double wTIC;
    double int_sum;
    double match_odds;
    vector< double > xcorrx;
    double xcorrx_max;
    vector< double > xcorrc;
    double xcorrc_max;
    Size matched_common_alpha;
    Size matched_common_beta;
    Size matched_xlink_alpha;
    Size matched_xlink_beta;

    Size peptide_id_index;
    PeptideIdentification *peptide_id;

    bool operator<(const CrossLinkSpectrumMatch& other) const
    {
      return score < other.score;
    }

    bool operator==(const CrossLinkSpectrumMatch& other) const
    {
      return cross_link.alpha == other.cross_link.alpha &&
           cross_link.beta == other.cross_link.beta &&
           cross_link.cross_link_position == other.cross_link.cross_link_position &&
           score == other.score &&
           pre_score == other.pre_score &&
           percTIC == other.percTIC &&
           wTIC == other.wTIC &&
           int_sum == other.int_sum &&
           match_odds == other.match_odds &&
           xcorrx == other.xcorrx &&
           xcorrx_max == other.xcorrx_max &&
           xcorrc == other.xcorrc &&
           xcorrc_max == other.xcorrc_max &&
           matched_common_alpha == other.matched_common_alpha &&
           matched_common_beta == other.matched_common_beta &&
           matched_xlink_alpha == other.matched_xlink_alpha &&
           matched_xlink_beta == other.matched_xlink_beta &&
           rank == other.rank;
    }

  };

  //TODO: make protected
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
  // TODO: make protected
  static String getxQuestBase64EncodedSpectrum_(const RichPeakSpectrum& spec, String header)
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
      if (spec[i].metaValueExists("z"))
      {
        s += String(spec[i].getMetaValue("z"));
      }
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

    registerInputFile_("consensus", "<file>", "", "Input file containing the linked mass peaks.");
    setValidFormats_("consensus", ListUtils::create<String>("consensusXML"));

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
    registerDoubleOption_("fragment:mass_tolerance", "<tolerance>", 0.2, "Fragment mass tolerance", false);
    registerDoubleOption_("fragment:mass_tolerance_xlinks", "<tolerance>", 0.3, "Fragment mass tolerance for cross-link ions", false);

    StringList fragment_mass_tolerance_unit_valid_strings;
    fragment_mass_tolerance_unit_valid_strings.push_back("ppm");
    fragment_mass_tolerance_unit_valid_strings.push_back("Da");

    registerStringOption_("fragment:mass_tolerance_unit", "<unit>", "Da", "Unit of fragment m", false, false);
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
    registerDoubleOption_("cross_linker:mass_light", "<mass>", 138.0680796, "Mass of the light cross-linker, linking two peptides", false);
    registerDoubleOption_("cross_linker:mass_iso_shift", "<mass>", 12.075321, "Mass of the isotopic shift between the light and heavy linkers", false);
    registerDoubleList_("cross_linker:mass_mono_link", "<mass>", ListUtils::create<double>("156.0786442, 155.0964278"), "Possible masses of the linker, when attached to only one peptide", false);
    // TODO a list of double numbers for different monolink weights
    // TODO different cross-linkers, also bi-functional

    registerTOPPSubsection_("algorithm", "Algorithm Options");
    registerStringOption_("algorithm:candidate_search", "<param>", "index", "Mode used to generate candidate peptides.", false, false);
    StringList candidate_search_modes_strings;
    candidate_search_modes_strings.push_back("index");
    candidate_search_modes_strings.push_back("enumeration");
    setValidStrings_("algorithm:candidate_search", candidate_search_modes_strings);

    // output file
    registerOutputFile_("out", "<file>", "", "Result file\n");
    setValidFormats_("out", ListUtils::create<String>("xml"));

    registerOutputFile_("out_idXML", "<file>", "", "output file ");
    setValidFormats_("out_idXML", ListUtils::create<String>("idXML"));
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

  void nlargest_filter_rich(RichPeakSpectrum & spectrum, Size peakcount)
  {
    if (spectrum.size() <= peakcount) return;

    // sort by reverse intensity
    spectrum.sortByIntensity(true);

    // keep the n largest peaks if more than n are present
    spectrum.resize(peakcount);
  }

  void preprocessSpectra_(PeakMap& exp, RichPeakMap& rich_exp)
  {
    // filter MS2 map
    // remove 0 intensities
    ThresholdMower threshold_mower_filter;
    threshold_mower_filter.filterPeakMap(exp);
 
    Normalizer normalizer;
    normalizer.filterPeakMap(exp);

    // sort by rt
    exp.sortSpectra(false);

    // filter settings
//    WindowMower window_mower_filter;
//    Param filter_param = window_mower_filter.getParameters();
//    filter_param.setValue("windowsize", 100.0, "The size of the sliding window along the m/z axis.");
//    filter_param.setValue("peakcount", 20, "The number of peaks that should be kept.");
//    filter_param.setValue("movetype", "jump", "Whether sliding window (one peak steps or jumping window window size steps) should be used.");
//    window_mower_filter.setParameters(filter_param);
    NLargest nlargest_filter = NLargest(250);   // De-noising in xQuest: Dynamic range = 1000, 250 most intense peaks?
  
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (SignedSize exp_index = 0; exp_index < static_cast<SignedSize>(exp.size()); ++exp_index)
    {
      // sort by mz, for deisotoping and window_mower
      exp[exp_index].sortByPosition();     
//      nlargest_filter.filterSpectrum(exp[exp_index]);

      // Deisotope, compute charges here (no 0 intensities, sorted by mz)
      //rich_exp.addSpectrum(deisotopeAndSingleChargeMSSpectrum(exp[exp_index], 1, 6, 0.3, true));
      //Params: spectrum, Int min_charge, Int max_charge, double fragment_tolerance, bool fragment_unit_ppm, bool keep_only_deisotoped = false, Size min_isopeaks = 3, Size max_isopeaks = 10, bool make_single_charged = true
  
      // remove noise, TODO window_mower not used in xQuest, is it necessary?
//      window_mower_filter.filterPeakSpectrum(exp[exp_index]);
      //nlargest_filter_rich(rich_exp[exp_index], 250);
      nlargest_filter.filterSpectrum(exp[exp_index]);
  
      // sort (nlargest changes order)
//      exp[exp_index].sortByPosition();
      rich_exp.addSpectrum(makeRichPeakSpectrum(exp[exp_index], true));
//      cout << "WHY? Exp Spec: " << exp[exp_index][50].getMZ() << "\t" << rich_exp[exp_index][50].getMZ() << endl;
      rich_exp[exp_index].sortByPosition();
     }
   }

  // xQuest, fast pre-Score for x-links (type 2)
  // required: numbers of peaks for each chain, and how many of them were matched
  float preScore(Size matchedAlpha, Size ionsAlpha, Size matchedBeta, Size ionsBeta)
  {
    if ( (ionsAlpha > 0) && (ionsBeta > 0) )
    {
      float result = sqrt((static_cast<float>(matchedAlpha) / static_cast<float>(ionsAlpha)) * (static_cast<float>(matchedBeta) / static_cast<float>(ionsBeta)));
      return result;
    } else
    {
      return 0.0;
    }
  }

  // xQuest, fast pre-Score for Mono links and Loop links (type 0 and type 1)
  float preScore(Size matchedAlpha, Size ionsAlpha)
  {
    if (ionsAlpha > 0)
    {
      float result = static_cast<float>(matchedAlpha) / static_cast<float>(ionsAlpha);
      return result;
    } else
    {
      return 0.0;
    }
  }

  // Statistics/Combinatorics functions for match-odds score
  // Standard cumulative binomial distribution
  double cumulativeBinomial(Size n, Size k, double p)
  {
    double p_cumul = 0.0;
    if (p < 1e-99) return static_cast<double>(k == 0); //  (not true/false, but 1/0 as probability)
    if (1 - p < 1e-99) return static_cast<double>(k != n); //
    if (k > n)  return 1.0;

    for (Size j = 0; j < k; j++)
    {
      double coeff = boost::math::binomial_coefficient<double>(static_cast<unsigned int>(n), static_cast<unsigned int>(j));
      p_cumul += coeff * pow(p,  j) * pow((1-p), (n-j));
    }

    // match-odds score becomes INFINITY for p_cumul = 1, p_cumul might reach 1 because of insufficient precision, solved by using largest value smaller than 1
    if (p_cumul >= 1.0)
    {
      p_cumul = nexttoward(1.0, 0.0);
    }

    return p_cumul;
  }


  // A priori probability of a random match given info about the theoretical spectrum, (range = maxmz - minmz)
  double aPrioriProb(double ms2_tolerance, Size n_peaks, double range, Size n_charges)
  {
    double a_priori_probability= (1 - ( pow( (1 - 2 * ms2_tolerance/ (0.5 * range)),  (n_peaks / n_charges))));
    return a_priori_probability;
  }
  
  //TODO more general spectra types?
  // Cross-correlation, with shifting the second spectrum from -maxshift to +maxshift of tolerance bins (Tolerance in Da, basically a constant binsize)
  vector< double > xCorrelation(const RichPeakSpectrum & spec1, const RichPeakSpectrum & spec2, Int maxshift, double tolerance)
  {
    double maxionsize = max(spec1[spec1.size()-1].getMZ(), spec2[spec2.size()-1].getMZ());
    //cout << "xcorr Maxionssize: " << maxionsize << endl;
    Int table_size = ceil(maxionsize / tolerance)+1;
    //cout << "xcorr table_size: " << table_size << endl;
    vector< double > ion_table1(table_size, 0);
    vector< double > ion_table2(table_size, 0);

    // Build tables of the same size, each bin has the size of the tolerance
    for (Size i = 0; i < spec1.size(); ++i)
    {
      Size pos = static_cast<Size>(ceil(spec1[i].getMZ() / tolerance));
      //cout << "xcorr Table1 pos: " << pos << endl;
      // TODO this line leads to using real intensities
//      ion_table1[pos] = spec1[i].getIntensity();

      // TODO this line leads to using intensities normalized to 10
      ion_table1[pos] = 10.0;
      //cout << "xcorr Table1 Inte: " << spec1[i].getIntensity() << endl;
    }
    for (Size i = 0; i < spec2.size(); ++i)
    {
      Size pos =static_cast<Size>(ceil(spec2[i].getMZ() / tolerance));
      //cout << "xcorr Table2 pos: " << pos << endl;
      // TODO this line leads to using real intensities
//      ion_table2[pos] = spec2[i].getIntensity();

      // TODO this line leads to using intensities normalized to 10
      ion_table2[pos] = 10.0;
      //cout << "xcorr Table2 Inte: " << spec2[i].getIntensity() << endl;
    }
    //cout << "xcorr Tables done" << endl;

    // Compute means
    double mean1 = (accumulate(ion_table1.begin(), ion_table1.end(), 0.0)) / table_size;
    double mean2 = (accumulate(ion_table2.begin(), ion_table2.end(), 0.0)) / table_size;
    //cout << "xcorr means: " << mean1 << "\t" << mean2 << endl;

    // Compute denominator
    double s1 = 0;
    double s2 = 0;
    for (Int i = 0; i < table_size; ++i)
    {
      s1 += pow((ion_table1[i] - mean1), 2);
      s2 += pow((ion_table2[i] - mean2), 2);
    }
    double denom = sqrt(s1 * s2);
    //cout << "xcorr Denominator: " << denom << endl;

    // Calculate correlation for each shift
    vector< double > results(maxshift * 2 + 1, 0);
    for (Int shift = -maxshift; shift <= maxshift; ++shift)
    {
      double s = 0;
      for (Int i = 0; i < table_size; ++i)
      {
        Int j = i + shift;
        //cout << "xcorr i: " << i << "\t shift: " << shift << "\t j: " << j << endl;
        if ( (j >= 0) && (j < table_size))
        {
          s += (ion_table1[i] - mean1) * (ion_table2[j] - mean2);
         // cout << "XCORR S: " << s << endl;
        }
      }
      if (denom > 0)
      {
        results[shift + maxshift] = s / denom;
      }
    }
    //cout << "xcorr s/denom vector: " << results << endl;
    return results;
  }
     

  struct PreprocessedPairSpectra_
  {
    // pre-initialize so we can simply std::swap the spectra (no synchronization in multi-threading context needed as we get no reallocation of the PeakMapreprocessed_pair_spectra. 
    RichPeakMap spectra_light_different; // peaks in light spectrum after common peaks have been removed
    RichPeakMap spectra_heavy_different; // peaks in heavy spectrum after common peaks have been removed
    RichPeakMap spectra_heavy_to_light; // heavy peaks transformed to light ones and after common peaks have been removed
    RichPeakMap spectra_common_peaks; // merge spectrum of common peaks (present in both spectra)
    RichPeakMap spectra_xlink_peaks; // Xlink peaks in the light spectrum (common peaks between spectra_light_different and spectra heavy_to_light)

    PreprocessedPairSpectra_(Size size)
    {
      for (Size i = 0; i != size; ++i)
      {
        spectra_light_different.addSpectrum(RichPeakSpectrum());
        spectra_heavy_different.addSpectrum(RichPeakSpectrum());
        spectra_heavy_to_light.addSpectrum(RichPeakSpectrum());
        spectra_common_peaks.addSpectrum(RichPeakSpectrum());
        spectra_xlink_peaks.addSpectrum(RichPeakSpectrum());
      }
    }  
  };


  // create common / shifted peak spectra for all pairs
  PreprocessedPairSpectra_ preprocessPairs_(const RichPeakMap& spectra, const vector< pair<Size, Size> >& spectrum_pairs, const double cross_link_mass_light, const double cross_link_mass_iso_shift)
  {
    PreprocessedPairSpectra_ preprocessed_pair_spectra(spectrum_pairs.size());

    Size max_charge_xlink = 6;
 
#ifdef _OPENMP
#pragma omp parallel for
#endif
//    for (SignedSize scan_index = 0; scan_index < static_cast<SignedSize>(spectra.size()); ++scan_index)
    for (Size pair_index = 0; pair_index < spectrum_pairs.size(); ++pair_index)
    {
      // assume that current MS2 corresponds to a light spectrum
      Size scan_index = spectrum_pairs[pair_index].first;
      const RichPeakSpectrum& spectrum_light = spectra[scan_index];

      // map light to heavy
//      map<Size, Size>::const_iterator scan_index_light_it = map_light_to_heavy.find(scan_index);

//      if (scan_index_light_it != map_light_to_heavy.end())
//      {
        const Size scan_index_heavy = spectrum_pairs[pair_index].second;
        const RichPeakSpectrum& spectrum_heavy = spectra[scan_index_heavy];
         vector< pair< Size, Size > > matched_fragments_without_shift;
        //ms2_aligner.getSpectrumAlignment(matched_fragments_without_shift, spectrum_light, spectrum_heavy);
//        cout << "Common matching, spectrum: " << scan_index << endl;
        getSpectrumIntensityMatching(matched_fragments_without_shift, spectrum_light, spectrum_heavy, 0.2, false, 0.3);

        // different fragments may carry light or heavy cross-linker.             
        RichPeakSpectrum spectrum_heavy_different;
        RichPeakSpectrum spectrum_light_different;

        // TODO: maybe speed this up - can be done in linear time
        for (Size i = 0; i != spectrum_light.size(); ++i)
        {
          bool found = false;
          for (Size j = 0; j != matched_fragments_without_shift.size(); ++j)
          {
            if (matched_fragments_without_shift[j].first == i) { found = true; break; }
          }
          if (!found)
          {
            spectrum_light_different.push_back(spectrum_light[i]);
          }
        }
        for (Size i = 0; i != spectrum_heavy.size(); ++i)
        {
          bool found = false;
          for (Size j = 0; j != matched_fragments_without_shift.size(); ++j)
          {
            if (matched_fragments_without_shift[j].second == i) { found = true; break; }
          }
          if (!found)
          {
            spectrum_heavy_different.push_back(spectrum_heavy[i]);
          }
        }

        // transform by m/z difference between unlabeled and labeled cross-link to make heavy and light comparable.
        // assume different charged MS2 fragments
        RichPeakSpectrum spectrum_heavy_to_light;
        RichPeakSpectrum xlink_peaks;
        for (Size charge = 2; charge <= max_charge_xlink; ++charge)
          {
            spectrum_heavy_to_light.clear(false);
            double mass_shift = cross_link_mass_iso_shift / charge;
            // transform heavy spectrum to light spectrum, collect peaks with the current charge into the spectrum
            for (Size i = 0; i != spectrum_heavy_different.size(); ++i)
            {
              RichPeak1D p = spectrum_heavy_different[i];
              p.setMZ(p.getMZ() - mass_shift);
              spectrum_heavy_to_light.push_back(p);
            }

            // align potentially shifted peaks from light MS2 with potentially shifted peaks from heavy (after transformation to resemble the light MS2)
            // matching fragments are potentially carrying the cross-linker
            vector< pair< Size, Size > > matched_fragments_with_shift;

            //ms2_aligner_xlinks.getSpectrumAlignment(matched_fragments_with_shift, spectrum_light_different, spectrum_heavy_to_light);
            spectrum_heavy_to_light.sortByPosition();
            if (spectrum_light_different.size() > 0 && spectrum_heavy_to_light.size() > 0) {
//            cout << "Xlink matching, spectrum: " << scan_index << "\tcharge: " << charge << endl;
              getSpectrumIntensityMatching(matched_fragments_with_shift, spectrum_light_different, spectrum_heavy_to_light, 0.3, false, 0.3);

              for (Size i = 0; i != matched_fragments_with_shift.size(); ++i)
              {
                xlink_peaks.push_back(spectrum_light_different[matched_fragments_with_shift[i].first]);
              }
              // make sure to not include the same peaks more than once with a different charge
              for (Size i = 0; i != matched_fragments_with_shift.size(); ++i)
              {
                spectrum_light_different.erase(spectrum_light_different.begin() + matched_fragments_with_shift[i].first);
                Size heavy_index = spectrum_heavy_different.findNearest(spectrum_heavy_to_light[matched_fragments_with_shift[i].second].getMZ() + mass_shift);
                spectrum_heavy_different.erase(spectrum_heavy_different.begin() + heavy_index);
              }
            }
          }


#ifdef DEBUG_XQUEST
        cout << "Common peaks: " << matched_fragments_without_shift.size() << " different peaks: " << spectrum_light.size() - matched_fragments_without_shift.size() << ", " << spectrum_heavy.size() - matched_fragments_without_shift.size() << endl;
        cout << "Matched shifted peaks: " << matched_fragments_with_shift.size() << " unexplained peaks: " << spectrum_light_different.size() - matched_fragments_with_shift.size() << ", " << spectrum_heavy_to_light.size() - matched_fragments_with_shift.size() << endl;
#endif
        // generate common peaks spectrum TODO: check if only light / only heavy or e.g. mean of both peaks are used here
        RichPeakSpectrum common_peaks;
        for (Size i = 0; i != matched_fragments_without_shift.size(); ++i)
        {
          common_peaks.push_back(spectrum_light[matched_fragments_without_shift[i].first]);
        }
#ifdef DEBUG_XQUEST
        cout << "Peaks to match: " << common_peaks.size() << endl;
#endif

        swap(preprocessed_pair_spectra.spectra_light_different[pair_index], spectrum_light_different);
        swap(preprocessed_pair_spectra.spectra_heavy_different[pair_index], spectrum_heavy_different);
        swap(preprocessed_pair_spectra.spectra_heavy_to_light[pair_index], spectrum_heavy_to_light);
        swap(preprocessed_pair_spectra.spectra_common_peaks[pair_index], common_peaks);
        swap(preprocessed_pair_spectra.spectra_xlink_peaks[pair_index], xlink_peaks);

        preprocessed_pair_spectra.spectra_common_peaks[pair_index].setPrecursors(spectrum_light.getPrecursors());
        preprocessed_pair_spectra.spectra_xlink_peaks[pair_index].setPrecursors(spectrum_light.getPrecursors());

        preprocessed_pair_spectra.spectra_light_different[pair_index].sortByPosition();
        preprocessed_pair_spectra.spectra_heavy_different[pair_index].sortByPosition();
        preprocessed_pair_spectra.spectra_heavy_to_light[pair_index].sortByPosition();
        preprocessed_pair_spectra.spectra_common_peaks[pair_index].sortByPosition();
        preprocessed_pair_spectra.spectra_xlink_peaks[pair_index].sortByPosition();

        // Debug support output
        /*
        cout << "spectrum_light_different: " << preprocessed_pair_spectra.spectra_light_different[pair_index].size() << endl;
        cout << "spectrum_heavy_different: " << preprocessed_pair_spectra.spectra_heavy_different[pair_index].size() << endl;
        cout << "spectrum_heavy_to_light alignment: " << preprocessed_pair_spectra.spectra_heavy_to_light[pair_index].size() << endl;
        cout << "spctrum_common_peaks: " << preprocessed_pair_spectra.spectra_common_peaks[pair_index].size() << endl;
        cout << "spectrum_xlink_peaks: " << preprocessed_pair_spectra.spectra_xlink_peaks[pair_index].size() << endl;
        */
//      }
    }
    return preprocessed_pair_spectra;
  }

  struct OPENMS_DLLAPI MatchedIonCount
  {
    // Note: code is optimize for theo_spectrum to contain less peaks than the exp_spectrum
    static Size compute(double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, const MSSpectrum<RichPeak1D>& theo_spectrum, const MSSpectrum<Peak1D>& exp_spectrum)
    {
      Size matches(0), i(0), j(0);
      const Size count_exp_spectrum = exp_spectrum.size();
      const Size count_theo_spectrum = theo_spectrum.size();
      do
      { 
        // advance i until in or right of tolerance window
        while (i < count_exp_spectrum && j < count_theo_spectrum)
        {
          const double tolerance_Th = fragment_mass_tolerance_unit_ppm ? theo_spectrum[j].getMZ() * 1e-6 * fragment_mass_tolerance : fragment_mass_tolerance;
          if (exp_spectrum[i].getMZ() < (theo_spectrum[j].getMZ() - tolerance_Th)) 
          {
            ++i;
          }
          else
          {
            break;
          }
        }

        if (i == count_exp_spectrum || j == count_theo_spectrum) return matches;

        // check if in tolerance window (could also be to the right of it) and count as match
        const double tolerance_Th = fragment_mass_tolerance_unit_ppm ? theo_spectrum[j].getMZ() * 1e-6 * fragment_mass_tolerance : fragment_mass_tolerance;
        if (exp_spectrum[i].getMZ() < (theo_spectrum[j].getMZ() + tolerance_Th))
        {
          // matched i and j (can't match again)
          ++matches; 
          ++j; 
          ++i; 
        } 

        if (i == count_exp_spectrum || j == count_theo_spectrum) return matches;

        ++j;
      } while (true);  
      return matches;
    }
  };

  struct OPENMS_DLLAPI HashGrid1D
  {
    // bucket size should be 2.0 * fragment mass tolerance
    // this ensures that two neighboring buckets cover all position possible given the fragment mass tolerance
    // Note: min and max are the absolute boundaries so no (position < min) or (position > max) is allowed
    HashGrid1D(double min, double max, double bucket_size) : 
      min_(min), 
      max_(max), 
      bucket_size_(bucket_size)
    {
      Size n_buckets = ceil((max - min) / bucket_size) + 1;
      h_.rehash(n_buckets);
    }

    void insert(double position, AASequence*& v)
    {     
      if (position < min_ || position > max_) 
      { 
        cerr << "Trying to add element left or right of allowed bounds. (min, max, position): " << min_ << ", " << max_ << ", " << position << endl;
        return;
      }

      const double bucket_index = (position - min_) / bucket_size_;
      h_.insert(make_pair(bucket_index, v));
    }

    vector<AASequence*> get(double position, Size max_elements = numeric_limits<Size>::max())
    {
      if (position < min_ || position > max_) 
      {
        cerr << "Trying to access element left or right of allowed bounds. (min, max, position): " << min_ << ", " << max_ << ", " << position << endl;
        if (position < min_) return vector<AASequence*>();
        if (position > max_) return vector<AASequence*>();
      }

      vector<AASequence*> elements;

      double bucket_pos = (position - min_) / bucket_size_;
      int bucket_index = static_cast<int>(bucket_pos);

      boost::unordered_multimap<int, AASequence*>::const_iterator it = h_.find(bucket_index);

      while (it != h_.end() && it->first == bucket_index && elements.size() < max_elements)
      {
        elements.push_back(it->second);
        ++it;
      }

      // add elements from neighboring buckets
//      if (bucket_pos - bucket_index <= 0.5)
//      {
//        it = h_.find(bucket_index - 1);
//        while (it != h_.end() && it->first == (bucket_index - 1) && elements.size() < max_elements)
//        {
//          elements.push_back(it->second);
//          ++it;
//        }
//      }
//      else
//      {
//        it = h_.find(bucket_index + 1);
//        while (it != h_.end() && it->first == (bucket_index + 1) && elements.size() < max_elements)
//        {
//          elements.push_back(it->second);
//          ++it;
//        }
//      }
      return elements;
    }

    boost::unordered_multimap<int, AASequence*> h_; // map bucket number to AASequences
    double min_;
    double max_;
    double bucket_size_;
  };

  // TODO MEMORY, use different method to save those AASequences? (Strings, Indices?)
  static multimap<double, pair<const AASequence*, const AASequence*> > enumerateCrossLinksAndMasses_(const multimap<StringView, AASequence>&  peptides, double cross_link_mass_light, DoubleList cross_link_mass_mono_link)
  {
    multimap<double, pair<const AASequence*, const AASequence*> > mass_to_candidates;
    for (map<StringView, AASequence>::const_iterator a = peptides.begin(); a != peptides.end(); ++a)
    {
      if (a->second.toUnmodifiedString().find("K") >= a->second.size()-1)
      {
        continue;
      }

      // generate mono-links
      for (Size i = 0; i < cross_link_mass_mono_link.size(); i++)
      {
        double cross_link_mass = a->second.getMonoWeight() + cross_link_mass_mono_link[i];
        // Make sure it is clear this is a monolink, (is an empty pointer an option?)
//        mass_to_candidates.insert(make_pair(cross_link_mass, make_pair<const AASequence*, const AASequence*>(&(a->second), ???)));
      }

      for (map<StringView, AASequence>::const_iterator b = a; b != peptides.end(); ++b)
      {

        if (b->second.toUnmodifiedString().find("K") >= b->second.size()-1)
        {
          continue;
        }

        // mass peptide1 + mass peptide2 + cross linker mass - cross link loss
        double cross_link_mass = a->second.getMonoWeight() + b->second.getMonoWeight() + cross_link_mass_light;
        mass_to_candidates.insert(make_pair(cross_link_mass, make_pair<const AASequence*, const AASequence*>(&(a->second), &(b->second))));
//        mass_to_candidates.insert(make_pair(cross_link_mass, make_pair<String, String>(a->second.toString(), b->second.toString())));
      }
    }

    return mass_to_candidates;
  }

  RichPeakSpectrum makeRichPeakSpectrum(PeakSpectrum spectrum, bool is_common_or_xlink_spectrum)
  {
    RichPeakSpectrum rich_spectrum;
    //vector<Precursor> precursors;
    //Precursor pc = spectrum.getPrecursors()[0];
    //pc.setMZ(spectrum.getPrecursors()[0].getMZ());
    //pc.setIntensity(spectrum.getPrecursors()[0].getIntensity());
    //pc.setCharge(getPrecursors()[0].getCharge());
    //pc.setMZ(15);
    //pc.setIntensity(30);
    //precursors.push_back(pc);
    //rich_spectrum.setPrecursors(precursors);

    rich_spectrum.setPrecursors(spectrum.getPrecursors());
    rich_spectrum.setRT(spectrum.getRT());

    for(Size i = 0; i < spectrum.size(); ++i)
    {
      Peak1D peak = spectrum[i];
      RichPeak1D p;
      p.setMZ(peak.getMZ());
      p.setIntensity(peak.getIntensity());
      if (is_common_or_xlink_spectrum)
      {
        p.setMetaValue("z", 0); // TODO Where to get the actual charge?
      }
      rich_spectrum.push_back(p);
    }
    return rich_spectrum;
  }

    // spectrum must not contain 0 intensity peaks and must be sorted by m/z
//    RichPeakSpectrum deisotopeAndSingleChargeMSSpectrum(PeakSpectrum& old_spectrum, Int min_charge, Int max_charge, double fragment_tolerance, bool fragment_unit_ppm, bool keep_only_deisotoped = false, Size min_isopeaks = 3, Size max_isopeaks = 10, bool make_single_charged = false)
//    {

//      // Input Spectrum originally called "in"
//      //PeakSpectrum old_spectrum = in;
//      RichPeakSpectrum out;
//      vector<Size> mono_isotopic_peak(old_spectrum.size(), 0);
//      if (old_spectrum.empty())
//      {
//        return out;
//      }


//      // determine charge seeds and extend them
//      vector<Int> features(old_spectrum.size(), -1);
//      Int feature_number = 0;

//      for (Size current_peak = 0; current_peak != old_spectrum.size(); ++current_peak)
//      {
//        double current_mz = old_spectrum[current_peak].getPosition()[0];

//        for (Int q = max_charge; q >= min_charge; --q)   // important: test charge hypothesis from high to low
//        {
//          // try to extend isotopes from mono-isotopic peak
//          // if extension larger then min_isopeaks possible:
//          //   - save charge q in mono_isotopic_peak[]
//          //   - annotate all isotopic peaks with feature number
//          if (features[current_peak] == -1)   // only process peaks which have no assigned feature number
//          {
//            bool has_min_isopeaks = true;
//            vector<Size> extensions;
//            for (Size i = 0; i < max_isopeaks; ++i)
//            {
//              double expected_mz = current_mz + i * Constants::C13C12_MASSDIFF_U / q;
//              Size p = old_spectrum.findNearest(expected_mz);
//              double tolerance_dalton = fragment_unit_ppm ? fragment_tolerance * old_spectrum[p].getPosition()[0] * 1e-6 : fragment_tolerance;
//              if (fabs(old_spectrum[p].getPosition()[0] - expected_mz) > tolerance_dalton)   // test for missing peak
//              {
//                if (i < min_isopeaks)
//                {
//                  has_min_isopeaks = false;
//                }
//                break;
//              }
//              else
//              {
//                // TODO: include proper averagine model filtering. for now start at the second peak to test hypothesis
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

//                // averagine check passed
//                extensions.push_back(p);
//              }
//            }

//            if (has_min_isopeaks)
//            {
//              //cout << "min peaks at " << current_mz << " " << " extensions: " << extensions.size() << endl;
//              mono_isotopic_peak[current_peak] = q;
//              for (Size i = 0; i != extensions.size(); ++i)
//              {
//                features[extensions[i]] = feature_number;
//              }
//              feature_number++;
//            }
//          }
//        }
//      }

//// OLD VERSION, creating deisotoped PeakSpectrum
////      in.clear(false);
////      for (Size i = 0; i != old_spectrum.size(); ++i)
////      {
////        Int z = mono_isotopic_peak[i];
////        if (keep_only_deisotoped)
////        {
////          if (z == 0)
////          {
////            continue;
////          }

////          // if already single charged or no decharging selected keep peak as it is
////          if (!make_single_charged)
////          {
////            in.push_back(old_spectrum[i]);
////          }
////          else
////          {
////            Peak1D p = old_spectrum[i];
////            p.setMZ(p.getMZ() * z - (z - 1) * Constants::PROTON_MASS_U);
////            in.push_back(p);
////          }
////        }
////        else
////        {
////          // keep all unassigned peaks
////          if (features[i] < 0)
////          {
////            in.push_back(old_spectrum[i]);
////            continue;
////          }

////          // convert mono-isotopic peak with charge assigned by deisotoping
////          if (z != 0)
////          {
////            if (!make_single_charged)
////            {
////              in.push_back(old_spectrum[i]);
////            }
////            else
////            {
////              Peak1D p = old_spectrum[i];
////              p.setMZ(p.getMZ() * z - (z - 1) * Constants::PROTON_MASS_U);
////              in.push_back(p);
////            }
////          }
////        }
////      }

//      // NEW VERSION, creating RichPeakSpectrum containing charges
//      //out.clear(false);

//      for (Size i = 0; i != old_spectrum.size(); ++i)
//      {
//        Int z = mono_isotopic_peak[i];
//        if (keep_only_deisotoped)
//        {
//          if (z == 0)
//          {
//            continue;
//          }

//          // if already single charged or no decharging selected keep peak as it is
//          if (!make_single_charged)
//          {
//            RichPeak1D p;
//            p.setMZ(old_spectrum[i].getMZ());
//            p.setIntensity(old_spectrum[i].getIntensity());
//            p.setMetaValue("z", z);
//            out.push_back(p);
//          }
//          else
//          {
//            RichPeak1D p;
//            p.setIntensity(old_spectrum[i].getIntensity());
//            p.setMZ(old_spectrum[i].getMZ() * z - (z - 1) * Constants::PROTON_MASS_U);
//            p.setMetaValue("z", 1);
//            out.push_back(p);
//          }
//        }
//        else
//        {
//          // keep all unassigned peaks
//          if (features[i] < 0)
//          {
//            RichPeak1D p;
//            p.setMZ(old_spectrum[i].getMZ());
//            p.setIntensity(old_spectrum[i].getIntensity());
//            p.setMetaValue("z", 0);
//            out.push_back(p);
//            continue;
//          }

//          // convert mono-isotopic peak with charge assigned by deisotoping
//          if (z != 0)
//          {
//            if (!make_single_charged)
//            {
//              RichPeak1D p;
//              p.setMZ(old_spectrum[i].getMZ());
//              p.setIntensity(old_spectrum[i].getIntensity());
//              p.setMetaValue("z", z);
//              out.push_back(p);
//            }
//            else
//            {
//              RichPeak1D p;
//              p.setMZ(old_spectrum[i].getMZ());
//              p.setIntensity(old_spectrum[i].getIntensity());
//              p.setMetaValue("z", z);
//              out.push_back(p);
//            }
//          }
//        }
//      }

////      vector<Precursor> precursors;
////      Precursor pc;
////      pc.setMZ(712.0402832);
////      pc.setCharge(3);
////      precursors.push_back(pc);
//      out.setPrecursors(old_spectrum.getPrecursors());
//      out.setRT(old_spectrum.getRT());
//      out.sortByPosition();
//      return out;
//    }

    // Spectrum Alignment function adapted from SpectrumAlignment.h, intensity_cutoff: 0 for not considering, 0.3 = lower intentity has to be at least 30% of higher intensity (< 70% difference)
    // TODO Copies spectra, so that nearest peaks with dissimilar intensity can be removed to try with 2. nearest until out of tolerance range, make more efficient?
    void getSpectrumAlignment(std::vector<std::pair<Size, Size> > & alignment, const RichPeakSpectrum & s1, RichPeakSpectrum s2, double tolerance, bool relative_tolerance, double intensity_cutoff = 0.0) const
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

                  // NEW: check for similar intensity values
                  double intensity1(s1[i - 1].getIntensity());
                  double intensity2(s2[j - 1].getIntensity());
                  bool diff_int_clear = (min(intensity1, intensity2) / max(intensity1, intensity2)) > intensity_cutoff;


                  //cout << "Intensity test: " << intensity1 << "\t" << intensity2 << "\t" << (max(intensity1, intensity2) / min(intensity1, intensity2)) << "\t" << diff_int_clear << endl;
                  if (score_align <= score_up && score_align <= score_left && diff_align < tolerance && diff_int_clear)
                    {
                      matrix[i][j] = score_align;
                      traceback[i][j] = std::make_pair(i - 1, j - 1);
                      last_i = i;
                      last_j = j;
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

              // NEW: check for similar intensity values
              double intensity1(s1[i - 1].getIntensity());
              double intensity2(s2[j - 1].getIntensity());
              bool diff_int_clear = (min(intensity1, intensity2) / max(intensity1, intensity2)) > intensity_cutoff;

              // found peak match
              if (std::abs(theo_mz - exp_mz) <= max_dist_dalton && diff_int_clear)
                {
                  alignment.push_back(std::make_pair(i, j));
                } else if (std::abs(theo_mz - exp_mz) <= max_dist_dalton && !diff_int_clear)
                {
                  s2.erase(s2.begin() + s2.findNearest(theo_mz));
                  if (i > 0)
                  {
                    i--;
                   }
                }
            }
        }
    }

    RichPeakSpectrum getToleranceWindowPeaks(RichPeakSpectrum spec, double mz, double tolerance, bool relative_tolerance) const
    {
//      RichPeakSpectrum spec = spectrum;
      RichPeakSpectrum peaks;
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
//        cout << "Peak_mz " << peak_mz <<  "\tin range: " << (dist < max_dist) << endl;
        if (dist <= max_dist)
        {
          peaks.push_back(spec[index]);
//          cout << "Added Peak at " << peak_mz << "\twith Intensity: " << spec[index].getIntensity() << endl;
          spec.erase(spec.begin() + index);
        } else
        {
          inside = false;
        }
      }
//      spec.clear(true);
      return peaks;
    }

    // Spectrum Alignment function adapted from SpectrumAlignment.h, intensity_cutoff: 0 for not considering, 0.3 = lower intentity has to be at least 30% of higher intensity (< 70% difference)
    void getSpectrumIntensityMatching(std::vector<std::pair<Size, Size> > & alignment, RichPeakSpectrum s1, RichPeakSpectrum s2, double tolerance, bool relative_tolerance, double intensity_cutoff = 0.0) const
    {
      if (!s2.isSorted())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Input to SpectrumAlignment is not sorted!");
      }
//      RichPeakSpectrum s1 = spec1;
//      RichPeakSpectrum s2 = spec2;
      RichPeakSpectrum spec1 = s1;

      // clear result
      alignment.clear();

      s1.sortByIntensity(true);
      for (Size i = 0; i != s1.size(); ++i)
      {
        const double& s1_mz = s1[i].getMZ();
//        cout << "Search peak: " << s1_mz << "\t with Intensity: " << s1[i].getIntensity() << endl;
        RichPeakSpectrum peaks = getToleranceWindowPeaks(s2, s1_mz, tolerance, relative_tolerance);
//        cout << "NOP: " << peaks.size() << endl;

        if (peaks.size() > 0)
        {
          peaks.sortByIntensity();
          double intensity1(s1[i].getIntensity());

          for (Size j = 0; j < peaks.size(); ++j)
          {
            // check for similar intensity values
            double intensity2(peaks[j].getIntensity());
            bool diff_int_clear = (min(intensity1, intensity2) / max(intensity1, intensity2)) > intensity_cutoff;

//            cout << "Most intense peak mz: " << peaks[j].getMZ() << "\twith Intensity: " << peaks[j].getIntensity() << "\tIntensity diff cleared: " << diff_int_clear << endl;
            // found peak match. if intensity similar enough, update alignment and remove peak, so that it will not get matched again to another peak
            if (diff_int_clear)
            {
              Size s2_index = s2.findNearest(peaks[j].getMZ());
              Size s1_index = spec1.findNearest(s1_mz);
//              if (spec1[s1_index] == s1[i] && s2[s2_index] == peaks[j])
//              {
                alignment.push_back(std::make_pair(s1_index, s2_index));
//              }
              s2.erase(s2.begin() + s2_index);
//              s2.sortByPosition();
              break;
            } else
            {
              // erase the most intense peak, because it does not fulfill the similarity constraint and try with the rest of the peaks
              peaks.erase(peaks.begin() + j);
              --j;
            }
          }
        }
//        cout << "######################################################################" << endl;
//        s1.sortByPosition();
      }
//        s1.clear(true);
//        s2.clear(true);
    }

    // Write xQuest.xml output
    void writeXQuestXML(const String& out_file, const vector< vector< CrossLinkSpectrumMatch > >& all_top_csms, const vector< PeptideIdentification >& peptide_ids, const RichPeakMap& spectra, String spec_xml_name, double cross_link_mass_light, vector< double > cross_link_mass_mono_link)
    {

      // TODO Missing Parameters/Information :

      ofstream xml_file;
      xml_file.open(out_file.c_str(), ios::trunc);
      // XML Header
      xml_file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
      xml_file << "<?xml-stylesheet type=\"text/xsl\" href=\"\"?>" << endl;

      // TODO write actual experiment data (from protein_ids[0] ?)
      xml_file << "<xquest_results xquest_version=\"xquest 2.1.1\" date=\"Fri Dec 18 12:28:23 2015\" author=\"Thomas Walzthoeni,Oliver Rinner\" homepage=\"http://proteomics.ethz.ch\" deffile=\"xquest.def\" tolerancemeasure_ms2=\"Da\" crosslinkername=\"DSS\" commonlossxcorrweigth=\"3\" poolisotopes=\"0\" xcorr_tolerance_window=\"0\" redundant_peps=\"0\" picktolerance=\"500\" monolinkmw=\"156.0786442,155.0964278\" search_maxcandidate_peps=\"250\" requiredmissed_cleavages=\"0\" picktolerance_measure=\"ppm\" database=\"/home/eugen/MSData/26S_testdataset/db/26Syeast.fasta\" fragmentresiduals=\"HASH(0x37998e0)\" xlinktypes=\"1111\" maxiontaghits=\"0\" AArequired=\"K\" miniontaghits=\"1\" cp_minpeaknumber=\"25\" cp_isotopediff=\"12.075321\" xkinkerID=\"DSS\" maxdigestlength=\"50\" y_ion=\"19.0183888\" cp_nhighest=\"100\" uselossionsformatching=\"0\" mindigestlength=\"5\" indexcharges_common=\"ARRAY(0x378c600)\" enzyme_num=\"1\" testionspick=\"intensity\" xlink_ms2tolerance=\"0.3\" waterloss=\"0\" database_dc=\"/home/eugen/MSData/26S_testdataset/db/26Syeast_decoy.fasta\" iontag_match_xlinkions=\"1\" averageMS2=\"0\" wTICweight=\"12.829\" minionsize=\"200\" x_ion=\"44.9976546\" commonxcorrweigth=\"10\" a_ion=\"-26.9870904\" drawlogscale=\"0\" fwd_ions=\"a|b|c\" experiment=\"BSAtest\" matchoddsweight=\"1.973\" nh3loss=\"0\" verbose=\"0\" enumeration_index_mode=\"smarthash\" maxionsize=\"2000\" outputpath=\"" << spec_xml_name << "\" search_intercrosslinks=\"1\" cp_tolerancemeasure=\"ppm\" reportnbesthits=\"5\" cp_dynamic_range=\"1000\" xlinkermw=\"138.0680796\" ms1tol_maxborder=\"10\" enumerate=\"0\" ionindexintprecision=\"10\" writetodiskaftern=\"100000000\" intsumweight=\"0.018\" xlinkxcorrweigth=\"10\" search_monolinks=\"1\" cp_peakratio=\"0.3\" rev_ions=\"x|y|z\" xcorrprecision=\"0.2\" RuntimeDecoys=\"1\" intprecision=\"10\" Iontag_charges_for_index=\"1\" Iontag_writeaftern=\"1200\" z_ion=\"2.9998388\" tryptic_termini=\"2\" cp_tolerance=\"400\" printpeptides=\"1\" ntestions=\"100\" b_ion=\"1.0078246\" normxcorr=\"1\" ioncharge_xlink=\"ARRAY(0x378c8d0)\" ionseries_array=\"ARRAY(0x378cd68)\" cp_scaleby=\"max\" printdigestpeps=\"0\" cp_tolerancexl=\"500\" search_intralinks=\"1\" minionintensity=\"1\" nvariable_mod=\"1\" missed_cleavages=\"2\" ntermxlinkable=\"0\" cp_threshold=\"1\" tolerancemeasure=\"ppm\" CID_match2ndisotope=\"1\" variable_mod=\"M,15.99491\" nocutatxlink=\"1\" minpepmr=\"550\" realintensities4xcorr=\"0\" printcandidatepeps=\"0\" xcorrdelay=\"5\" xcorrxweight=\"2.488\" minhits=\"1\" ms1tol_minborder=\"-10\" drawspectra=\"0\" cp_scaleintensity=\"1\" ms2tolerance=\"0.2\" maxpepmr=\"5500\" printtables=\"0\" ionseries=\"HASH(0x378cc00)\" usenprescores=\"100\" search_intracrosslinks=\"1\" Iontagmode=\"1\" ms1tolerance=\"10\" ioncharge_common=\"ARRAY(0x378c720)\" xcorrbweight=\"21.279\" c_ion=\"18.0343724\" copydb2resdir=\"1\" Hatom=\"1.007825032\" >" << endl;



      for (vector< vector< CrossLinkSpectrumMatch > >::const_iterator top_csms_spectrum = all_top_csms.begin(); top_csms_spectrum != all_top_csms.end(); ++top_csms_spectrum)
      {

        vector< CrossLinkSpectrumMatch > top_vector = (*top_csms_spectrum);

        if (!top_vector.empty())
        {
        // Spectrum Data, for each spectrum
        Size scan_index_light = top_vector[0].scan_index_light;
        Size scan_index_heavy = top_vector[0].scan_index_heavy;
//        cout << "Scan indices: " << scan_index_light << "\t" << scan_index_heavy << endl;
        const RichPeakSpectrum& spectrum_light = spectra[scan_index_light];
        double precursor_charge = spectrum_light.getPrecursors()[0].getCharge();

        double precursor_mz = spectrum_light.getPrecursors()[0].getMZ();
        double precursor_rt = spectrum_light.getRT();
        double precursor_mass = precursor_mz * static_cast<double>(precursor_charge) - static_cast<double>(precursor_charge) * Constants::PROTON_MASS_U;

//        double precursor_charge_heavy = spectra[scan_index_heavy].getPrecursors()[0].getCharge();
        double precursor_mz_heavy = spectra[scan_index_heavy].getPrecursors()[0].getMZ();
        double precursor_rt_heavy = spectra[scan_index_heavy].getRT();

        // print information about new peak to file (starts with <spectrum_search..., ends with </spectrum_search>
        // TODO what to do with useless information? leave =0 or delete?
        String spectrum_light_name = String("spectrumlight") + scan_index_light;
        String spectrum_heavy_name = String("spectrumheavy") + scan_index_heavy;
        String spectrum_name = spectrum_light_name + String("_") + spectrum_heavy_name;
        String rt_scans = String(precursor_rt) + ":" + String(precursor_rt_heavy);
        String mz_scans = String(precursor_mz) + ":" + String(precursor_mz_heavy);

        // Mean ion intensity (light spectrum, TODO add heavy spectrum?)
        double mean_intensity= 0;
        for (SignedSize j = 0; j < static_cast<SignedSize>(spectrum_light.size()); ++j) mean_intensity += spectrum_light[j].getIntensity();
        for (SignedSize j = 0; j < static_cast<SignedSize>(spectra[scan_index_heavy].size()); ++j) mean_intensity += spectra[scan_index_heavy][j].getIntensity();
        mean_intensity = mean_intensity / (spectrum_light.size() + spectra[scan_index_heavy].size());

        xml_file << "<spectrum_search spectrum=\"" << spectrum_name << "\" mean_ionintensity=\"" << mean_intensity << "\" ionintensity_stdev=\"" << "TODO" << "\" addedMass=\"" << "TODO" << "\" iontag_ncandidates=\"" << "TODO"
          << "\"  apriori_pmatch_common=\"" << "TODO" << "\" apriori_pmatch_xlink=\"" << "TODO" << "\" ncommonions=\"" << "TODO" << "\" nxlinkions=\"" << "TODO" << "\" mz_precursor=\"" << precursor_mz
          << "\" scantype=\"" << "light_heavy" << "\" charge_precursor=\"" << precursor_charge << "\" Mr_precursor=\"" << precursor_mass <<  "\" rtsecscans=\"" << rt_scans << "\" mzscans=\"" << mz_scans << "\" >" << endl;


        for (vector< CrossLinkSpectrumMatch>::const_iterator top_csm = top_csms_spectrum->begin(); top_csm != top_csms_spectrum->end(); ++top_csm)
        {
          String xltype = "monolink";
          String structure = top_csm->cross_link.alpha.toUnmodifiedString();
          String topology = String("K") + (top_csm->cross_link.cross_link_position.first+1);

           // TODO track or otherwise find out, which kind of mono-link it was (if there are several possibilities for the weigths)
           double weight = top_csm->cross_link.alpha.getMonoWeight() + cross_link_mass_mono_link[1];

           bool is_monolink = (top_csm->cross_link.cross_link_position.second == -1);
           int alpha_pos = top_csm->cross_link.cross_link_position.first + 1;
           int beta_pos = top_csm->cross_link.cross_link_position.second + 1;

           if (!is_monolink)
           {
            xltype = "xlink";
            structure =top_csm->cross_link.alpha.toUnmodifiedString() + "-" + top_csm->cross_link.beta.toUnmodifiedString();
            topology = String("a") +  alpha_pos + String("-b") + beta_pos;
            weight = top_csm->cross_link.alpha.getMonoWeight() + top_csm->cross_link.beta.getMonoWeight() + cross_link_mass_light;
           }
           String id = structure + "-" + topology;

           // Error calculation
           double cl_mz = (weight + static_cast<double>(precursor_charge) * Constants::PROTON_MASS_U) / static_cast<double>(precursor_charge);
           double error = abs(cl_mz - precursor_mz);
           double rel_error = error / precursor_mz / 1e-6;
//           cout << "Relative error (ppm): " << rel_error << endl;

           // Protein Accessions
           String prot_alpha;

//           cout << "Size and index: " << peptide_ids.size() << "\t | \t" << top_csm->peptide_id_index << endl;

           PeptideIdentification pep_id = peptide_ids[top_csm->peptide_id_index];
//           PeptideIdentification *pep_id_test = top_csm->peptide_id;

//           cout << "MZ, pointer, current adress: " << pep_id.getMZ() << "\t | \t" << pep_id_test << "\t | \t" << &peptide_ids[top_csm->peptide_id_index] << endl;

           vector< PeptideHit > pep_hits = pep_id.getHits();
//           vector< PeptideHit > pep_hits_test = pep_id_test->getHits();

//           cout << "n of hits, 1. seq: " << pep_hits.size() << "\t | \t" << pep_hits_test.size() << "\t | \t" << pep_hits[0].getSequence().toString() << endl;

           prot_alpha = pep_hits[0].getPeptideEvidences()[0].getProteinAccession();

           String prot_beta = "-";
//           const vector< PeptideEvidence > &evidences_beta =ph_beta.getPeptideEvidences();
           if (pep_hits.size() > 1)
           {
              prot_beta= pep_hits[1].getPeptideEvidences()[0].getProteinAccession();
//              cout << "Proteins: " << pep_hits[0].getPeptideEvidences()[0].getProteinAccession() << "\t" << pep_hits_test[1].getPeptideEvidences()[0].getProteinAccession() << endl;
           }

           // Hit Data, for each cross-link to Spectrum Hit (e.g. top 5 per spectrum)
           xml_file << "<search_hit search_hit_rank=\"" <<top_csm->rank << "\" id=\"" << id << "\" type=\"" << xltype << "\" structure=\"" << structure << "\" seq1=\"" << top_csm->cross_link.alpha.toUnmodifiedString() << "\" seq2=\"" << top_csm->cross_link.beta.toUnmodifiedString()
                << "\" prot1=\"" << prot_alpha << "\" prot2=\"" << prot_beta << "\" topology=\"" << topology << "\" xlinkposition=\"" << (top_csm->cross_link.cross_link_position.first+1) << "," << (top_csm->cross_link.cross_link_position.second+1)
                << "\" Mr=\"" << weight << "\" mz=\"" << cl_mz << "\" charge=\"" << precursor_charge << "\" xlinkermass=\"" << cross_link_mass_light << "\" measured_mass=\"" << precursor_mass << "\" error=\"" << error
                << "\" error_rel=\"" << rel_error << "\" xlinkions_matched=\"" << (top_csm->matched_xlink_alpha + top_csm->matched_xlink_beta) << "\" backboneions_matched=\"" << (top_csm->matched_common_alpha + top_csm->matched_common_beta)
                << "\" weighted_matchodds_mean=\"" << "TODO" << "\" weighted_matchodds_sum=\"" << "TODO" << "\" match_error_mean=\"" << "TODO" << "\" match_error_stdev=\"" << "TODO" << "\" xcorrx=\"" << top_csm->xcorrx_max << "\" xcorrb=\"" << top_csm->xcorrc_max << "\" match_odds=\"" <<top_csm->match_odds << "\" prescore=\"" << top_csm->pre_score
                << "\" prescore_alpha=\"" << "TODO" << "\" prescore_beta=\"" << "TODO" << "\" match_odds_alphacommon=\"" << "TODO" << "\" match_odds_betacommon=\"" << "TODO" << "\" match_odds_alphaxlink=\"" << "TODO"
                << "\" match_odds_betaxlink=\"" << "TODO" << "\" num_of_matched_ions_alpha=\"" << (top_csm->matched_common_alpha + top_csm->matched_xlink_alpha) << "\" num_of_matched_ions_beta=\"" << (top_csm->matched_common_beta + top_csm->matched_xlink_beta) << "\" num_of_matched_common_ions_alpha=\"" << top_csm->matched_common_alpha
                << "\" num_of_matched_common_ions_beta=\"" << top_csm->matched_common_beta << "\" num_of_matched_xlink_ions_alpha=\"" << top_csm->matched_xlink_alpha << "\" num_of_matched_xlink_ions_beta=\"" << top_csm->matched_xlink_beta << "\" xcorrall=\"" << "TODO" << "\" TIC=\"" << top_csm->percTIC
                << "\" TIC_alpha=\"" << "TODO" << "\" TIC_beta=\"" << "TODO" << "\" wTIC=\"" << top_csm->wTIC << "\" intsum=\"" << top_csm->int_sum * 100 << "\" apriori_match_probs=\"" << "TODO" << "\" apriori_match_probs_log=\"" << "TODO"
                << "\" series_score_mean=\"" << "TODO" << "\" annotated_spec=\"" << "" << "\" score=\"" << top_csm->score << "\" >" << endl;
           xml_file << "</search_hit>" << endl;
          }
          // Closing tag for Spectrum
          xml_file << "</spectrum_search>" << endl;


        }

      }

       // Closing tag for results (end of file)
      xml_file << "</xquest_results>" << endl;
      xml_file.close();

      return;
    }

    void writeXQuestXMLSpec(String spec_xml_filename, const RichPeakMap& spectra, const PreprocessedPairSpectra_& preprocessed_pair_spectra, const vector< pair<Size, Size> >& spectrum_pairs, const vector< vector< CrossLinkSpectrumMatch > >& all_top_csms)
    {

      // XML Header
      ofstream spec_xml_file;
      spec_xml_file.open(spec_xml_filename.c_str(), ios::trunc); // ios::app = append to file, ios::trunc = overwrites file
      spec_xml_file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?><xquest_spectra compare_peaks_version=\"3.4\" date=\"Tue Nov 24 12:41:18 2015\" author=\"Thomas Walzthoeni,Oliver Rinner\" homepage=\"http://proteomics.ethz.ch\" resultdir=\"aleitner_M1012_004_matched\" deffile=\"xquest.def\" >" << endl;

      for (Size i = 0; i < spectrum_pairs.size(); ++i)
      {
        if (!all_top_csms[i].empty())
        {
          Size scan_index_light = spectrum_pairs[i].first;
          Size scan_index_heavy = spectrum_pairs[i].second;
          String spectrum_light_name = String("spectrumlight") + scan_index_light;
          String spectrum_heavy_name = String("spectrumheavy") + scan_index_heavy;
          String spectrum_name = spectrum_light_name + String("_") + spectrum_heavy_name;

          // 4 Spectra resulting from a light/heavy spectra pair.  Write for each spectrum, that is written to xquest.xml (should be all considered pairs, or better only those with at least one sensible Hit, meaning a score was computed)
          spec_xml_file << "<spectrum filename=\"" << spectrum_light_name << ".dta" << "\" type=\"light\">" << endl;
          spec_xml_file << TOPPxQuest::getxQuestBase64EncodedSpectrum_(spectra[scan_index_light], String(""));
          spec_xml_file << "</spectrum>" << endl;

          spec_xml_file << "<spectrum filename=\"" << spectrum_heavy_name << ".dta" << "\" type=\"heavy\">" << endl;
          spec_xml_file << TOPPxQuest::getxQuestBase64EncodedSpectrum_(spectra[scan_index_heavy], String(""));
          spec_xml_file << "</spectrum>" << endl;

          String spectrum_common_name = spectrum_name + String("_common.txt");
          spec_xml_file << "<spectrum filename=\"" << spectrum_common_name << "\" type=\"common\">" << endl;
          spec_xml_file << TOPPxQuest::getxQuestBase64EncodedSpectrum_(preprocessed_pair_spectra.spectra_common_peaks[i], spectrum_light_name + ".dta," + spectrum_heavy_name + ".dta");
          spec_xml_file << "</spectrum>" << endl;

          String spectrum_xlink_name = spectrum_name + String("_xlinker.txt");
          spec_xml_file << "<spectrum filename=\"" << spectrum_xlink_name << "\" type=\"xlinker\">" << endl;
          spec_xml_file <<TOPPxQuest::getxQuestBase64EncodedSpectrum_(preprocessed_pair_spectra.spectra_xlink_peaks[i], spectrum_light_name + ".dta," + spectrum_heavy_name + ".dta");
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


    /*
    // TEMPORARY TESTCODE ###########
    cout.precision(20);
    cout << "cumulBinom = " << cumulativeBinomial(33, 16, 0.00777026) << endl;

    long double p_cumul = 0.0;
    Size n = 33;
    Size k = 16;
    double prob = 0.00777026;
    for (Size j = 0; j < k; j++)
    {
      double coeff = boost::math::binomial_coefficient<double>((unsigned int)n, (unsigned int)j);
      p_cumul += coeff * pow(prob,  j) * pow((1-prob), (n-j));
      cout << "p_cumul = " << p_cumul << "\t pow(prob,  j) = " << pow(prob,  j)  << "\t pow((1-prob), (n-j)) = " << pow((1-prob), (n-j)) << endl;
    }

    cout << "Smallest double: " << numeric_limits<double>::min() << "Double closest to 1: " <<  nexttoward(1.0, 0.0) <<endl;

    return EXECUTION_OK;
    //##########################
    */


    const string in_mzml(getStringOption_("in"));
    const string in_fasta(getStringOption_("database"));
    const string in_decoy_fasta(getStringOption_("decoy_database"));
    const string in_consensus(getStringOption_("consensus"));
    const string out_idxml(getStringOption_("out_idXML"));
    const string out_xquest = getStringOption_("out");

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
    double cross_link_mass_light = getDoubleOption_("cross_linker:mass_light");
    double cross_link_mass_iso_shift = getDoubleOption_("cross_linker:mass_iso_shift");
    DoubleList cross_link_mass_mono_link = getDoubleList_("cross_linker:mass_mono_link");

    StringList fixedModNames = getStringList_("modifications:fixed");
    set<String> fixed_unique(fixedModNames.begin(), fixedModNames.end());

    Size peptide_min_size = getIntOption_("peptide:min_size");

    bool ion_index_mode = (getStringOption_("algorithm:candidate_search") == "index");

    if (fixed_unique.size() != fixedModNames.size())
    {
      cout << "duplicate fixed modification provided." << endl;
      return ILLEGAL_PARAMETERS;
    }

    StringList varModNames = getStringList_("modifications:variable");
    set<String> var_unique(varModNames.begin(), varModNames.end());
    if (var_unique.size() != varModNames.size())
    {
      cout << "duplicate variable modification provided." << endl;
      return ILLEGAL_PARAMETERS;
    }

    vector<ResidueModification> fixed_modifications = getModifications_(fixedModNames);
    vector<ResidueModification> variable_modifications = getModifications_(varModNames);
    Size max_variable_mods_per_peptide = getIntOption_("modifications:variable_max_per_peptide");
    
    // load MS2 map
    PeakMap unprocessed_spectra;
    RichPeakMap spectra;
    MzMLFile f;
    f.setLogType(log_type_);

    PeakFileOptions options;
    options.clearMSLevels();
    options.addMSLevel(2);
    f.getOptions() = options;
    f.load(in_mzml, unprocessed_spectra);
    unprocessed_spectra.sortSpectra(true);

    // filter noise
    progresslogger.startProgress(0, 1, "Filtering spectra...");
    preprocessSpectra_(unprocessed_spectra, spectra);
    progresslogger.endProgress();

    // load linked features
    ConsensusMap cfeatures;
    ConsensusXMLFile cf;
    cf.load(in_consensus, cfeatures); 

    // load fasta database    
    progresslogger.startProgress(0, 1, "Load database from FASTA file...");
    FASTAFile fastaFile;
    vector<FASTAFile::FASTAEntry> fasta_db;
    fastaFile.load(in_fasta, fasta_db);

    progresslogger.endProgress();
    
    const Size missed_cleavages = getIntOption_("peptide:missed_cleavages");
    EnzymaticDigestion digestor;
    digestor.setEnzyme(getStringOption_("peptide:enzyme"));
    digestor.setMissedCleavages(missed_cleavages);

//    progresslogger.startProgress(0, static_cast<Size>(fasta_db.end() - fasta_db.begin()), "Scoring peptide models against spectra...");

    // lookup for processed peptides. must be defined outside of omp section and synchronized
    multimap<StringView, AASequence> processed_peptides;
    
    // set minimum size of peptide after digestion
    Size min_peptide_length = getIntOption_("peptide:min_size");

    // build multimap of precursor mass to scan index
    multimap<double, Size> multimap_mass_2_scan_index;

    vector<PeptideIdentification> pseudo_ids; // used to map precursor positions to consensus features
    for (RichPeakMap::ConstIterator s_it = spectra.begin(); s_it != spectra.end(); ++s_it)
    {
      int scan_index = s_it - spectra.begin();
      vector<Precursor> precursor = s_it->getPrecursors();
    
      // there should only one precursor and MS2 should contain at least a few peaks to be considered (e.g. at least for every AA in the peptide)
      if (precursor.size() == 1 && s_it->size() >= peptide_min_size)
      {
        int precursor_charge = precursor[0].getCharge();
        if (precursor_charge < min_precursor_charge || precursor_charge > max_precursor_charge)
        {
          continue;
        }
    
        double precursor_mz = precursor[0].getMZ();
        double precursor_mass = static_cast<double>(precursor_charge) * precursor_mz - static_cast<double>(precursor_charge) * Constants::PROTON_MASS_U;
    
        multimap_mass_2_scan_index.insert(make_pair(precursor_mass, scan_index));
        PeptideIdentification temp_pi;
        temp_pi.setRT(s_it->getRT());
        temp_pi.setMZ(precursor_mz);
        temp_pi.setMetaValue("scan_index", scan_index);
        vector<PeptideHit> temp_hits;
        PeptideHit temp_ph;
        temp_ph.setCharge(precursor_charge);
        temp_hits.push_back(temp_ph);
        temp_pi.setHits(temp_hits);
        pseudo_ids.push_back(temp_pi);
      }
    }

    IDMapper idmapper;
    Param p = idmapper.getParameters();
    p.setValue("rt_tolerance", 30.0);
    p.setValue("mz_tolerance", precursor_mass_tolerance);
    String mz_measure = precursor_mass_tolerance_unit_ppm ? "ppm" : "Da";
    p.setValue("mz_measure", mz_measure);
    p.setValue("mz_reference", "precursor");
    p.setValue("ignore_charge", "false");
    idmapper.setParameters(p);

    vector< pair<Size, Size> > spectrum_pairs;
    {
      vector<ProteinIdentification> protein_ids;
      // protein identifications (leave as is...)
      protein_ids = vector<ProteinIdentification>(1);
      protein_ids[0].setDateTime(DateTime::now());
      protein_ids[0].setSearchEngineVersion(VersionInfo::getVersion());
      protein_ids[0].setMetaValue("SpectrumIdentificationProtocol", DataValue("MS:1002494")); // cross-linking search = MS:1002494

      idmapper.annotate(cfeatures, pseudo_ids, protein_ids, true, true);

    // maps the index of a light precursor peptide to its corresponding heavier partner
//    map<Size, Size> map_light_to_heavy;
      for (ConsensusMap::const_iterator cit = cfeatures.begin(); cit != cfeatures.end(); ++cit)
      {
//      if (cit->getFeatures().size() > 2 || cit->getPeptideIdentifications().size() > 2)
//      {
//        cout << "Multiple PeptideIDs, Features: " << cit->getFeatures().size() << "\t PeptideIDs: " << cit->getPeptideIdentifications().size() << endl;
//      }

      // if x == light && y == heavy, then make pair
        if (cit->getFeatures().size() == 2 && cit->getPeptideIdentifications().size() >= 2)
        {
          for (Size x = 0; x < cit->getPeptideIdentifications().size(); ++x)
          {
//          cout << "MetaValue map index: " <<  cit->getPeptideIdentifications()[x].getMetaValue("map index") << endl;
//          cout << "Scan index 1: " << cit->getPeptideIdentifications()[x].getMetaValue("scan_index") << endl;
            if (static_cast<int>(cit->getPeptideIdentifications()[x].getMetaValue("map index")) == 0)
            {
              for (Size y = 0; y < cit->getPeptideIdentifications().size(); ++y)  {
//            cout << "Scan index 2: " << cit->getPeptideIdentifications()[y].getMetaValue("scan_index") << endl;
                if (static_cast<int>(cit->getPeptideIdentifications()[y].getMetaValue("map index")) == 1)
                {
                  const PeptideIdentification& pi_0 = cit->getPeptideIdentifications()[x];
                  const PeptideIdentification& pi_1 = cit->getPeptideIdentifications()[y];
//                map_light_to_heavy[pi_0.getMetaValue("scan_index")] = pi_1.getMetaValue("scan_index");
                  spectrum_pairs.push_back(make_pair(pi_0.getMetaValue("scan_index"), pi_1.getMetaValue("scan_index")));

                }
              }
            }
          }
        }
      }
    }

    //cout << "Number of MS2 pairs connceted by consensus feature: " << map_light_to_heavy.size() << endl;

    // create common peak / shifted peak spectra for all pairs
    progresslogger.startProgress(0, 1, "Preprocessing Spectra Pairs...");
//    PreprocessedPairSpectra_ preprocessed_pair_spectra = preprocessPairs_(spectra, map_light_to_heavy, cross_link_mass_light, cross_link_mass_heavy);
    PreprocessedPairSpectra_ preprocessed_pair_spectra = preprocessPairs_(spectra, spectrum_pairs, cross_link_mass_light, cross_link_mass_iso_shift);
    progresslogger.endProgress();
 
    Size count_proteins = 0;
    Size count_peptides = 0;


    // one identification run
    vector<ProteinIdentification> protein_ids(1);
    protein_ids[0].setDateTime(DateTime::now());
    protein_ids[0].setSearchEngine("OpenMSxQuest");
    protein_ids[0].setSearchEngineVersion(VersionInfo::getVersion());
    protein_ids[0].setPrimaryMSRunPath(spectra.getPrimaryMSRunPath());

    vector<PeptideIdentification> peptide_ids;

#ifdef _OPENMP
#pragma omp parallel for
#endif
    // digest and filter database
    for (SignedSize fasta_index = 0; fasta_index < static_cast<SignedSize>(fasta_db.size()); ++fasta_index)
    {
#ifdef _OPENMP
#pragma omp atomic
#endif
      ++count_proteins;

      IF_MASTERTHREAD
      {
        progresslogger.setProgress(static_cast<SignedSize>(fasta_index) * NUMBER_OF_THREADS);
      }

      // store vector of substrings pointing in fasta database (bounded by pairs of begin, end iterators)    
      vector<StringView> current_digest;
      digestor.digestUnmodifiedString(fasta_db[fasta_index].sequence, current_digest, min_peptide_length);

      for (vector<StringView>::iterator cit = current_digest.begin(); cit != current_digest.end(); ++cit)
      {
        // skip peptides with invalid AAs
        if (cit->getString().has('B') || cit->getString().has('O') || cit->getString().has('U') || cit->getString().has('X') || cit->getString().has('Z')) continue;

        // skip if no K
        if (!cit->getString().has('K')) continue; 

        bool already_processed = false;
#ifdef _OPENMP
#pragma omp critical (processed_peptides_access)
#endif
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
        if (cit->getString().find('K') >= cit->getString().size()-1)
        {
          continue;
        }
            
#ifdef _OPENMP
#pragma omp atomic
#endif
        ++count_peptides;

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

#ifdef _OPENMP
#pragma omp critical (processed_peptides_access)
#endif
          {
            processed_peptides.insert(pair<StringView, AASequence>(*cit, candidate));
          }
        }
      }
    }


    // digest and filter decoy database
    if (!in_decoy_fasta.empty())
    {
      vector<FASTAFile::FASTAEntry> fasta_decoys;
      fastaFile.load(in_decoy_fasta, fasta_decoys);

      for (SignedSize fasta_index = 0; fasta_index < static_cast<SignedSize>(fasta_decoys.size()); ++fasta_index)
      {
  #ifdef _OPENMP
  #pragma omp atomic
  #endif
        ++count_proteins;

        IF_MASTERTHREAD
        {
          progresslogger.setProgress(static_cast<SignedSize>(fasta_index) * NUMBER_OF_THREADS);
        }

        // store vector of substrings pointing in fasta database (bounded by pairs of begin, end iterators)
        vector<StringView> current_digest;
        digestor.digestUnmodifiedString(fasta_decoys[fasta_index].sequence, current_digest, min_peptide_length);

        for (vector<StringView>::iterator cit = current_digest.begin(); cit != current_digest.end(); ++cit)
        {
          // skip peptides with invalid AAs
          if (cit->getString().has('B') || cit->getString().has('O') || cit->getString().has('U') || cit->getString().has('X') || cit->getString().has('Z')) continue;

          // skip if no K
          if (!cit->getString().has('K')) continue;

          bool already_processed = false;
  #ifdef _OPENMP
  #pragma omp critical (processed_peptides_access)
  #endif
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
          if (cit->getString().find('K') >= cit->getString().size()-1)
          {
            continue;
          }

  #ifdef _OPENMP
  #pragma omp atomic
  #endif
          ++count_peptides;

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

  #ifdef _OPENMP
  #pragma omp critical (processed_peptides_access)
  #endif
            {
              processed_peptides.insert(pair<StringView, AASequence>(*cit, candidate));
              fasta_db.reserve(fasta_db.size() + fasta_decoys.size());
              fasta_db.insert(fasta_db.end(), fasta_decoys.begin(), fasta_decoys.end());

              // TODO doeas this actually save space or something? or is the object removed automatically, since it is not used after this
              fasta_decoys.clear();
            }
          }
        }
      }
    }

    // Initialize output to file
//    String out_file = getStringOption_("out");
//    String spec_xml_name = getStringOption_("in").prefix(getStringOption_("in").size()-7) + "_matched";
//    String spec_xml_filename = spec_xml_name + ".spec.xml";
//    ofstream xml_file;
//    xml_file.open(out_file.c_str(), ios::trunc);
//    xml_file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
//    xml_file << "<?xml-stylesheet type=\"text/xsl\" href=\"\"?>" << endl;
//    xml_file << "<xquest_results xquest_version=\"xquest 2.1.1\" date=\"Fri Dec 18 12:28:23 2015\" author=\"Thomas Walzthoeni,Oliver Rinner\" homepage=\"http://proteomics.ethz.ch\" deffile=\"xquest.def\" tolerancemeasure_ms2=\"Da\" crosslinkername=\"DSS\" commonlossxcorrweigth=\"3\" poolisotopes=\"0\" xcorr_tolerance_window=\"0\" redundant_peps=\"0\" picktolerance=\"500\" monolinkmw=\"156.0786442,155.0964278\" search_maxcandidate_peps=\"250\" requiredmissed_cleavages=\"0\" picktolerance_measure=\"ppm\" database=\"/home/eugen/MSData/26S_testdataset/db/26Syeast.fasta\" fragmentresiduals=\"HASH(0x37998e0)\" xlinktypes=\"1111\" maxiontaghits=\"0\" AArequired=\"K\" miniontaghits=\"1\" cp_minpeaknumber=\"25\" cp_isotopediff=\"12.075321\" xkinkerID=\"DSS\" maxdigestlength=\"50\" y_ion=\"19.0183888\" cp_nhighest=\"100\" uselossionsformatching=\"0\" mindigestlength=\"5\" indexcharges_common=\"ARRAY(0x378c600)\" enzyme_num=\"1\" testionspick=\"intensity\" xlink_ms2tolerance=\"0.3\" waterloss=\"0\" database_dc=\"/home/eugen/MSData/26S_testdataset/db/26Syeast_decoy.fasta\" iontag_match_xlinkions=\"1\" averageMS2=\"0\" wTICweight=\"12.829\" minionsize=\"200\" x_ion=\"44.9976546\" commonxcorrweigth=\"10\" a_ion=\"-26.9870904\" drawlogscale=\"0\" fwd_ions=\"a|b|c\" experiment=\"BSAtest\" matchoddsweight=\"1.973\" nh3loss=\"0\" verbose=\"0\" enumeration_index_mode=\"smarthash\" maxionsize=\"2000\" outputpath=\"" << spec_xml_name << "\" search_intercrosslinks=\"1\" cp_tolerancemeasure=\"ppm\" reportnbesthits=\"5\" cp_dynamic_range=\"1000\" xlinkermw=\"138.0680796\" ms1tol_maxborder=\"10\" enumerate=\"0\" ionindexintprecision=\"10\" writetodiskaftern=\"100000000\" intsumweight=\"0.018\" xlinkxcorrweigth=\"10\" search_monolinks=\"1\" cp_peakratio=\"0.3\" rev_ions=\"x|y|z\" xcorrprecision=\"0.2\" RuntimeDecoys=\"1\" intprecision=\"10\" Iontag_charges_for_index=\"1\" Iontag_writeaftern=\"1200\" z_ion=\"2.9998388\" tryptic_termini=\"2\" cp_tolerance=\"400\" printpeptides=\"1\" ntestions=\"100\" b_ion=\"1.0078246\" normxcorr=\"1\" ioncharge_xlink=\"ARRAY(0x378c8d0)\" ionseries_array=\"ARRAY(0x378cd68)\" cp_scaleby=\"max\" printdigestpeps=\"0\" cp_tolerancexl=\"500\" search_intralinks=\"1\" minionintensity=\"1\" nvariable_mod=\"1\" missed_cleavages=\"2\" ntermxlinkable=\"0\" cp_threshold=\"1\" tolerancemeasure=\"ppm\" CID_match2ndisotope=\"1\" variable_mod=\"M,15.99491\" nocutatxlink=\"1\" minpepmr=\"550\" realintensities4xcorr=\"0\" printcandidatepeps=\"0\" xcorrdelay=\"5\" xcorrxweight=\"2.488\" minhits=\"1\" ms1tol_minborder=\"-10\" drawspectra=\"0\" cp_scaleintensity=\"1\" ms2tolerance=\"0.2\" maxpepmr=\"5500\" printtables=\"0\" ionseries=\"HASH(0x378cc00)\" usenprescores=\"100\" search_intracrosslinks=\"1\" Iontagmode=\"1\" ms1tolerance=\"10\" ioncharge_common=\"ARRAY(0x378c720)\" xcorrbweight=\"21.279\" c_ion=\"18.0343724\" copydb2resdir=\"1\" Hatom=\"1.007825032\" >" << endl;

//    //String spec_xml_filename = getStringOption_("in") + "spec.xml";
//    ofstream spec_xml_file;
//    spec_xml_file.open(spec_xml_filename.c_str(), ios::trunc); // ios::app = append to file, ios::trunc = overwrites file
//    spec_xml_file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?><xquest_spectra compare_peaks_version=\"3.4\" date=\"Tue Nov 24 12:41:18 2015\" author=\"Thomas Walzthoeni,Oliver Rinner\" homepage=\"http://proteomics.ethz.ch\" resultdir=\"aleitner_M1012_004_matched\" deffile=\"xquest.def\" >" << endl;



    // create spectrum generator
    TheoreticalSpectrumGenerator spectrum_generator;
    TheoreticalSpectrumGeneratorXLinks specGen;

    // TODO constant binsize for HashGrid and aPrioriProb computation
    double tolerance_binsize = 0.2;

    cout << "Peptide candidates: " << processed_peptides.size() << endl;

    //TODO refactor, so that only the used mode is initialized and the pre-scoring code only appears once
    // Initialize enumeration mode
    multimap<double, pair<const AASequence*, const AASequence*> > enumerated_cross_link_masses;

    //TODO remove, adapt to ppm
    HashGrid1D hg(0.0, 20000.0, tolerance_binsize);

    if (!ion_index_mode)
    {
      progresslogger.startProgress(0, 1, "Enumerating cross-links...");
      enumerated_cross_link_masses = enumerateCrossLinksAndMasses_(processed_peptides, cross_link_mass_light, cross_link_mass_mono_link);
      progresslogger.endProgress();
      cout << "Enumerated cross-links: " << enumerated_cross_link_masses.size() << endl;
    }
    else
    {
      cout << "Adding peaks to hash map ...";
      for (map<StringView, AASequence>::iterator a = processed_peptides.begin(); a != processed_peptides.end(); ++a)
      {
        //create theoretical spectrum
        MSSpectrum<RichPeak1D> theo_spectrum = MSSpectrum<RichPeak1D>();
        // cout << a->second.toString() << endl;
        AASequence * seq = &(a->second);
        // generate common ions
        specGen.getCommonIonSpectrum(theo_spectrum, *seq, 3); // TODO check which charge and which ion series are used for ion index
      
        //sort by mz (is done in getCommonIonSpectrum)
//        theo_spectrum.sortByPosition();

        for (Size i = 0; i != theo_spectrum.size(); ++i)
        {
          hg.insert(theo_spectrum[i].getMZ(), seq);
        }
      }
      cout << " finished."  << endl;
    }

    // TODO collect peptide Hits and stuff to write out to mzIdentML
    ProteinIdentification prot_id_run;


    // TODO test variable, can be removed
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

//    for (SignedSize scan_index = 0; scan_index < static_cast<SignedSize>(spectra.size()); ++scan_index)
    for (Size pair_index = 0; pair_index < spectrum_pairs.size(); ++pair_index)
    {
      // assume that current MS2 corresponds to a light spectrum
      Size scan_index = spectrum_pairs[pair_index].first;
      Size scan_index_heavy = spectrum_pairs[pair_index].second;
      cout << "Scan indices: " << scan_index << "\t" << scan_index_heavy << endl;
      const RichPeakSpectrum& spectrum_light = spectra[scan_index];
      const double precursor_charge = spectrum_light.getPrecursors()[0].getCharge();
      const double precursor_mz = spectrum_light.getPrecursors()[0].getMZ();
      const double precursor_mass = precursor_mz * static_cast<double>(precursor_charge) - static_cast<double>(precursor_charge) * Constants::PROTON_MASS_U;


       // print information about new peak to file (starts with <spectrum_search..., ends with </spectrum_search>
       // TODO what to do with useless information? leave =0 or delete?
       String spectrum_light_name = String("spectrumlight") + scan_index;
       String spectrum_heavy_name = String("spectrumheavy") + scan_index_heavy;
       String spectrum_name = spectrum_light_name + String("_") + spectrum_heavy_name;

       // Mean ion intensity (light spectrum, TODO add heavy spectrum?)
       double mean_intensity= 0;
       for (SignedSize j = 0; j < static_cast<SignedSize>(spectrum_light.size()); ++j) mean_intensity += spectrum_light[j].getIntensity();
       mean_intensity = mean_intensity / spectrum_light.size();

//       xml_file << "<spectrum_search spectrum=\"" << spectrum_name << "\" mean_ionintensity=\"" << mean_intensity << "\" ionintensity_stdev=\"" << "TODO" << "\" addedMass=\"" << "TODO" << "\" iontag_ncandidates=\"" << "TODO"
//         << "\"  apriori_pmatch_common=\"" << "TODO" << "\" apriori_pmatch_xlink=\"" << "TODO" << "\" ncommonions=\"" << "TODO" << "\" nxlinkions=\"" << "TODO" << "\" mz_precursor=\"" << precursor_mz
//         << "\" scantype=\"" << "light_heavy" << "\" charge_precursor=\"" << precursor_charge << "\" Mr_precursor=\"" << precursor_mass <<  "\" rtsecscans=\"" << "TODO" << "\" mzscans=\"" << "0:1" << "\" >" << endl;

   

      // map light to heavy
//      map<Size, Size>::const_iterator scan_index_light_it = map_light_to_heavy.find(scan_index);
      //cout << "? " << scan_index << "\t" << preprocessed_pair_spectra.spectra_common_peaks[scan_index] << endl << spectrum_light[70].getMZ() << endl;
//      if (scan_index_light_it != map_light_to_heavy.end())

      //cout << " << scan_index << "\t" << " (mz, charge, mass) " << precursor_mz << "," << precursor_charge << "," << precursor_mass <<  endl;

//      // Write light, heavy, common, xlink spectra at position peak_index[i] to input-file.spec.xml for spectra visualization in the xQuest results manager
//      //String spectrum_light_name = String("spectrum_light_") + scan_index;
//      spec_xml_file << "<spectrum filename=\"" << spectrum_light_name << ".dta" << "\" type=\"light\">" << endl;
////      RichPeakSpectrum rich_light_spectrum = makeRichPeakSpectrum(spectrum_light, false);
////      spec_xml_file << TOPPxQuest::getxQuestBase64EncodedSpectrum_(rich_light_spectrum, String(""));
//      spec_xml_file << TOPPxQuest::getxQuestBase64EncodedSpectrum_(spectrum_light, String(""));
//      spec_xml_file << "</spectrum>" << endl;

//      //String spectrum_heavy_name = String("spectrum_heavy_") + scan_index;
//      spec_xml_file << "<spectrum filename=\"" << spectrum_heavy_name << ".dta" << "\" type=\"heavy\">" << endl;
//      //map<Size, Size>::const_iterator scan_index_light_it = map_light_to_heavy.find(scan_index);
////      const Size scan_index_heavy = scan_index_light_it->second;
////      RichPeakSpectrum rich_heavy_spectrum = makeRichPeakSpectrum(spectra[scan_index_heavy], false);
////      spec_xml_file << TOPPxQuest::getxQuestBase64EncodedSpectrum_(rich_heavy_spectrum, String(""));
//      spec_xml_file << TOPPxQuest::getxQuestBase64EncodedSpectrum_(spectra[scan_index_heavy], String(""));
//      spec_xml_file << "</spectrum>" << endl;

//      //String spectrum_name = spectrum_light_name + String("_") + spectrum_heavy_name;
//      String spectrum_common_name = spectrum_name + String("_common.txt");
//      spec_xml_file << "<spectrum filename=\"" << spectrum_common_name << "\" type=\"common\">" << endl;
////      RichPeakSpectrum rich_common_spectrum = makeRichPeakSpectrum(preprocessed_pair_spectra.spectra_common_peaks[scan_index], true);
////      spec_xml_file << TOPPxQuest::getxQuestBase64EncodedSpectrum_(rich_common_spectrum, spectrum_light_name + ".dta," + spectrum_heavy_name + ".dta");
//      spec_xml_file << TOPPxQuest::getxQuestBase64EncodedSpectrum_(preprocessed_pair_spectra.spectra_common_peaks[pair_index], spectrum_light_name + ".dta," + spectrum_heavy_name + ".dta");
//      spec_xml_file << "</spectrum>" << endl;

//      String spectrum_xlink_name = spectrum_name + String("_xlinker.txt");
//      spec_xml_file << "<spectrum filename=\"" << spectrum_xlink_name << "\" type=\"xlinker\">" << endl;
////      RichPeakSpectrum rich_xlink_spectrum = makeRichPeakSpectrum(preprocessed_pair_spectra.spectra_xlink_peaks[scan_index], true);
////      spec_xml_file <<TOPPxQuest::getxQuestBase64EncodedSpectrum_(rich_xlink_spectrum, spectrum_light_name + ".dta," + spectrum_heavy_name + ".dta");
//      spec_xml_file <<TOPPxQuest::getxQuestBase64EncodedSpectrum_(preprocessed_pair_spectra.spectra_xlink_peaks[pair_index], spectrum_light_name + ".dta," + spectrum_heavy_name + ".dta");
//      spec_xml_file << "</spectrum>" << endl;


        const RichPeakSpectrum& common_peaks = preprocessed_pair_spectra.spectra_common_peaks[pair_index];

        vector< CrossLinkSpectrumMatch > top_csms_spectrum;

        vector< double > top_preScore;
        vector< double > top_TIC;
        vector< double > top_wTIC;
        vector< double > top_intSum;
        vector< double > top_matchOdds;
        vector< vector< double > > top_xcorrx;
        vector< double > top_xcorrx_max;
        vector< vector< double > > top_xcorrc;
        vector< double > top_xcorrc_max;
        vector< double > top_score;
        vector< TheoreticalSpectrumGeneratorXLinks::ProteinProteinCrossLink > top_candidate_data;
        vector< Size > top_matched_spec_common_alpha;
        vector< Size > top_matched_spec_common_beta;
        vector< Size > top_matched_spec_xlink_alpha;
        vector< Size > top_matched_spec_xlink_beta;
        vector< Size > top_index;

        if(common_peaks.size() > 0 || preprocessed_pair_spectra.spectra_xlink_peaks[pair_index].size() > 0) // TODO: check if this is done in xQuest?
        {
          // determine candidates
          vector<pair<const AASequence*, const AASequence*> > candidates;
          if (ion_index_mode)
          {
            // Use 50 most intense common peaks of exp. spectrum, consider all peptides that produce any of these as theor. common ions
            //NLargest nlargest_filter = NLargest(50);
            RichPeakSpectrum common_peaks_50 = common_peaks;
            nlargest_filter_rich(common_peaks_50, 50);

            // old version with one most intensive peak
//            double most_intensive_peak_mz(0);
//            double most_intensive_peak_int(-1);

//            for (Size i = 0; i != common_peaks.size(); ++i)
//            {
//              double current_intensity = common_peaks[i].getIntensity();
//              double current_mz = common_peaks[i].getMZ();
//              if (current_intensity > most_intensive_peak_int)
//              {
//                most_intensive_peak_int = current_intensity;
//                most_intensive_peak_mz = current_mz;
//              }
//            }
//            const vector<AASequence*> ion_tag_candidates = hg.get(most_intensive_peak_mz, 5000);  //TODO: check if reduction to 5000 is already performed here or later after filtering


            vector<AASequence*> ion_tag_candidates;
            for (Size i = 0; i != common_peaks_50.size(); ++i)
            {
              const vector<AASequence*> new_ion_tag_candidates = hg.get(common_peaks_50[i].getMZ(), 5000);
              ion_tag_candidates.insert(ion_tag_candidates.end(), new_ion_tag_candidates.begin(), new_ion_tag_candidates.end());
            }
            sort(ion_tag_candidates.begin(), ion_tag_candidates.end());
            vector<AASequence*>::iterator last_unique = unique(ion_tag_candidates.begin(), ion_tag_candidates.end());
            ion_tag_candidates.erase(last_unique, ion_tag_candidates.end());
            cout << "Ion tag candidates before mass filtering: " << ion_tag_candidates.size() << endl;

            //vector<pair<AASequence, AASequence> > candidates;
            for (Size i = 0; i != ion_tag_candidates.size(); ++i)
            {
              const AASequence* peptide_a = ion_tag_candidates[i];
              
              for (Size j = i + 1; j < ion_tag_candidates.size(); ++j)
              {
                const AASequence* peptide_b = ion_tag_candidates[j];

                double cross_link_mass = peptide_a->getMonoWeight() + peptide_b->getMonoWeight() + cross_link_mass_light; //TODO: find a way to precalculate individual peptide masses
                double error_Da = abs(cross_link_mass - precursor_mass);
                if (error_Da < precursor_mass_tolerance)
                {
                  // TODO adapt to new candidates structure
                  candidates.push_back(make_pair(peptide_a, peptide_b));
                }                
              } 
            } 
            cout << "Ion tag candidates after mass filtering: " << candidates.size() << endl;
          } else // enumeration mode
          {
            //cout << "Number of common peaks, xlink peaks: " << preprocessed_pair_spectra.spectra_common_peaks[scan_index].size() << "\t" << preprocessed_pair_spectra.spectra_xlink_peaks[scan_index].size();
//            if (preprocessed_pair_spectra.spectra_common_peaks[scan_index].size() < 1 || preprocessed_pair_spectra.spectra_xlink_peaks[scan_index].size() < 1)
//            {
//              continue;
//            }
            // determine MS2 precursors that match to the current peptide mass
            multimap<double, pair<const AASequence*, const AASequence*> >::const_iterator low_it;
            multimap<double, pair<const AASequence*, const AASequence*> >::const_iterator up_it;

            if (precursor_mass_tolerance_unit_ppm) // ppm
            {
              low_it = enumerated_cross_link_masses.lower_bound(precursor_mass - (precursor_mass * precursor_mass_tolerance * 1e-6));
              up_it = enumerated_cross_link_masses.upper_bound(precursor_mass + (precursor_mass * precursor_mass_tolerance * 1e-6));
            }
            else // Dalton
            {
              low_it = enumerated_cross_link_masses.lower_bound(precursor_mass - precursor_mass_tolerance);
              up_it =  enumerated_cross_link_masses.upper_bound(precursor_mass + precursor_mass_tolerance);
            }

            if (!(low_it == up_it)) // no matching precursor in data
            {
              for (; low_it != up_it; ++low_it)
              {
                candidates.push_back(low_it->second);
              }
            }
          }

          // Find all positions of lysine (K) in the peptides (possible scross-linking sites), create cross_link_candidates with all combinations
          vector <TheoreticalSpectrumGeneratorXLinks::ProteinProteinCrossLink> cross_link_candidates;
          for (Size i = 0; i != candidates.size(); ++i)
          {
            pair<const AASequence*, const AASequence*> candidate = candidates[i];
            vector <SignedSize> K_pos_first;
            vector <SignedSize> K_pos_second;
            AASequence peptide_first = *candidate.first;
            AASequence peptide_second = *candidate.second;
//            AASequence peptide_first = AASequence::fromString(candidate.first);
//            AASequence peptide_second = AASequence::fromString(candidate.second);
            String seq_first = peptide_first.toUnmodifiedString();
            String seq_second =  peptide_second.toUnmodifiedString();




            for (Size k = 0; k < seq_first.size()-1; ++k)
            {
              if (seq_first[k] == 'K') K_pos_first.push_back(k);
            }
            if (seq_second.size() > 0)
            {
              for (Size k = 0; k < seq_second.size()-1; ++k)
              {
                if (seq_second[k] == 'K') K_pos_second.push_back(k);
              }
            } else
            {
              K_pos_second.push_back(-1);
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

            for (Size x = 0; x < K_pos_first.size(); ++x)
            {
              for (Size y = 0; y < K_pos_second.size(); ++y)
              {
                TheoreticalSpectrumGeneratorXLinks::ProteinProteinCrossLink cross_link_candidate;
                if (alpha_first)
                {
                  cross_link_candidate.alpha = peptide_first;
                  cross_link_candidate.beta = peptide_second;
                  cross_link_candidate.cross_link_position.first = K_pos_first[x];
                  cross_link_candidate.cross_link_position.second = K_pos_second[y];
                } else
                {
                  cross_link_candidate.alpha = peptide_second;
                  cross_link_candidate.beta = peptide_first;
                  cross_link_candidate.cross_link_position.first = K_pos_second[y];
                  cross_link_candidate.cross_link_position.second = K_pos_first[x];
                }
                cross_link_candidates.push_back(cross_link_candidate);
              }
            }
          }

          // lists for one spectrum, to determine best match to the spectrum
          vector< CrossLinkSpectrumMatch > all_csms_spectrum;

          vector< double > candidate_preScore;
          vector< double > candidate_TIC;
          vector< double > candidate_wTIC;
          vector< double > candidate_intSum;
          vector< double > candidate_matchOdds;
          vector< vector< double > > candidate_xcorrx;
          vector< double > candidate_xcorrx_max;
          vector< vector< double > > candidate_xcorrc;
          vector< double > candidate_xcorrc_max;
          vector< double > candidate_score;
          vector< TheoreticalSpectrumGeneratorXLinks::ProteinProteinCrossLink > candidate_data;
          vector< Size > candidate_matched_spec_common_alpha;
          vector< Size > candidate_matched_spec_common_beta;
          vector< Size > candidate_matched_spec_xlink_alpha;
          vector< Size > candidate_matched_spec_xlink_beta;
          vector< Size > candidate_index;

          // TODO variables for benchmarking and testing purposes
          if (cross_link_candidates.size() > maxMatchCount) maxMatchCount = cross_link_candidates.size();
          sumMatchCount += cross_link_candidates.size();

          for (Size i = 0; i != cross_link_candidates.size(); ++i)
          {
            TheoreticalSpectrumGeneratorXLinks::ProteinProteinCrossLink cross_link_candidate = cross_link_candidates[i];
            double candidate_mz = (cross_link_candidate.alpha.getMonoWeight() + cross_link_candidate.beta.getMonoWeight() +  cross_link_mass_light + (static_cast<double>(precursor_charge) * Constants::PROTON_MASS_U)) / precursor_charge;

            cout << "Pair: " << cross_link_candidate.alpha.toUnmodifiedString() << "-" << cross_link_candidate.beta.toUnmodifiedString() << " matched to light spectrum " << scan_index << "\t and heavy spectrum " << scan_index_heavy
              << " with m/z: " << spectrum_light.getPrecursors()[0].getMZ() << "\t" << "and candidate m/z: " << candidate_mz << "\tK Positions: " << cross_link_candidate.cross_link_position.first << "\t" << cross_link_candidate.cross_link_position.second << endl;
//          cout << a->second.getMonoWeight() << ", " << b->second.getMonoWeight() << " cross_link_mass_light: " <<  cross_link_mass_light <<  endl;

	    CrossLinkSpectrumMatch csm;
	    csm.cross_link = cross_link_candidate;

	    RichPeakSpectrum theoretical_spec_beta;
	    RichPeakSpectrum theoretical_spec_alpha;
	    RichPeakSpectrum theoretical_spec_xlinks_alpha;
	    RichPeakSpectrum theoretical_spec_xlinks_beta;

            //specGen.getSpectrum(theoretical_spec, cross_link_candidate, 1);
            //getXLinkIonSpectrum(theoretical_spec_xlinks , cross_link_candidate, 1)
            specGen.getCommonIonSpectrum(theoretical_spec_alpha, cross_link_candidate.alpha, 3);
            if (cross_link_candidate.cross_link_position.second != -1)
              {
              specGen.getCommonIonSpectrum(theoretical_spec_beta, cross_link_candidate.beta, 3);
              specGen.getXLinkIonSpectrum(theoretical_spec_xlinks_alpha, theoretical_spec_xlinks_beta, cross_link_candidate, 2, 7, cross_link_mass_light);
            } else
            {
              for (Size k = 0; k < cross_link_mass_mono_link.size(); ++k)
              {
                specGen.getXLinkIonSpectrum(theoretical_spec_xlinks_alpha, theoretical_spec_xlinks_beta, cross_link_candidate, 2, 7, cross_link_mass_mono_link[k]);
              }
            }

            vector< pair< Size, Size > > matched_spec_alpha;
            vector< pair< Size, Size > > matched_spec_beta;
            vector< pair< Size, Size > > matched_spec_xlinks_alpha;
            vector< pair< Size, Size > > matched_spec_xlinks_beta;

            if (preprocessed_pair_spectra.spectra_common_peaks[pair_index].size() > 0)
            {
              getSpectrumAlignment(matched_spec_alpha, theoretical_spec_alpha, preprocessed_pair_spectra.spectra_common_peaks[pair_index], fragment_mass_tolerance, false);
              getSpectrumAlignment(matched_spec_beta, theoretical_spec_beta, preprocessed_pair_spectra.spectra_common_peaks[pair_index], fragment_mass_tolerance, false);
            }
            if (preprocessed_pair_spectra.spectra_xlink_peaks[pair_index].size() > 0)
            {
              getSpectrumAlignment(matched_spec_xlinks_alpha, theoretical_spec_xlinks_alpha, preprocessed_pair_spectra.spectra_xlink_peaks[pair_index], fragment_mass_tolerance_xlinks, false);
              getSpectrumAlignment(matched_spec_xlinks_beta, theoretical_spec_xlinks_beta, preprocessed_pair_spectra.spectra_xlink_peaks[pair_index], fragment_mass_tolerance_xlinks, false);
            }

              // Pre-Score calculations
              Size matched_alpha_count = matched_spec_alpha.size() + matched_spec_xlinks_alpha.size();
              Size theor_alpha_count = theoretical_spec_alpha.size() + theoretical_spec_xlinks_alpha.size();
              Size matched_beta_count = matched_spec_beta.size() + matched_spec_xlinks_beta.size();
              Size theor_beta_count = theoretical_spec_beta.size() + theoretical_spec_xlinks_beta.size();

              if (matched_alpha_count + matched_beta_count > 0)
              {
                  // Simplified pre-Score
                  //float pre_score = preScore(matched_fragments_theor_spec.size(), theoretical_spec.size());
                  float pre_score = preScore(matched_alpha_count, theor_alpha_count, matched_beta_count, theor_beta_count);
                  //cout << "Number of matched peaks to theor. spectrum: " << matched_alpha_count << "\t" << matched_beta_count << endl;
                  //cout << "Number of theoretical ions: " << theor_alpha_count << "\t" << theor_beta_count << endl;
                  //cout << "Pre Score: " << pre_score << endl;
                  //cout << "Peptide size: " << a->second.size() << "\t" << b->second.size() << "\t" << "K Pos:" << K_pos_a[i] << "\t" << K_pos_b[i] << endl;
                  if (pre_score > pScoreMax) pScoreMax = pre_score;

                  // %TIC score calculations
                  double matched_current_alpha = 0;
                  double matched_current_beta = 0;
                  for (SignedSize j = 0; j < static_cast<SignedSize>(matched_spec_alpha.size()); ++j)
                  {
                    matched_current_alpha += preprocessed_pair_spectra.spectra_common_peaks[pair_index][matched_spec_alpha[j].second].getIntensity();
                   }
                  for (SignedSize j = 0; j < static_cast<SignedSize>(matched_spec_beta.size()); ++j)
                  {
                    matched_current_beta += preprocessed_pair_spectra.spectra_common_peaks[pair_index][matched_spec_beta[j].second].getIntensity();
                  }
                  for (SignedSize j = 0; j < static_cast<SignedSize>(matched_spec_xlinks_alpha.size()); ++j)
                  {
                    matched_current_alpha += preprocessed_pair_spectra.spectra_xlink_peaks[pair_index][matched_spec_xlinks_alpha[j].second].getIntensity();
                  }
                  for (SignedSize j = 0; j < static_cast<SignedSize>(matched_spec_xlinks_beta.size()); ++j)
                  {
                    matched_current_beta += preprocessed_pair_spectra.spectra_xlink_peaks[pair_index][matched_spec_xlinks_beta[j].second].getIntensity();
                  }
                  double matched_current = matched_current_alpha + matched_current_beta;

                  double total_current = 0;
                  for (SignedSize j = 0; j < static_cast<SignedSize>(spectrum_light.size()); ++j)
                  {
                    total_current += spectrum_light[j].getIntensity();
                   }

                  double TIC = matched_current / total_current;
                  //cout << "matched current: " << matched_current << "\t total_current: " << total_current << endl;
                  //cout << "%TIC: " << TIC << "%";
                  if (TIC > TICMax) TICMax = TIC;

                  double aatotal = cross_link_candidate.alpha.size() + cross_link_candidate.beta.size();
                  // TODO from xquest.def. is this the right parameter, should we use something different? Why digest length, what does it tell us about this specific case?
                  double maxdigestlength = 50;
                  double mindigestlength = 5;
//                  double invMax = 1 / (min(cross_link_candidate.alpha.size(), cross_link_candidate.beta.size()) / aatotal);
                  double invMax = 1 / (mindigestlength / (mindigestlength + maxdigestlength));
                  double invFrac_alpha = 1 / (cross_link_candidate.alpha.size() / aatotal);
                  double invFrac_beta = 1 / (cross_link_candidate.beta.size() / aatotal);
                  double TIC_weight_alpha = invFrac_alpha / invMax;
                  double TIC_weight_beta = invFrac_beta / invMax;
                  // (alpha weight * TIC_alpha) + (beta weight * TIC_beta)
                  double wTIC = TIC_weight_alpha * (matched_current_alpha / total_current ) + TIC_weight_beta * (matched_current_beta / total_current);
                  //cout << "\t wTIC: " << wTIC;
                  if (wTIC > wTICMax) wTICMax = wTIC;

                  //cout << "\t Intsum: " << matched_current << endl;
                  if (matched_current > intsumMax) intsumMax = matched_current;

                  // match-odds score
                  double range_c_alpha = theoretical_spec_alpha[theoretical_spec_alpha.size()-1].getMZ() -  theoretical_spec_alpha[0].getMZ();
                  double range_x_alpha= theoretical_spec_xlinks_alpha[theoretical_spec_xlinks_alpha.size()-1].getMZ() -  theoretical_spec_xlinks_alpha[0].getMZ();
                  double range_c_beta = 0;
                  double range_x_beta = 0;

                  if (cross_link_candidate.cross_link_position.second != -1)
                  {
                    range_c_beta = theoretical_spec_beta[theoretical_spec_beta.size()-1].getMZ() -  theoretical_spec_beta[0].getMZ();
                    range_x_beta = theoretical_spec_xlinks_beta[theoretical_spec_xlinks_beta.size()-1].getMZ() -  theoretical_spec_xlinks_beta[0].getMZ();
                  } else
                  {
                    range_c_beta = 0;
                    range_x_beta = 0;
                  }

                  double a_priori_ca = aPrioriProb(tolerance_binsize, theoretical_spec_alpha.size(), range_c_alpha, 3);
                  double a_priori_xa = aPrioriProb(tolerance_binsize, theoretical_spec_xlinks_alpha.size(), range_x_alpha, 3);
                  double a_priori_cb = aPrioriProb(tolerance_binsize, theoretical_spec_beta.size(), range_c_beta, 3);
                  double a_priori_xb = aPrioriProb(tolerance_binsize, theoretical_spec_xlinks_beta.size(), range_x_beta, 3);

                  double match_odds_c_alpha = 0;
                  double match_odds_x_alpha = 0;
                  double match_odds_c_beta = 0;
                  double match_odds_x_beta = 0;

                  if (matched_spec_alpha.size() == theoretical_spec_alpha.size())
                  {
                    //match_odds_c_alpha = 100;
                  } else if (matched_spec_alpha.size() > 0 && theoretical_spec_alpha.size() > 0)
                  {
                    match_odds_c_alpha = -log(1 - cumulativeBinomial(theoretical_spec_alpha.size(), matched_spec_alpha.size(), a_priori_ca) );
                  }

                  if (matched_spec_xlinks_alpha.size() == theoretical_spec_xlinks_alpha.size())
                  {
                    //match_odds_x_alpha = 100;
                  } else if (matched_spec_xlinks_alpha.size() > 0 && theoretical_spec_xlinks_alpha.size() > 0)
                  {
                    match_odds_x_alpha = -log(1 - cumulativeBinomial(theoretical_spec_xlinks_alpha.size(), matched_spec_xlinks_alpha.size(), a_priori_xa) );
                  }

                  if (matched_spec_beta.size() == theoretical_spec_beta.size())
                  {
                    //match_odds_c_beta = 100;
                  } else if (matched_spec_beta.size() > 0 && theoretical_spec_beta.size() > 0)
                  {
                    match_odds_c_beta = -log(1 - cumulativeBinomial(theoretical_spec_beta.size(), matched_spec_beta.size(), a_priori_cb) );
                  }

                  if (matched_spec_xlinks_beta.size() == theoretical_spec_xlinks_beta.size())
                  {
                    //match_odds_x_beta = 100;
                  } else if (matched_spec_xlinks_beta.size() > 0 && theoretical_spec_xlinks_beta.size() > 0)
                  {
                    match_odds_x_beta = -log(1 - cumulativeBinomial(theoretical_spec_xlinks_beta.size(), matched_spec_xlinks_beta.size(), a_priori_xb) );
                  }

                  double match_odds = (match_odds_c_alpha + match_odds_x_alpha + match_odds_c_beta + match_odds_x_beta) / 4;

                  // Debug output
                  if (match_odds == INFINITY)
                  {
                    cout << "INFINITY Match-odds (bug, should not happen) \n";
                    cout << "Match-odds Scores" << match_odds_c_alpha << "\t" << match_odds_x_alpha << "\t" << match_odds_c_beta << "\t" << match_odds_x_beta << endl;
                    cout << "Matched, Theo sizes: " << matched_spec_beta.size() << "\t" << theoretical_spec_beta.size() << "\t cumulBinom: " << cumulativeBinomial(theoretical_spec_beta.size(), matched_spec_beta.size(), a_priori_cb) << "\t a priori: " << a_priori_cb << endl;
                    return EXECUTION_OK;
                  }

                  //cout << "Range Alpha: " << theoretical_spec_alpha[theoretical_spec_alpha.size()-1].getMZ()  << " - " << (double) theoretical_spec_alpha[0].getMZ() << " = " << range_c_alpha << "\t A Priori probaility common alpha: " << a_priori_ca << endl;
                  //cout << "Match-Odds Scores: " << match_odds_c_alpha << "\t" << match_odds_x_alpha << "\t" << match_odds_c_beta << "\t" << match_odds_x_beta << "\t Match-Odds Final: " << match_odds << endl;
                  if (match_odds > matchOddsMax) matchOddsMax = match_odds;


                  //Cross-correlation
                  RichPeakSpectrum theoretical_spec_common;
                  RichPeakSpectrum theoretical_spec_xlinks;
                  //cout << "Build theor. Spec." << endl;

                  theoretical_spec_common.reserve(theoretical_spec_alpha.size() + theoretical_spec_beta.size());
                  theoretical_spec_xlinks.reserve( theoretical_spec_xlinks_alpha.size() + theoretical_spec_xlinks_beta.size());
                  theoretical_spec_common.insert(theoretical_spec_common.end(), theoretical_spec_alpha.begin(), theoretical_spec_alpha.end());
                  theoretical_spec_common.insert(theoretical_spec_common.end(), theoretical_spec_beta.begin(), theoretical_spec_beta.end());
                  theoretical_spec_xlinks.insert(theoretical_spec_xlinks.end(), theoretical_spec_xlinks_alpha.begin(), theoretical_spec_xlinks_alpha.end());
                  theoretical_spec_xlinks.insert(theoretical_spec_xlinks.end(), theoretical_spec_xlinks_beta.begin(), theoretical_spec_xlinks_beta.end());
                  theoretical_spec_common.sortByPosition();
                  theoretical_spec_xlinks.sortByPosition();

                  //cout << "Compute xCorr" << endl;
                  // xCorr parameters: spec1, spec2, max_shift, binsize
                  vector< double > xcorrx = xCorrelation(spectrum_light, theoretical_spec_xlinks, 5, 0.3);
                  vector< double > xcorrc = xCorrelation(spectrum_light, theoretical_spec_common, 5, 0.2);
                  double xcorrx_max = *max_element(xcorrx.begin(), xcorrx.end());
                  double xcorrc_max = *max_element(xcorrc.begin(), xcorrc.end());
                  //cout << "XLink Cross correlation score: " << xcorrx_max << "\t Common Cross correlation score: " << xcorrc_max << endl;
                  if (xcorrx_max > xcorrxMax) xcorrxMax = xcorrx_max;
                  if (xcorrc_max > xcorrcMax) xcorrcMax = xcorrc_max;

                  // Compute score from the 4 scores and 4 weights
                  double xcorrx_weight = 2.488;
                  double xcorrc_weight = 21.279;
                  double match_odds_weight = 1.973;
                  double wTIC_weight = 12.829;
                  double intsum_weight = 18;

                  double score = xcorrx_weight * xcorrx_max + xcorrc_weight * xcorrc_max + match_odds_weight * match_odds + wTIC_weight * wTIC + intsum_weight * matched_current;

                  csm.score = score;
                  csm.pre_score = pre_score;
                  csm.percTIC = TIC;
                  csm.wTIC = wTIC;
                  csm.int_sum = matched_current;
                  csm.match_odds = match_odds;
                  csm.xcorrx = xcorrx;
                  csm.xcorrx_max = xcorrx_max;
                  csm.xcorrc = xcorrc;
                  csm.xcorrc_max = xcorrc_max;
                  csm.matched_common_alpha = matched_spec_alpha.size();
                  csm.matched_common_beta = matched_spec_beta.size();
                  csm.matched_xlink_alpha = matched_spec_xlinks_alpha.size();
                  csm.matched_xlink_beta = matched_spec_xlinks_beta.size();
                  csm.scan_index_light = scan_index;
                  csm.scan_index_heavy = scan_index_heavy;
                  all_csms_spectrum.push_back(csm);

                  candidate_preScore.push_back(pre_score);
                  candidate_TIC.push_back(TIC);
                  candidate_wTIC.push_back(wTIC);
                  candidate_intSum.push_back(matched_current);
                  candidate_matchOdds.push_back(match_odds);
                  candidate_xcorrx.push_back(xcorrx);
                  candidate_xcorrx_max.push_back((xcorrx_max));
                  candidate_xcorrc.push_back(xcorrc);
                  candidate_xcorrc_max.push_back((xcorrc_max));
                  candidate_score.push_back(score);
                  candidate_data.push_back(cross_link_candidate);
                  candidate_matched_spec_common_alpha.push_back(matched_spec_alpha.size());
                  candidate_matched_spec_common_beta.push_back(matched_spec_beta.size());
                  candidate_matched_spec_xlink_alpha.push_back(matched_spec_xlinks_alpha.size());
                  candidate_matched_spec_xlink_beta.push_back(matched_spec_xlinks_beta.size());
                  candidate_index.push_back(scan_index);

                  //cout << "Next candidate: " << endl;
              }
            } // candidates for peak finished, determine best matching candidate

            Int top = 0;

            // collect top 5 matches to spectrum
            while(!all_csms_spectrum.empty() && top < 5)
            {
              top++;

              //double max_score = *max_element(candidate_score.begin(), candidate_score.end());
              Int max_position = distance(all_csms_spectrum.begin(), max_element(all_csms_spectrum.begin(), all_csms_spectrum.end()));
//              cout << "Max_position_csms: " << max_position << endl;
//              max_position = distance(candidate_score.begin(), max_element(candidate_score.begin(), candidate_score.end()));
//              cout << "Max_position_vectors: " << max_position << endl;

              all_csms_spectrum[max_position].rank = top;
              cout << "all_csms_spectrum[max] score: " << all_csms_spectrum[max_position].score << endl;
              top_csms_spectrum.push_back(all_csms_spectrum[max_position]);
              all_csms_spectrum.erase(all_csms_spectrum.begin() + max_position);

              top_preScore.push_back(candidate_preScore[max_position]);
              top_TIC.push_back(candidate_TIC[max_position]);
              top_wTIC.push_back(candidate_wTIC[max_position]);
              top_intSum.push_back(candidate_intSum[max_position]);
              top_matchOdds.push_back(candidate_matchOdds[max_position]);
              top_xcorrx.push_back(candidate_xcorrx[max_position]);
              top_xcorrx_max.push_back(candidate_xcorrx_max[max_position]);
              top_xcorrc.push_back(candidate_xcorrc[max_position]);
              top_xcorrc_max.push_back(candidate_xcorrc_max[max_position]);
              top_score.push_back(candidate_score[max_position]);
              top_candidate_data.push_back(candidate_data[max_position]);
              top_matched_spec_common_alpha.push_back(candidate_matched_spec_common_alpha[max_position]);
              top_matched_spec_common_beta.push_back(candidate_matched_spec_common_beta[max_position]);
              top_matched_spec_xlink_alpha.push_back(candidate_matched_spec_xlink_alpha[max_position]);
              top_matched_spec_xlink_beta.push_back(candidate_matched_spec_xlink_beta[max_position]);
              top_index.push_back(candidate_index[max_position]);

              candidate_preScore.erase(candidate_preScore.begin() + max_position);
              candidate_TIC.erase(candidate_TIC.begin() + max_position);
              candidate_wTIC.erase(candidate_wTIC.begin() + max_position);
              candidate_intSum.erase(candidate_intSum.begin() + max_position);
              candidate_matchOdds.erase(candidate_matchOdds.begin() + max_position);
              candidate_xcorrx.erase(candidate_xcorrx.begin() + max_position);
              candidate_xcorrx_max.erase(candidate_xcorrx_max.begin() + max_position);
              candidate_xcorrc.erase(candidate_xcorrc.begin() + max_position);
              candidate_xcorrc_max.erase(candidate_xcorrc_max.begin() + max_position);
              candidate_score.erase(candidate_score.begin() + max_position);
              candidate_data.erase(candidate_data.begin() + max_position);
              candidate_matched_spec_common_alpha.erase(candidate_matched_spec_common_alpha.begin() + max_position);
              candidate_matched_spec_common_beta.erase(candidate_matched_spec_common_beta.begin() + max_position);
              candidate_matched_spec_xlink_alpha.erase(candidate_matched_spec_xlink_alpha.begin() + max_position);
              candidate_matched_spec_xlink_beta.erase(candidate_matched_spec_xlink_beta.begin() + max_position);
              candidate_index.erase(candidate_index.begin() + max_position);
            }
            all_top_csms.push_back(top_csms_spectrum);
            if (top_csms_spectrum.size() > 0)
            {
              cout << "Spectrum top score (csm):  " << top_csms_spectrum[0].score << "\t" << all_top_csms[all_top_csms.size()-1][0].score << endl;
            }

            cout << "Spectrum scores: " << top_score << endl;
            cout << "Spectrum preScores: " << top_preScore << endl;
            cout << "Spectrum wTICs: " << top_wTIC << endl;
            cout << "Spectrum xcorrx: " << top_xcorrx_max << endl;
            cout << "Spectrum xcorrc: " << top_xcorrc_max << endl;
            cout << "Spectrum match-odds: " << top_matchOdds << endl;
            cout << "Spectrum Intsum: " << top_intSum << endl;
            cout << "Spectrum matched ions calpha , cbeta , xalpha , xbeta: " << top_matched_spec_common_alpha << "\t" << top_matched_spec_common_beta <<
              "\t" << top_matched_spec_xlink_alpha << "\t" << top_matched_spec_xlink_beta << endl;

            // Write top 5 hits to file
            for (Size i = 0; i < top_csms_spectrum.size(); ++i)
            {
              PeptideIdentification peptide_id;

              String xltype = "monolink";
              String structure = top_csms_spectrum[i].cross_link.alpha.toUnmodifiedString();
              String topology = String("K") + (top_csms_spectrum[i].cross_link.cross_link_position.first+1);

              // TODO track or otherwise find out, which kind of mono-link it was (if there are several possibilities for the weigths)
              double weight = top_csms_spectrum[i].cross_link.alpha.getMonoWeight() + cross_link_mass_mono_link[1];

              bool is_monolink = (top_csms_spectrum[i].cross_link.cross_link_position.second == -1);
              int alpha_pos = top_csms_spectrum[i].cross_link.cross_link_position.first + 1;
              int beta_pos = top_csms_spectrum[i].cross_link.cross_link_position.second + 1;

              if (!is_monolink)
              {
                xltype = "xlink";
                structure = top_csms_spectrum[i].cross_link.alpha.toUnmodifiedString() + "-" + top_csms_spectrum[i].cross_link.beta.toUnmodifiedString();
                topology = String("a") +  alpha_pos + String("-b") + beta_pos;
                weight = top_csms_spectrum[i].cross_link.alpha.getMonoWeight() + top_csms_spectrum[i].cross_link.beta.getMonoWeight() + cross_link_mass_light;
              }
              String id = structure + "-" + topology;

              vector<PeptideHit> phs;
              PeptideHit ph_alpha, ph_beta;
              ph_alpha.setSequence(top_csms_spectrum[i].cross_link.alpha);
              ph_alpha.setCharge(precursor_charge);
              ph_alpha.setScore(top_csms_spectrum[i].score);
              ph_alpha.setMetaValue("xl_chain", "MS:1002510");  // receiver
              ph_alpha.setMetaValue("xl_pos", DataValue(alpha_pos));
              phs.push_back(ph_alpha);

              if (!is_monolink)
              {
                ph_beta.setSequence(top_csms_spectrum[i].cross_link.beta);
                ph_beta.setCharge(precursor_charge);
                ph_beta.setScore(top_csms_spectrum[i].score);
                ph_beta.setMetaValue("xl_chain", "MS:1002509"); // donor 
                ph_beta.setMetaValue("xl_pos", DataValue(beta_pos));
                phs.push_back(ph_beta);
              }

              peptide_id.setRT(spectrum_light.getRT());
              peptide_id.setMZ(precursor_mz);
              peptide_id.setMetaValue("spec_heavy_RT", spectra[scan_index_heavy].getRT());
              peptide_id.setMetaValue("spec_heavy_MZ", spectra[scan_index_heavy].getPrecursors()[0].getMZ());
              peptide_id.setMetaValue("xl_type", xltype); // TODO: needs CV term
              peptide_id.setMetaValue("xl_rank", DataValue(i + 1)); 

//            peptide_id.setMetaValue("xl_relation", ); //TODO: needs CV term
              peptide_id.setHits(phs);
              cout << "peptide_ids number: " << peptide_ids.size() << endl;
              peptide_ids.push_back(peptide_id);
              all_top_csms[all_top_csms.size()-1][i].peptide_id_index = peptide_ids.size()-1;
              all_top_csms[all_top_csms.size()-1][i].peptide_id = &peptide_ids[peptide_ids.size()-1];

//              xml_file << "<search_hit search_hit_rank=\"" << i+1 << "\" id=\"" << id << "\" type=\"" << xltype << "\" structure=\"" << structure << "\" seq1=\"" << top_csms_spectrum[i].cross_link.alpha.toUnmodifiedString() << "\" seq2=\"" << top_csms_spectrum[i].cross_link.beta.toUnmodifiedString()
//                << "\" prot1=\"" << "TODO" << "\" prot2=\"" << "TODO" << "\" topology=\"" << topology << "\" xlinkposition=\"" << (top_csms_spectrum[i].cross_link.cross_link_position.first+1) << "," << (top_csms_spectrum[i].cross_link.cross_link_position.second+1)
//                << "\" Mr=\"" << weight << "\" mz=\"" << precursor_mz << "\" charge=\"" << precursor_charge << "\" xlinkermass=\"" << cross_link_mass_light << "\" measured_mass=\"" << precursor_mass << "\" error=\"" << "TODO"
//                << "\" error_rel=\"" << "TODO" << "\" xlinkions_matched=\"" << (top_csms_spectrum[i].matched_xlink_alpha + top_csms_spectrum[i].matched_xlink_beta) << "\" backboneions_matched=\"" << (top_csms_spectrum[i].matched_common_alpha + top_csms_spectrum[i].matched_common_beta)
//                << "\" weighted_matchodds_mean=\"" << "TODO" << "\" weighted_matchodds_sum=\"" << "TODO" << "\" match_error_mean=\"" << "TODO" << "\" match_error_stdev=\"" << "TODO" << "\" xcorrx=\"" << top_csms_spectrum[i].xcorrx_max << "\" xcorrb=\"" << top_csms_spectrum[i].xcorrc_max << "\" match_odds=\"" << top_csms_spectrum[i].match_odds << "\" prescore=\"" << top_csms_spectrum[i].pre_score
//                << "\" prescore_alpha=\"" << "TODO" << "\" prescore_beta=\"" << "TODO" << "\" match_odds_alphacommon=\"" << "TODO" << "\" match_odds_betacommon=\"" << "TODO" << "\" match_odds_alphaxlink=\"" << "TODO"
//                << "\" match_odds_betaxlink=\"" << "TODO" << "\" num_of_matched_ions_alpha=\"" << (top_csms_spectrum[i].matched_common_alpha + top_csms_spectrum[i].matched_xlink_alpha) << "\" num_of_matched_ions_beta=\"" << (top_csms_spectrum[i].matched_common_beta + top_csms_spectrum[i].matched_xlink_beta) << "\" num_of_matched_common_ions_alpha=\"" << top_csms_spectrum[i].matched_common_alpha
//                << "\" num_of_matched_common_ions_beta=\"" << top_csms_spectrum[i].matched_common_beta << "\" num_of_matched_xlink_ions_alpha=\"" << top_matched_spec_xlink_alpha[i] << "\" num_of_matched_xlink_ions_beta=\"" << top_csms_spectrum[i].matched_xlink_beta << "\" xcorrall=\"" << "TODO" << "\" TIC=\"" << top_csms_spectrum[i].percTIC
//                << "\" TIC_alpha=\"" << "TODO" << "\" TIC_beta=\"" << "TODO" << "\" wTIC=\"" << top_csms_spectrum[i].wTIC << "\" intsum=\"" << top_csms_spectrum[i].int_sum << "\" apriori_match_probs=\"" << "TODO" << "\" apriori_match_probs_log=\"" << "TODO"
//                << "\" series_score_mean=\"" << "TODO" << "\" annotated_spec=\"" << "" << "\" score=\"" << top_csms_spectrum[i].score << "\" >" << endl;
//              xml_file << "</search_hit>" << endl;

            }


            cout << "Next Spectrum ################################## \n";
          }
//      xml_file << "</spectrum_search>" << endl;
    }
    for (Size i = 0; i < all_top_csms.size(); ++i)
    {
      vector< CrossLinkSpectrumMatch > top5 = all_top_csms[i];
      cout << "Pair Index: " << i << endl;
        for (Size k = 0; k < top5.size(); ++k)
        {
          cout << "Peptide_index: " << all_top_csms[i][k].peptide_id_index <<  "\t vector pos: " << i << "\t" << k << endl;
          cout << "Peptide_mz: " << all_top_csms[i][k].peptide_id->getMZ() << endl;
          cout << "Peptide_index: " << top5[k].peptide_id_index << endl;
          cout << "Peptide_mz: " << top5[k].peptide_id->getMZ() << endl;
        }
    }

    // end of matching / scoring
    progresslogger.endProgress();

//    xml_file << "</xquest_results>" << endl;
//    xml_file.close();

//    spec_xml_file << "</xquest_spectra>" << endl;
//    spec_xml_file.close();

    cout << "Pre Score maximum: " << pScoreMax << "\t TIC maximum: " << TICMax << "\t wTIC maximum: " << wTICMax << "\t Match-Odds maximum: " << matchOddsMax << endl;
    cout << "XLink Cross-correlation maximum: " << xcorrxMax << "\t Common Cross-correlation maximum: " << xcorrcMax << "\t Intsum maximum: " << intsumMax << endl;
    cout << "Total number of matched candidates: " << sumMatchCount << "\t Maximum number of matched candidates to one spectrum pair: " << maxMatchCount << "\t Average: " << sumMatchCount / spectra.size() << endl;
    //cout << "Random Charge: " << spectra[10][15].getMetaValue("z") << endl;

    // Add protein identifications
    PeptideIndexing pep_indexing;
    Param indexing_param = pep_indexing.getParameters();

    // TODO add additional parameters of PeptideIndexing to this tool (enzyme etc.)
    String d_prefix = decoy_prefix ? "true" : "false";
    indexing_param.setValue("prefix", d_prefix, "If set, protein accessions in the database contain 'decoy_string' as prefix.");
    indexing_param.setValue("decoy_string", decoy_string, "String that was appended (or prefixed - see 'prefix' flag below) to the accessions in the protein database to indicate decoy proteins.");
    pep_indexing.setParameters(indexing_param);

//    cout << "peptide_ids size 1 : " << peptide_ids.size() << endl;
    pep_indexing.run(fasta_db, protein_ids, peptide_ids);


    // write cross-links to IdXML and xquest.xml, write spectra for xquest visualization
    String spec_xml_name = getStringOption_("in").prefix(getStringOption_("in").size()-7) + "_matched";
    String spec_xml_filename = spec_xml_name + ".spec.xml";

    progresslogger.startProgress(0, 1, "Writing output...");
    IdXMLFile().store(out_idxml, protein_ids, peptide_ids);
    writeXQuestXML(out_xquest, all_top_csms, peptide_ids, spectra, spec_xml_name, cross_link_mass_light, cross_link_mass_mono_link);
    writeXQuestXMLSpec(spec_xml_filename, spectra, preprocessed_pair_spectra, spectrum_pairs, all_top_csms);
    progresslogger.endProgress();

 
    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{

  TOPPxQuest tool;
  
  return tool.main(argc, argv);
}

