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

// preprocessing and filtering
#include <OpenMS/FILTERING/TRANSFORMERS/ThresholdMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/NLargest.h>
#include <OpenMS/FILTERING/TRANSFORMERS/WindowMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignment.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGeneratorXLinks.h>

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

    registerTOPPSubsection_("precursor", "Precursor (Parent Ion) Options");
    registerDoubleOption_("precursor:mass_tolerance", "<tolerance>", 10.0, "Width of precursor mass tolerance window", false);

    StringList precursor_mass_tolerance_unit_valid_strings;
    precursor_mass_tolerance_unit_valid_strings.push_back("ppm");
    precursor_mass_tolerance_unit_valid_strings.push_back("Da");

    registerStringOption_("precursor:mass_tolerance_unit", "<unit>", "ppm", "Unit of precursor mass tolerance.", false, false);
    setValidStrings_("precursor:mass_tolerance_unit", precursor_mass_tolerance_unit_valid_strings);

    registerIntOption_("precursor:min_charge", "<num>", 3, "Minimum precursor charge to be considered.", false, true);
    registerIntOption_("precursor:max_charge", "<num>", 6, "Maximum precursor charge to be considered.", false, true);

    registerTOPPSubsection_("fragment", "Fragments (Product Ion) Options");
    registerDoubleOption_("fragment:mass_tolerance", "<tolerance>", 800, "Fragment mass tolerance", false);
    registerDoubleOption_("fragment:mass_tolerance_xlinks", "<tolerance>", 1000, "Fragment mass tolerance for cross-link ions", false);

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
    registerIntOption_("peptide:min_size", "<num>", 7, "Minimum size a peptide must have after digestion to be considered in the search.", false, true);
    registerIntOption_("peptide:missed_cleavages", "<num>", 1, "Number of missed cleavages.", false, false);
    vector<String> all_enzymes;
    EnzymesDB::getInstance()->getAllNames(all_enzymes);
    registerStringOption_("peptide:enzyme", "<cleavage site>", "Trypsin", "The enzyme used for peptide digestion.", false);
    setValidStrings_("peptide:enzyme", all_enzymes);


    registerTOPPSubsection_("cross_linker", "Cross Linker Options");
    registerDoubleOption_("cross_linker:mass_light", "<mass>", 156.078644, "Mass of the light cross-linker", false);
    registerDoubleOption_("cross_linker:mass_heavy", "<mass>", 168.153965, "Mass of the heavy cross-linker", false);
    registerDoubleOption_("cross_linker:mass_loss_type2", "<mass>", 18.01056, "Mass difference observed in an intra or inter peptide link", false);

    registerTOPPSubsection_("algorithm", "Algorithm Options");
    registerStringOption_("algorithm:candidate_search", "<param>", "index", "Mode used to generate candidate peptides.", false, false);
    StringList candidate_search_modes_strings;
    candidate_search_modes_strings.push_back("index");
    candidate_search_modes_strings.push_back("enumeration");
    setValidStrings_("algorithm:candidate_search", candidate_search_modes_strings);

    // output file
    registerOutputFile_("out", "<file>", "", "Result file\n");
    setValidFormats_("out", ListUtils::create<String>("xml"));
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



  void preprocessSpectra_(PeakMap& exp)
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
    WindowMower window_mower_filter;
    Param filter_param = window_mower_filter.getParameters();
    filter_param.setValue("windowsize", 100.0, "The size of the sliding window along the m/z axis.");
    filter_param.setValue("peakcount", 20, "The number of peaks that should be kept.");
    filter_param.setValue("movetype", "jump", "Whether sliding window (one peak steps or jumping window window size steps) should be used.");
    window_mower_filter.setParameters(filter_param);
    NLargest nlargest_filter = NLargest(250);   // De-noising in xQuest: Dynamic range = 1000, 250 most intense peaks
  
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (SignedSize exp_index = 0; exp_index < static_cast<SignedSize>(exp.size()); ++exp_index)
    {
      // sort by mz, for window_mower
      //exp[exp_index].sortByPosition();
  
      // remove noise, TODO window_mower not used in xQuest, is it necessary?
      //window_mower_filter.filterPeakSpectrum(exp[exp_index]);
      nlargest_filter.filterPeakSpectrum(exp[exp_index]);
  
      // sort (nlargest changes order)
      exp[exp_index].sortByPosition();
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
  std::vector< double > xCorrelation(const PeakSpectrum & spec1, const RichPeakSpectrum & spec2, Int maxshift, double tolerance)
  {
    double maxionsize = std::max(spec1[spec1.size()-1].getMZ(), spec2[spec2.size()-1].getMZ());
    //cout << "xcorr Maxionssize: " << maxionsize << endl;
    Int table_size = std::ceil(maxionsize / tolerance)+1;
    //cout << "xcorr table_size: " << table_size << endl;
    std::vector< double > ion_table1(table_size, 0);
    std::vector< double > ion_table2(table_size, 0);

    // Build tables of the same size, each bin has the size of the tolerance
    for (Size i = 0; i < spec1.size(); ++i)
    {
      Size pos = static_cast<Size>(std::ceil(spec1[i].getMZ() / tolerance));
      //cout << "xcorr Table1 pos: " << pos << endl;
      ion_table1[pos] = spec1[i].getIntensity() * 100;
      //cout << "xcorr Table1 Inte: " << spec1[i].getIntensity() << endl;
    }
    for (Size i = 0; i < spec2.size(); ++i)
    {
      Size pos =static_cast<Size>(std::ceil(spec2[i].getMZ() / tolerance));
      //cout << "xcorr Table2 pos: " << pos << endl;
      ion_table2[pos] = spec2[i].getIntensity() * 100;
      //cout << "xcorr Table2 Inte: " << spec2[i].getIntensity() << endl;
    }
    //cout << "xcorr Tables done" << endl;

    // Compute means
    double mean1 = (std::accumulate(ion_table1.begin(), ion_table1.end(), 0.0)) / table_size;
    double mean2 = (std::accumulate(ion_table2.begin(), ion_table2.end(), 0.0)) / table_size;
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
    std::vector< double > results(maxshift * 2 + 1, 0);
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
    PeakMap spectra_light_different; // peaks in light spectrum after common peaks have been removed
    PeakMap spectra_heavy_different; // peaks in heavy spectrum after common peaks have been removed
    PeakMap spectra_heavy_to_light; // heavy peaks transformed to light ones and after common peaks have been removed
    PeakMap spectra_common_peaks; // merge spectrum of common peaks (present in both spectra)
    PeakMap spectra_xlink_peaks; // Xlink peaks in the light spectrum (common peaks between spectra_light_different and spectra heavy_to_light)

    PreprocessedPairSpectra_(Size size)
    {
      for (Size i = 0; i != size; ++i)
      {
        spectra_light_different.addSpectrum(PeakSpectrum());
        spectra_heavy_different.addSpectrum(PeakSpectrum());
        spectra_heavy_to_light.addSpectrum(PeakSpectrum());
        spectra_common_peaks.addSpectrum(PeakSpectrum());
        spectra_xlink_peaks.addSpectrum(PeakSpectrum());
      }
    }  
  };


  // create common / shifted peak spectra for all pairs
  PreprocessedPairSpectra_ preprocessPairs_(const PeakMap& spectra, const map<Size, Size>& map_light_to_heavy, const SpectrumAlignment& ms2_aligner, const double cross_link_mass_light, const double cross_link_mass_heavy)
  {
    PreprocessedPairSpectra_ preprocessed_pair_spectra(spectra.size());
 
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (SignedSize scan_index = 0; scan_index < static_cast<SignedSize>(spectra.size()); ++scan_index)
    {
      // assume that current MS2 corresponds to a light spectrum
      const PeakSpectrum& spectrum_light = spectra[scan_index];

      // map light to heavy
      map<Size, Size>::const_iterator scan_index_light_it = map_light_to_heavy.find(scan_index);

      if (scan_index_light_it != map_light_to_heavy.end())
      {
        const Size scan_index_heavy = scan_index_light_it->second;
        const PeakSpectrum& spectrum_heavy = spectra[scan_index_heavy];
        std::vector< std::pair< Size, Size > > matched_fragments_without_shift;
        ms2_aligner.getSpectrumAlignment(matched_fragments_without_shift, spectrum_light, spectrum_heavy);

        // different fragments may carry light or heavy cross-linker.             
        PeakSpectrum spectrum_heavy_different;
        PeakSpectrum spectrum_light_different;

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
        // TODO: important: assume different charged MS2 fragments. Now only single charged ones are assumed

        // transform heavy spectrum to light spectrum
        PeakSpectrum spectrum_heavy_to_light;
        for (Size i = 0; i != spectrum_heavy_different.size(); ++i)
        {
          Peak1D p = spectrum_heavy_different[i];
          p.setMZ(p.getMZ() - (cross_link_mass_heavy - cross_link_mass_light));
          spectrum_heavy_to_light.push_back(p); 
        }

        // align potentially shifted peaks from light MS2 with potentially shifted peaks from heavy (after transformation to resemble the light MS2)
        // matching fragments are potentially carrying the cross-linker
        std::vector< std::pair< Size, Size > > matched_fragments_with_shift;

        ms2_aligner.getSpectrumAlignment(matched_fragments_with_shift, spectrum_light_different, spectrum_heavy_to_light);


        PeakSpectrum xlink_peaks;
        for (Size i = 0; i != matched_fragments_with_shift.size(); ++i)
        {
          xlink_peaks.push_back(spectrum_light[matched_fragments_with_shift[i].first]);
        }


#ifdef DEBUG_XQUEST
        cout << "Common peaks: " << matched_fragments_without_shift.size() << " different peaks: " << spectrum_light.size() - matched_fragments_without_shift.size() << ", " << spectrum_heavy.size() - matched_fragments_without_shift.size() << endl;
        cout << "Matched shifted peaks: " << matched_fragments_with_shift.size() << " unexplained peaks: " << spectrum_light_different.size() - matched_fragments_with_shift.size() << ", " << spectrum_heavy_to_light.size() - matched_fragments_with_shift.size() << endl;
#endif
        // generate common peaks spectrum TODO: check if only light / only heavy or e.g. mean of both peaks are used here
        PeakSpectrum common_peaks;
        for (Size i = 0; i != matched_fragments_without_shift.size(); ++i)
        {
          common_peaks.push_back(spectrum_light[matched_fragments_without_shift[i].first]);
        }
#ifdef DEBUG_XQUEST
        cout << "Peaks to match: " << common_peaks.size() << endl;
#endif

        std::swap(preprocessed_pair_spectra.spectra_light_different[scan_index], spectrum_light_different);
        std::swap(preprocessed_pair_spectra.spectra_heavy_different[scan_index], spectrum_heavy_different);
        std::swap(preprocessed_pair_spectra.spectra_heavy_to_light[scan_index], spectrum_heavy_to_light);
        std::swap(preprocessed_pair_spectra.spectra_common_peaks[scan_index], common_peaks);
        std::swap(preprocessed_pair_spectra.spectra_xlink_peaks[scan_index], xlink_peaks);
        preprocessed_pair_spectra.spectra_light_different[scan_index].sortByPosition();
        preprocessed_pair_spectra.spectra_heavy_different[scan_index].sortByPosition();
        preprocessed_pair_spectra.spectra_heavy_to_light[scan_index].sortByPosition();
        preprocessed_pair_spectra.spectra_common_peaks[scan_index].sortByPosition();
        preprocessed_pair_spectra.spectra_xlink_peaks[scan_index].sortByPosition();

        // Debug support output
        /*
        cout << "spectrum_light_different: " << preprocessed_pair_spectra.spectra_light_different[scan_index].size() << endl;
        cout << "spectrum_heavy_different: " << preprocessed_pair_spectra.spectra_heavy_different[scan_index].size() << endl;
        cout << "spectrum_heavy_to_light alignment: " << preprocessed_pair_spectra.spectra_heavy_to_light[scan_index].size() << endl;
        cout << "spctrum_common_peaks: " << preprocessed_pair_spectra.spectra_common_peaks[scan_index].size() << endl;
        cout << "spectrum_xlink_peaks: " << preprocessed_pair_spectra.spectra_xlink_peaks[scan_index].size() << endl;
        */
      }
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
        std::cerr << "Trying to add element left or right of allowed bounds. (min, max, position): " << min_ << ", " << max_ << ", " << position << std::endl;
        return;
      }

      const double bucket_index = (position - min_) / bucket_size_;
      h_.insert(make_pair(bucket_index, v));
    }

    vector<AASequence*> get(double position, Size max_elements = std::numeric_limits<Size>::max()) 
    {
      if (position < min_ || position > max_) 
      {
        std::cerr << "Trying to access element left or right of allowed bounds. (min, max, position): " << min_ << ", " << max_ << ", " << position << std::endl;
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


  static multimap<double, pair<AASequence, AASequence> > enumerateCrossLinksAndMasses_(const multimap<StringView, AASequence>&  peptides, double cross_link_mass_light, double cross_link_mass_loss_type2)
  {
    multimap<double, pair<AASequence, AASequence> > mass_to_candidates;
    for (map<StringView, AASequence>::const_iterator a = peptides.begin(); a != peptides.end(); ++a)
    {
      if (a->second.toUnmodifiedString().find("K") >= a->second.size()-1)
      {
        continue;
      }

      for (map<StringView, AASequence>::const_iterator b = a; b != peptides.end(); ++b)
      {

        if (b->second.toUnmodifiedString().find("K") >= b->second.size()-1)
        {
          continue;
        }
/*
        // Find all positions of lysine (K) in the peptides (possible scross-linking sites)
        vector <Size> K_pos_a;
        vector <Size> K_pos_b;
        for (Size i = 0; i < a->second.size()-1; ++i)
        {
          if (a->second.toUnmodifiedString()[i] == 'K') K_pos_a.push_back(i);
        }
        for (Size i = 0; i < b->second.size()-1; ++i)
        {
          if (b->second.toUnmodifiedString()[i] == 'K') K_pos_b.push_back(i);
        }
*/
        // mass peptide1 + mass peptide2 + cross linker mass - cross link loss
        double cross_link_mass = a->second.getMonoWeight() + b->second.getMonoWeight() + cross_link_mass_light - cross_link_mass_loss_type2;
        mass_to_candidates.insert(make_pair(cross_link_mass, make_pair<AASequence, AASequence>(a->second, b->second)));
      }
    }

    return mass_to_candidates;
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

    cout << "Smallest double: " << std::numeric_limits<double>::min() << "Double closest to 1: " <<  nexttoward(1.0, 0.0) <<endl;

    return EXECUTION_OK;
    //##########################
    */


    const string in_mzml(getStringOption_("in"));
    const string in_fasta(getStringOption_("database"));
    const string in_consensus(getStringOption_("consensus"));
    const string out_idxml(getStringOption_("out"));

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

    double cross_link_mass_light = getDoubleOption_("cross_linker:mass_light");
    double cross_link_mass_heavy = getDoubleOption_("cross_linker:mass_heavy");
    double cross_link_mass_loss_type2 = getDoubleOption_("cross_linker:mass_loss_type2");

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
    PeakMap spectra;
    MzMLFile f;
    f.setLogType(log_type_);

    PeakFileOptions options;
    options.clearMSLevels();
    options.addMSLevel(2);
    f.getOptions() = options;
    f.load(in_mzml, spectra);
    spectra.sortSpectra(true);

    // filter noise
    progresslogger.startProgress(0, 1, "Filtering spectra...");
    preprocessSpectra_(spectra);
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

    progresslogger.startProgress(0, static_cast<Size>(fasta_db.end() - fasta_db.begin()), "Scoring peptide models against spectra...");

    // lookup for processed peptides. must be defined outside of omp section and synchronized
    multimap<StringView, AASequence> processed_peptides;
    
    // set minimum size of peptide after digestion
    Size min_peptide_length = getIntOption_("peptide:min_size");

    // build multimap of precursor mass to scan index
    multimap<double, Size> multimap_mass_2_scan_index;

    vector<PeptideIdentification> pseudo_ids; // used to map precursor positions to consensus features
    for (PeakMap::ConstIterator s_it = spectra.begin(); s_it != spectra.end(); ++s_it)
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

    vector<ProteinIdentification> protein_ids;
    // protein identifications (leave as is...)
    protein_ids = vector<ProteinIdentification>(1);
    protein_ids[0].setDateTime(DateTime::now());
    protein_ids[0].setSearchEngineVersion(VersionInfo::getVersion());
    idmapper.annotate(cfeatures, pseudo_ids, protein_ids, true, true);

    // maps the index of a light precursor peptide to its corresponding heavier partner
    map<Size, Size> map_light_to_heavy;
    for (ConsensusMap::const_iterator cit = cfeatures.begin(); cit != cfeatures.end(); ++cit)
    {
      if (cit->getFeatures().size() == 2 && cit->getPeptideIdentifications().size() == 2)
      {
        const PeptideIdentification& pi_0 = cit->getPeptideIdentifications()[0];
        const PeptideIdentification& pi_1 = cit->getPeptideIdentifications()[1];
        map_light_to_heavy[pi_0.getMetaValue("scan_index")] = pi_1.getMetaValue("scan_index");
      }
    }
   
    //cout << "Number of MS2 pairs connceted by consensus feature: " << map_light_to_heavy.size() << endl;

    // create common peak / shifted peak spectra for all pairs
    PreprocessedPairSpectra_ preprocessed_pair_spectra = preprocessPairs_(spectra, map_light_to_heavy, ms2_aligner, cross_link_mass_light, cross_link_mass_heavy);
 
    Size count_proteins = 0;
    Size count_peptides = 0;
    
#ifdef _OPENMP
#pragma omp parallel for
#endif
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

    // TODO Initialize output to file
    String out_file = getStringOption_("out");
    ofstream xml_file;
    xml_file.open(out_file.c_str(), ios::trunc);
    xml_file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
    xml_file << "<?xml-stylesheet type=\"text/xsl\" href=\"\"?>" << endl;
    xml_file << "<xquest_results xquest_version=\"xquest 2.1.1\" date=\"Fri Dec 18 12:28:23 2015\" author=\"Thomas Walzthoeni,Oliver Rinner\" homepage=\"http://proteomics.ethz.ch\" deffile=\"xquest.def\" tolerancemeasure_ms2=\"Da\" crosslinkername=\"DSS\" commonlossxcorrweigth=\"3\" poolisotopes=\"0\" xcorr_tolerance_window=\"0\" redundant_peps=\"0\" picktolerance=\"500\" monolinkmw=\"156.0786442,155.0964278\" search_maxcandidate_peps=\"250\" requiredmissed_cleavages=\"0\" picktolerance_measure=\"ppm\" database=\"/home/eugen/MSData/26S_testdataset/db/26Syeast.fasta\" fragmentresiduals=\"HASH(0x37998e0)\" xlinktypes=\"1111\" maxiontaghits=\"0\" AArequired=\"K\" miniontaghits=\"1\" cp_minpeaknumber=\"25\" cp_isotopediff=\"12.075321\" xkinkerID=\"DSS\" maxdigestlength=\"50\" y_ion=\"19.0183888\" cp_nhighest=\"100\" uselossionsformatching=\"0\" mindigestlength=\"5\" indexcharges_common=\"ARRAY(0x378c600)\" enzyme_num=\"1\" testionspick=\"intensity\" xlink_ms2tolerance=\"0.3\" waterloss=\"0\" database_dc=\"/home/eugen/MSData/26S_testdataset/db/26Syeast_decoy.fasta\" iontag_match_xlinkions=\"1\" averageMS2=\"0\" wTICweight=\"12.829\" minionsize=\"200\" x_ion=\"44.9976546\" commonxcorrweigth=\"10\" a_ion=\"-26.9870904\" drawlogscale=\"0\" fwd_ions=\"a|b|c\" experiment=\"BSAtest\" matchoddsweight=\"1.973\" nh3loss=\"0\" verbose=\"0\" enumeration_index_mode=\"smarthash\" maxionsize=\"2000\" outputpath=\"aleitner_M1012_004_matched\" search_intercrosslinks=\"1\" cp_tolerancemeasure=\"ppm\" reportnbesthits=\"5\" cp_dynamic_range=\"1000\" xlinkermw=\"138.0680796\" ms1tol_maxborder=\"10\" enumerate=\"0\" ionindexintprecision=\"10\" writetodiskaftern=\"100000000\" intsumweight=\"0.018\" xlinkxcorrweigth=\"10\" search_monolinks=\"1\" cp_peakratio=\"0.3\" rev_ions=\"x|y|z\" xcorrprecision=\"0.2\" RuntimeDecoys=\"1\" intprecision=\"10\" Iontag_charges_for_index=\"1\" Iontag_writeaftern=\"1200\" z_ion=\"2.9998388\" tryptic_termini=\"2\" cp_tolerance=\"400\" printpeptides=\"1\" ntestions=\"100\" b_ion=\"1.0078246\" normxcorr=\"1\" ioncharge_xlink=\"ARRAY(0x378c8d0)\" ionseries_array=\"ARRAY(0x378cd68)\" cp_scaleby=\"max\" printdigestpeps=\"0\" cp_tolerancexl=\"500\" search_intralinks=\"1\" minionintensity=\"1\" nvariable_mod=\"1\" missed_cleavages=\"2\" ntermxlinkable=\"0\" cp_threshold=\"1\" tolerancemeasure=\"ppm\" CID_match2ndisotope=\"1\" variable_mod=\"M,15.99491\" nocutatxlink=\"1\" minpepmr=\"550\" realintensities4xcorr=\"0\" printcandidatepeps=\"0\" xcorrdelay=\"5\" xcorrxweight=\"2.488\" minhits=\"1\" ms1tol_minborder=\"-10\" drawspectra=\"0\" cp_scaleintensity=\"1\" ms2tolerance=\"0.2\" maxpepmr=\"5500\" printtables=\"0\" ionseries=\"HASH(0x378cc00)\" usenprescores=\"100\" search_intracrosslinks=\"1\" Iontagmode=\"1\" ms1tolerance=\"10\" ioncharge_common=\"ARRAY(0x378c720)\" xcorrbweight=\"21.279\" c_ion=\"18.0343724\" copydb2resdir=\"1\" Hatom=\"1.007825032\" >" << endl;

    // create spectrum generator
    TheoreticalSpectrumGenerator spectrum_generator;
    TheoreticalSpectrumGeneratorXLinks specGen;

    // TODO constant binsize for HashGrid and aPrioriProb computation
    double tolerance_binsize = 0.2;

    cout << "Peptide candidates: " << processed_peptides.size() << endl;

    //TODO refactor, so that only the used mode is initialized and the pre-scoring code only appears once
    // Initialize enumeration mode
    multimap<double, pair<AASequence, AASequence> > enumerated_cross_link_masses;

    //TODO remove, adapt to ppm
    HashGrid1D hg(0.0, 20000.0, tolerance_binsize);

    if (!ion_index_mode)
    {
      enumerated_cross_link_masses = enumerateCrossLinksAndMasses_(processed_peptides, cross_link_mass_light, cross_link_mass_loss_type2);
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
        // TODO Spectra should be: a until last K, b from first K
        specGen.getCommonIonSpectrum(theo_spectrum, *seq, 2); // TODO check which charge and which ion series are used for ion index
      
        //sort by mz
        theo_spectrum.sortByPosition();

        for (Size i = 0; i != theo_spectrum.size(); ++i)
        {
          hg.insert(theo_spectrum[i].getMZ(), seq);
        }
      }
      cout << " finished."  << endl;
    }

    // TODO test variable, can be removed
    double pScoreMax = 0;
    double TICMax = 0;
    double wTICMax = 0;
    double intsumMax = 0;
    double matchOddsMax = 0;
    double xcorrxMax = 0;
    double xcorrcMax = 0;

    // iterate over all spectra
    for (SignedSize scan_index = 0; scan_index < static_cast<SignedSize>(spectra.size()); ++scan_index)
    {
      // assume that current MS2 corresponds to a light spectrum
      const PeakSpectrum& spectrum_light = spectra[scan_index];
      const double precursor_charge = spectrum_light.getPrecursors()[0].getCharge();
      const double precursor_mz = spectrum_light.getPrecursors()[0].getMZ();
      const double precursor_mass = precursor_mz * static_cast<double>(precursor_charge) - static_cast<double>(precursor_charge) * Constants::PROTON_MASS_U;

       // TODO print information about new peak to file (starts with <spectrum_search..., ends with </spectrum_search>
       // TODO what to do with useless information? leave =0 or delete?
       // Mean ion intensity (light spectrum, TODO add heavy spectrum)
       double mean_intensity= 0;
       for (SignedSize j = 0; j < static_cast<SignedSize>(spectrum_light.size()); ++j) mean_intensity += spectrum_light[j].getIntensity();
       mean_intensity = mean_intensity / spectrum_light.size();

       xml_file << "<spectrum_search spectrum=\"" << spectrum_light.getName() << "\" mean_ionintensity=\"" << mean_intensity << "\" ionintensity_stdev=\"" << "TODO" << "\" addedMass=\"" << "TODO" << "\" iontag_ncandidates=\"" << "TODO"
         << "\"  apriori_pmatch_common=\"" << "TODO" << "\" apriori_pmatch_xlink=\"" << "TODO" << "\" ncommonions=\"" << "TODO" << "\" nxlinkions=\"" << "TODO" << "\" mz_precursor=\"" << precursor_mz
         << "\" scantype=\"" << "light_heavy" << "\" charge_precursor=\"" << precursor_charge << "\" Mr_precursor=\"" << precursor_mass <<  "\" rtsecscans=\"" << "TODO" << "\" mzscans=\"" << "0:1" << "\" >" << endl;


      // map light to heavy
      map<Size, Size>::const_iterator scan_index_light_it = map_light_to_heavy.find(scan_index);
   
      if (scan_index_light_it != map_light_to_heavy.end())
      { 
        // cout << "Pair: " << scan_index << ", " << scan_index_heavy << " (mz, charge, mass) " << precursor_mz << "," << precursor_charge << "," << precursor_mass <<  endl;

        const PeakSpectrum& common_peaks = preprocessed_pair_spectra.spectra_common_peaks[scan_index];

        std::vector< double > peak_preScore;
        std::vector< double > peak_wTIC;
        std::vector< double > peak_intSum;
        std::vector< double > peak_matchOdds;
        std::vector< std::vector< double > > peak_xcorrx;
        std::vector< double > peak_xcorrx_max;
        std::vector< std::vector< double > > peak_xcorrc;
        std::vector< double > peak_xcorrc_max;
        std::vector< double > peak_score;
        std::vector< TheoreticalSpectrumGeneratorXLinks::ProteinProteinCrossLink > peak_candidate_data;
        std::vector< Size > peak_matched_spec_common_alpha;
        std::vector< Size > peak_matched_spec_common_beta;
        std::vector< Size > peak_matched_spec_xlink_alpha;
        std::vector< Size > peak_matched_spec_xlink_beta;

        if(common_peaks.size() > 0 || preprocessed_pair_spectra.spectra_xlink_peaks[scan_index].size() > 0) // TODO: check if this is done in xQuest
        {
          // determine candidates
          vector<pair<AASequence, AASequence> > candidates;
          if (ion_index_mode)
          {
            // TODO Use 50 most intense common peaks of exp. spectrum, consider all peptides that produce any of these as theor. common ions
            NLargest nlargest_filter = NLargest(50);
            PeakSpectrum common_peaks_50 = common_peaks;
            nlargest_filter.filterPeakSpectrum(common_peaks_50);

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
            std::sort(ion_tag_candidates.begin(), ion_tag_candidates.end());
            vector<AASequence*>::iterator last_unique = std::unique(ion_tag_candidates.begin(), ion_tag_candidates.end());
            ion_tag_candidates.erase(last_unique, ion_tag_candidates.end());
            cout << "Ion tag candidates before mass filtering: " << ion_tag_candidates.size() << endl;

            //vector<pair<AASequence, AASequence> > candidates;
            for (Size i = 0; i != ion_tag_candidates.size(); ++i)
            {
              const AASequence* peptide_a = ion_tag_candidates[i];
              
              for (Size j = i + 1; j < ion_tag_candidates.size(); ++j)
              {
                const AASequence* peptide_b = ion_tag_candidates[j];

                double cross_link_mass = peptide_a->getMonoWeight() + peptide_b->getMonoWeight() + cross_link_mass_light - cross_link_mass_loss_type2; //TODO: find a way to precalculate individual peptide masses
                double error_Da = abs(cross_link_mass - precursor_mass);
                if (error_Da < precursor_mass_tolerance)
                {
                  candidates.push_back(make_pair(*peptide_a, *peptide_b));
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
            multimap<double, pair<AASequence, AASequence> >::const_iterator low_it;
            multimap<double, pair<AASequence, AASequence> >::const_iterator up_it;

            if (precursor_mass_tolerance_unit_ppm) // ppm
            {
              low_it = enumerated_cross_link_masses.lower_bound(precursor_mass - precursor_mass * precursor_mass_tolerance * 1e-6);
              up_it = enumerated_cross_link_masses.upper_bound(precursor_mass + precursor_mass * precursor_mass_tolerance * 1e-6);
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

//            for (Size i = 0; i != candidates.size(); ++i)
//            {
//              const pair<AASequence, AASequence>& candidate = filtered_candidates[i];
//              RichPeakSpectrum theoretical_spec;
//              TheoreticalSpectrumGeneratorXLinks::ProteinProteinCrossLink cross_link_candidate;
//              cross_link_candidate.alpha = candidate.first;
//              cross_link_candidate.beta = candidate.second;
//              cross_link_candidate.cross_link_position.first = candidate.first.toUnmodifiedString().find('K');
//              cross_link_candidate.cross_link_position.second = candidate.second.toUnmodifiedString().find('K');

//              specGen.getSpectrum(theoretical_spec, cross_link_candidate, 1);

//              std::vector< std::pair< Size, Size > > matched_fragments_theor_spec;
//              ms2_aligner.getSpectrumAlignment(matched_fragments_theor_spec, theoretical_spec, spectrum_light);

//              // Simplified pre-Score, as Alpha and Beta ions are not yet tracked
//              float pre_score = preScore(matched_fragments_theor_spec.size(), theoretical_spec.size());
//              cout << "Number of matched peaks to theor. spectrum: " << matched_fragments_theor_spec.size() << endl;
//              cout << "Number of theoretical ions: " << theoretical_spec.size() << endl;
//              cout << "Pre Score: " << pre_score << endl;
//              //cout << "Peptide size: " << a->second.size() << "\t" << b->second.size() << "\t" << "K Pos:" << K_pos_a[i] << "\t" << K_pos_b[i] << endl;
//              if (pre_score > pScoreMax) pScoreMax = pre_score;
//            }

//          else // enumeration mode
//          {
//            //cout << "Number of common peaks, xlink peaks: " << preprocessed_pair_spectra.spectra_common_peaks[scan_index].size() << "\t" << preprocessed_pair_spectra.spectra_xlink_peaks[scan_index].size();
//            if (preprocessed_pair_spectra.spectra_common_peaks[scan_index].size() < 1 || preprocessed_pair_spectra.spectra_xlink_peaks[scan_index].size() < 1)
//            {
//              continue;
//            }
//            // determine MS2 precursors that match to the current peptide mass
//            multimap<double, pair<AASequence, AASequence> >::const_iterator low_it;
//            multimap<double, pair<AASequence, AASequence> >::const_iterator up_it;

//            if (precursor_mass_tolerance_unit_ppm) // ppm
//            {
//              low_it = enumerated_cross_link_masses.lower_bound(precursor_mass - precursor_mass * precursor_mass_tolerance * 1e-6);
//              up_it = enumerated_cross_link_masses.upper_bound(precursor_mass + precursor_mass * precursor_mass_tolerance * 1e-6);
//            }
//            else // Dalton
//            {
//              low_it = enumerated_cross_link_masses.lower_bound(precursor_mass - precursor_mass_tolerance);
//              up_it =  enumerated_cross_link_masses.upper_bound(precursor_mass + precursor_mass_tolerance);
//            }

//            if (low_it == up_it) continue; // no matching precursor in data

//          }
            // lists for one peak, to determine best match to the peak
          std::vector< double > candidate_preScore;
          std::vector< double > candidate_wTIC;
          std::vector< double > candidate_intSum;
          std::vector< double > candidate_matchOdds;
          std::vector< std::vector< double > > candidate_xcorrx;
          std::vector< double > candidate_xcorrx_max;
          std::vector< std::vector< double > > candidate_xcorrc;
          std::vector< double > candidate_xcorrc_max;
          std::vector< double > candidate_score;
          std::vector< TheoreticalSpectrumGeneratorXLinks::ProteinProteinCrossLink > candidate_data;
          std::vector< Size > candidate_matched_spec_common_alpha;
          std::vector< Size > candidate_matched_spec_common_beta;
          std::vector< Size > candidate_matched_spec_xlink_alpha;
          std::vector< Size > candidate_matched_spec_xlink_beta;

            // loop over cross-link candidates
//            for (; low_it != up_it; ++low_it)
//            {
//              const pair<AASequence, AASequence>& candidate = low_it->second;


          for (Size i = 0; i != candidates.size(); ++i)
          {
            const pair<AASequence, AASequence>& candidate = candidates[i];

            cout << "Pair: " << candidate.first.toString() << ", " << candidate.second.toString() << " matched to light spectrum " << scan_index << " with m/z: " << spectrum_light.getPrecursors()[0].getMZ() <<  endl;
//          cout << a->second.getMonoWeight() << ", " << b->second.getMonoWeight() << " cross_link_mass: " <<  cross_link_mass <<  endl;

	    RichPeakSpectrum theoretical_spec_beta;
	    RichPeakSpectrum theoretical_spec_alpha;
	    RichPeakSpectrum theoretical_spec_xlinks_alpha;
	    RichPeakSpectrum theoretical_spec_xlinks_beta;

            // Determine larger peptide (alpha)
            bool alpha_first = true;
            if (candidate.second.toUnmodifiedString().size() > candidate.first.toUnmodifiedString().size())
            {
              alpha_first = false;
            } else if (candidate.second.toUnmodifiedString().size() == candidate.first.toUnmodifiedString().size() && candidate.second.getMonoWeight() > candidate.first.getMonoWeight())
            {
              alpha_first = false;
            }

            TheoreticalSpectrumGeneratorXLinks::ProteinProteinCrossLink cross_link_candidate;
            if (alpha_first)
            {
              cross_link_candidate.alpha = candidate.first;
              cross_link_candidate.beta = candidate.second;
              cross_link_candidate.cross_link_position.first = candidate.first.toUnmodifiedString().find('K');
              cross_link_candidate.cross_link_position.second = candidate.second.toUnmodifiedString().find('K');
            } else
            {
              cross_link_candidate.alpha = candidate.second;
              cross_link_candidate.beta = candidate.first;
              cross_link_candidate.cross_link_position.first = candidate.second.toUnmodifiedString().find('K');
              cross_link_candidate.cross_link_position.second = candidate.first.toUnmodifiedString().find('K');
            }

            //specGen.getSpectrum(theoretical_spec, cross_link_candidate, 1);
            //getXLinkIonSpectrum(theoretical_spec_xlinks , cross_link_candidate, 1)
            specGen.getCommonIonSpectrum(theoretical_spec_alpha, cross_link_candidate.alpha, 3);
            specGen.getCommonIonSpectrum(theoretical_spec_beta, cross_link_candidate.beta, 3);
            specGen.getXLinkIonSpectrum(theoretical_spec_xlinks_alpha, theoretical_spec_xlinks_beta, cross_link_candidate, 2, 5);

            std::vector< std::pair< Size, Size > > matched_spec_alpha;
            std::vector< std::pair< Size, Size > > matched_spec_beta;
            std::vector< std::pair< Size, Size > > matched_spec_xlinks_alpha;
            std::vector< std::pair< Size, Size > > matched_spec_xlinks_beta;

            ms2_aligner.getSpectrumAlignment(matched_spec_alpha, theoretical_spec_alpha, preprocessed_pair_spectra.spectra_common_peaks[scan_index]);
            ms2_aligner.getSpectrumAlignment(matched_spec_beta, theoretical_spec_beta, preprocessed_pair_spectra.spectra_common_peaks[scan_index]);
            ms2_aligner_xlinks.getSpectrumAlignment(matched_spec_xlinks_alpha, theoretical_spec_xlinks_alpha, preprocessed_pair_spectra.spectra_xlink_peaks[scan_index]);
            ms2_aligner_xlinks.getSpectrumAlignment(matched_spec_xlinks_beta, theoretical_spec_xlinks_beta, preprocessed_pair_spectra.spectra_xlink_peaks[scan_index]);

              /*
              ms2_aligner.getSpectrumAlignment(matched_spec_alpha, theoretical_spec_alpha, spectrum_light);
              ms2_aligner.getSpectrumAlignment(matched_spec_beta, theoretical_spec_beta, spectrum_light);
              ms2_aligner.getSpectrumAlignment(matched_spec_xlinks_alpha, theoretical_spec_xlinks_alpha, spectrum_light);
              ms2_aligner.getSpectrumAlignment(matched_spec_xlinks_beta, theoretical_spec_xlinks_beta, spectrum_light);
              */

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
                  for (SignedSize j = 0; j < static_cast<SignedSize>(matched_spec_alpha.size()); ++j) matched_current_alpha += preprocessed_pair_spectra.spectra_common_peaks[scan_index][matched_spec_alpha[j].second].getIntensity() * 100;
                  for (SignedSize j = 0; j < static_cast<SignedSize>(matched_spec_beta.size()); ++j) matched_current_beta += preprocessed_pair_spectra.spectra_common_peaks[scan_index][matched_spec_beta[j].second].getIntensity() * 100;
                  for (SignedSize j = 0; j < static_cast<SignedSize>(matched_spec_xlinks_alpha.size()); ++j) matched_current_alpha += preprocessed_pair_spectra.spectra_xlink_peaks[scan_index][matched_spec_xlinks_alpha[j].second].getIntensity() * 100;
                  for (SignedSize j = 0; j < static_cast<SignedSize>(matched_spec_xlinks_beta.size()); ++j) matched_current_beta += preprocessed_pair_spectra.spectra_xlink_peaks[scan_index][matched_spec_xlinks_beta[j].second].getIntensity() * 100;
                  double matched_current = matched_current_alpha + matched_current_beta;

                  double total_current = 0;
                  for (SignedSize j = 0; j < static_cast<SignedSize>(spectrum_light.size()); ++j) total_current += spectrum_light[j].getIntensity() * 100;

                  double TIC = matched_current / total_current;
                  //cout << "matched current: " << matched_current << "\t total_current: " << total_current << endl;
                  //cout << "%TIC: " << TIC << "%";
                  if (TIC > TICMax) TICMax = TIC;

                  double aatotal = candidate.first.size() + candidate.second.size();
                  double invMax = 1 / (std::min(candidate.first.size(), candidate.second.size()) / aatotal);
                  double TIC_weight_alpha = (1 / (candidate.first.size() / aatotal)) / invMax;
                  double TIC_weight_beta = (1 / (candidate.second.size() / aatotal)) / invMax;
                  double wTIC = TIC_weight_alpha * (matched_current_alpha / total_current ) + TIC_weight_beta * (matched_current_beta / total_current);
                  //cout << "\t wTIC: " << wTIC;
                  if (wTIC > wTICMax) wTICMax = wTIC;

                  //cout << "\t Intsum: " << matched_current << endl;
                  if (matched_current > intsumMax) intsumMax = matched_current;

                  // match-odds score
                  double range_c_alpha = theoretical_spec_alpha[theoretical_spec_alpha.size()-1].getMZ() -  theoretical_spec_alpha[0].getMZ();
                  double range_x_alpha= theoretical_spec_xlinks_alpha[theoretical_spec_xlinks_alpha.size()-1].getMZ() -  theoretical_spec_xlinks_alpha[0].getMZ();
                  double range_c_beta = theoretical_spec_beta[theoretical_spec_beta.size()-1].getMZ() -  theoretical_spec_beta[0].getMZ();
                  double range_x_beta = theoretical_spec_xlinks_beta[theoretical_spec_xlinks_beta.size()-1].getMZ() -  theoretical_spec_xlinks_beta[0].getMZ();

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
                  double xcorrprecision = 0.2;
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
                  std::vector< double > xcorrx = xCorrelation(spectrum_light, theoretical_spec_common, 5, xcorrprecision);
                  std::vector< double > xcorrc = xCorrelation(spectrum_light, theoretical_spec_xlinks, 5, xcorrprecision);
                  double xcorrx_max = *std::max_element(xcorrx.begin(), xcorrx.end());
                  double xcorrc_max = *std::max_element(xcorrc.begin(), xcorrc.end());
                  //cout << "XLink Cross correlation score: " << xcorrx_max << "\t Common Cross correlation score: " << xcorrc_max << endl;
                  if (xcorrx_max > xcorrxMax) xcorrxMax = xcorrx_max;
                  if (xcorrc_max > xcorrcMax) xcorrcMax = xcorrc_max;

                  // Compute score from the 4 scores and 4 weights
                  double xcorrx_weight = 2.488;
                  double xcorrc_weight = 21.279;
                  double match_odds_weight = 1.973;
                  double wTIC_weight = 12.829;
                  double intsum_weight = 0.018;

                  double score = xcorrx_weight * xcorrx_max + xcorrc_weight * xcorrc_max + match_odds_weight * match_odds + wTIC_weight * wTIC + intsum_weight * matched_current;

                  candidate_preScore.push_back(pre_score);
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

                  //cout << "Next candidate: " << endl;
              }
            } // candidates for peak finished, determine best matching candidate

            Int top = 0;

            while(!candidate_score.empty() && top < 5)
            {
              top++;

              double max_score = *std::max_element(candidate_score.begin(), candidate_score.end());
              Int max_position = std::distance(candidate_score.begin(), std::max_element(candidate_score.begin(), candidate_score.end()));

              peak_preScore.push_back(candidate_preScore[max_position]);
              peak_wTIC.push_back(candidate_wTIC[max_position]);
              peak_intSum.push_back(candidate_intSum[max_position]);
              peak_matchOdds.push_back(candidate_matchOdds[max_position]);
              peak_xcorrx.push_back(candidate_xcorrx[max_position]);
              peak_xcorrx_max.push_back(candidate_xcorrx_max[max_position]);
              peak_xcorrc.push_back(candidate_xcorrc[max_position]);
              peak_xcorrc_max.push_back(candidate_xcorrc_max[max_position]);
              peak_score.push_back(max_score);
              peak_candidate_data.push_back(candidate_data[max_position]);
              peak_matched_spec_common_alpha.push_back(candidate_matched_spec_common_alpha[max_position]);
              peak_matched_spec_common_beta.push_back(candidate_matched_spec_common_beta[max_position]);
              peak_matched_spec_xlink_alpha.push_back(candidate_matched_spec_xlink_alpha[max_position]);
              peak_matched_spec_xlink_beta.push_back(candidate_matched_spec_xlink_beta[max_position]);

              candidate_preScore.erase(candidate_preScore.begin() + max_position);
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
            }
            cout << "Peak scores: " << peak_score << endl;
            cout << "Peak preScores: " << peak_preScore << endl;
            cout << "Peak wTICs: " << peak_wTIC << endl;
            cout << "Peak xcorrx: " << peak_xcorrx_max << endl;
            cout << "Peak xcorrc: " << peak_xcorrc_max << endl;
            cout << "Peak match-odds: " << peak_matchOdds << endl;
            cout << "Peak Intsum: " << peak_intSum << endl;

            // TODO Write top 5 hits to file
            for (Size i = 0; i < peak_score.size(); ++i)
            {
              String structure = peak_candidate_data[i].alpha.toString() + "-" + peak_candidate_data[i].beta.toString();
              String topology = String("a") +  peak_candidate_data[i].cross_link_position.first + String("-b") + peak_candidate_data[i].cross_link_position.second;
              String id = structure + "-" + topology;
              double weight = peak_candidate_data[i].alpha.getMonoWeight() + peak_candidate_data[i].beta.getMonoWeight() + (cross_link_mass_light - cross_link_mass_loss_type2);

              xml_file << "<search_hit search_hit_rank=\"" << i+1 << "\" id=\"" << id << "\" type=\"xlink\"" << " structure=\"" << structure << "\" seq1=\"" << peak_candidate_data[i].alpha.toString() << "\" seq2=\"" << peak_candidate_data[i].beta.toString()
                << "\" prot1=\"" << "TODO" << "\" prot2=\"" << "TODO" << "\" topology=\"" << topology << "\" xlinkposition=\"" << peak_candidate_data[i].cross_link_position.first << "," << peak_candidate_data[i].cross_link_position.second
                << "\" Mr=\"" << weight << "\" mz=\"" << "TODO" << "\" charge=\"" << "TODO" << "\" xlinkermass=\"" << (cross_link_mass_light - cross_link_mass_loss_type2) << "\" measured_mass=\"" << precursor_mass << "\" error=\"" << "TODO"
                << "\" error_rel=\"" << "TODO" << "\" xlinkions_matched=\"" << (peak_matched_spec_xlink_alpha[i] + peak_matched_spec_xlink_beta[i]) << "\" backboneions_matched=\"" << (peak_matched_spec_common_alpha[i] + peak_matched_spec_common_beta[i])
                << "\" weighted_matchodds_mean=\"" << "TODO" << "\" weighted_matchodds_sum=\"" << "TODO" << "\" match_error_mean=\"" << "TODO" << "\" match_error_stdev=\"" << "TODO" << "\" xcorrx=\"" << peak_xcorrx_max[i] << "\" xcorrb=\"" << peak_xcorrc_max[i] << "\" match_odds=\"" << peak_matchOdds[i] << "\" prescore=\"" << peak_preScore[i]
                << "\" prescore_alpha=\"" << "TODO" << "\" prescore_beta=\"" << "TODO" << "\" match_odds_alphacommon=\"" << "TODO" << "\" match_odds_betacommon=\"" << "TODO" << "\" match_odds_alphaxlink=\"" << "TODO"
                << "\" match_odds_betaxlink=\"" << "TODO" << "\" num_of_matched_ions_alpha=\"" << (peak_matched_spec_common_alpha[i] + peak_matched_spec_xlink_alpha[i]) << "\" num_of_matched_ions_beta=\"" << (peak_matched_spec_common_beta[i] + peak_matched_spec_xlink_beta[i]) << "\" num_of_matched_common_ions_alpha=\"" << peak_matched_spec_common_alpha[i]
                << "\" num_of_matched_common_ions_beta=\"" << peak_matched_spec_common_beta[i] << "\" num_of_matched_xlink_ions_alpha=\"" << peak_matched_spec_xlink_alpha[i] << "\" num_of_matched_xlink_ions_beta=\"" << peak_matched_spec_xlink_beta[i] << "\" xcorrall=\"" << "TODO" << "\" TIC=\"" << "TODO"
                << "\" TIC_alpha=\"" << "TODO" << "\" TIC_beta=\"" << "TODO" << "\" wTIC=\"" << peak_wTIC[i] << "\" intsum=\"" << peak_intSum[i] << "\" apriori_match_probs=\"" << "TODO" << "\" apriori_match_probs_log=\"" << "TODO"
                << "\" series_score_mean=\"" << "TODO" << "\" annotated_spec=\"" << "" << "\" score=\"" << peak_score[i] << "\" >" << endl;
              xml_file << "</search_hit>" << endl;
            }


            cout << "Next Peak ################################## \n";
          }

        }
      xml_file << "</spectrum_search>" << endl;
    }
    xml_file << "</xquest_results>" << endl;
    xml_file.close();
    cout << "Pre Score maximum: " << pScoreMax << "\t TIC maximum: " << TICMax << "\t wTIC maximum: " << wTICMax << "\t Match-Odds maximum: " << matchOddsMax << endl;
    cout << "XLink Cross-correlation maximum: " << xcorrxMax << "\t Common Cross-correlation maximum: " << xcorrcMax << "\t Intsum maximum: " << intsumMax << endl;
 
    return EXECUTION_OK;
  }


};

int main(int argc, const char** argv)
{
/*
    in_strings.push_back("GUA1372-S14-A-LRRK2_DSS_1A3.03873.03873.3.dta,GUA1372-S14-A-LRRK2_DSS_1A3.03863.03863.3.dta\n712.0402832\n3\n202.067596435547\t1.00478214074753\t0\n225.138473510742\t10.5088150882632\t0\n240.306930541992\t1.97429014972974\t0\n246.150848388672\t2.11073807368531\t0\n253.116088867188\t0.788884157284357\t0\n254.077026367188\t1.28374306224133\t0\n302.334899902344\t1.01611426130302\t0\n363.247711181641\t1.21699536030682\t0\n364.359741210938\t3.07224448353778\t0\n365.233764648438\t1.41221257514478\t0\n382.492919921875\t1.44150675140929\t0\n454.173889160156\t1.55469310937435\t0\n468.17431640625\t1.25361597022867\t0\n476.370422363281\t1.18878053554465\t0\n478.18017578125\t1.37806417103704\t0\n503.177612304688\t1.00284247589148\t0\n505.511840820312\t5.09655338578696\t0\n508.137390136719\t2.11522771360041\t0\n511.180541992188\t2.15057904382722\t0\n571.420043945312\t6.90119656302518\t0\n583.662292480469\t4.66461443760661\t0\n587.343444824219\t2.87846978117167\t0\n589.50537109375\t7.13272916133555\t0\n594.251281738281\t4.93514552968795\t0\n595.426513671875\t3.44382536510506\t0\n596.37744140625\t11.9700389831397\t0\n610.205810546875\t4.63570482669097\t0\n661.512268066406\t2.14856175242209\t0\n666.367126464844\t4.27393999349033\t0\n669.846130371094\t3.10418295682631\t0\n673.281494140625\t1.31777402765441\t0\n675.322631835938\t3.47432917419416\t0\n716.352233886719\t8.67972468511475\t0\n725.564880371094\t12.897237882998\t0\n727.754760742188\t2.45193672303815\t0\n773.124877929688\t14.9709238344739\t0\n783.431457519531\t1.8729551112151\t0\n790.326721191406\t0.925976281710217\t0\n807.462829589844\t2.61764400546783\t0\n814.522888183594\t2.61683748845419\t0\n838.480285644531\t5.61032298854543\t0\n899.440734863281\t1.11527454352849\t0\n935.585876464844\t5.27870693712852\t0\n936.455749511719\t6.41740507779915\t0\n953.285339355469\t9.95025553746705\t0\n954.578002929688\t4.59257243378958\t0\n963.523315429688\t1.12443797236897\t0\n995.837219238281\t18.0752277973914\t0\n1025.43994140625\t3.67398836852921\t0\n1032.45239257812\t1.3500013173682\t0\n1050.49426269531\t19.376810410944\t0\n1051.4091796875\t5.43402852930932\t0\n1068.53137207031\t14.0416735388994\t0\n1069.45922851562\t8.71293686167004\t0\n1093.45495605469\t0.981511228178035\t0\n1105.52648925781\t1.47770612314566\t0\n1111.71630859375\t1.59983951167604\t0\n1226.41040039062\t1.80091881523256\t0\n1323.61682128906\t2.26763131809332\t0\n1520.763671875\t1.19462550075974\t0\n");
    cout << out << endl;
  */  
/*
  RichPeakSpectrum s;
  vector<Precursor> precursors;
  Precursor pc;
  pc.setMZ(712.0402832);
  pc.setCharge(3);
  precursors.push_back(pc);
  s.setPrecursors(precursors);
  String header = "GUA1372-S14-A-LRRK2_DSS_1A3.03873.03873.3.dta,GUA1372-S14-A-LRRK2_DSS_1A3.03863.03863.3.dta";  // needs to be provided for common and xlinker spectra
  RichPeak1D p;
  p.setMZ(202.067596435547);  
  p.setIntensity(1.00478214074753);
  p.setMetaValue("z", 0); // only set this metavalue for common and xlinker spectra (not for heavy or light)
  s.push_back(p);
  p.setMZ(225.138473510742);  
  p.setIntensity(10.5088150882632);
  p.setMetaValue("z", 0);
  s.push_back(p);
  p.setMZ(240.306930541992);  
  p.setIntensity(1.97429014972974);
  p.setMetaValue("z", 0);
  s.push_back(p);
  cout << TOPPxQuest::getxQuestBase64EncodedSpectrum_(s, header) << endl;
*/
  TOPPxQuest tool;
  
  return tool.main(argc, argv);
}

