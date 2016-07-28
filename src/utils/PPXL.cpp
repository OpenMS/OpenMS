
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
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/XLMS/OpenXQuestScores.h>
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

#ifdef _OPENMP
#include <omp.h>
#define NUMBER_OF_THREADS (omp_get_num_threads())
#else
#define NUMBER_OF_THREADS (1)
#endif

using namespace std;
using namespace OpenMS;

/**
    @page UTILS_xQuest xQuest

    @brief Perform protein-protein cross-linking experiment search.

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
    registerDoubleOption_("cross_linker:mass", "<mass>", 138.0680796, "Mass of the light cross-linker, linking two residues on one or two peptides", false);
    //registerDoubleOption_("cross_linker:mass_iso_shift", "<mass>", 12.075321, "Mass of the isotopic shift between the light and heavy linkers", false);
    registerDoubleList_("cross_linker:mass_mono_link", "<mass>", ListUtils::create<double>("156.0786442, 155.0964278"), "Possible masses of the linker, when attached to only one peptide", false);
    registerStringList_("cross_linker:names", "<list of strings>", ListUtils::create<String>("Xlink:DSS, Xlink:DSS!Hydrolyzed, Xlink:DSS!Amidated"), "Names of the searched cross-links, first the cross-link and then the mono-links in the same order as their masses", false);

    registerTOPPSubsection_("algorithm", "Algorithm Options");
    registerStringOption_("algorithm:candidate_search", "<param>", "index", "Mode used to generate candidate peptides.", false, false);
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
  void preprocessSpectra_(PeakMap& exp)
  {
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
    WindowMower window_mower_filter;
    Param filter_param = window_mower_filter.getParameters();
    filter_param.setValue("windowsize", 100.0, "The size of the sliding window along the m/z axis.");
    filter_param.setValue("peakcount", 20, "The number of peaks that should be kept.");
    filter_param.setValue("movetype", "jump", "Whether sliding window (one peak steps or jumping window window size steps) should be used.");
    window_mower_filter.setParameters(filter_param);
    NLargest nlargest_filter = NLargest(500);   // De-noising in xQuest: Dynamic range = 1000, 250 most intense peaks?

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (Size exp_index = 0; exp_index < exp.size(); ++exp_index)
    {
//      // sort by mz
      exp[exp_index].sortByPosition();
      window_mower_filter.filterPeakSpectrum(exp[exp_index]);
      nlargest_filter.filterSpectrum(exp[exp_index]);
//      // sort (nlargest changes order)
//      exp[exp_index].sortByPosition();
     }
   }

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

//  // Enumerates all possible combinations containing a cross-link, without specific cross-link positions. (There are cases where multiple positions are possible, but they have the same precursor mass)
//  // At this point the only difference between mono-links and loop-links is the added cross-link mass
//  static multimap<double, pair<const AASequence*, const AASequence*> > enumerateCrossLinksAndMasses_(const multimap<StringView, AASequence>&  peptides, double cross_link_mass, const DoubleList& cross_link_mass_mono_link, const StringList& cross_link_residue1, const StringList& cross_link_residue2)
//  {
//    multimap<double, pair<const AASequence*, const AASequence*> > mass_to_candidates;
//    Size countA = 0;

//    vector<const StringView*> peptide_SVs;
//    vector<const AASequence*> peptide_AASeqs;
//    // preparing vectors compatible with openmp multi-threading, TODO this should only be a temporary fix (too much overhead?)
//    for (map<StringView, AASequence>::const_iterator a = peptides.begin(); a != peptides.end(); ++a)
//    {
//      peptide_SVs.push_back(&(a->first));
//      peptide_AASeqs.push_back(&(a->second));
//    }

//    //for (map<StringView, AASequence>::const_iterator a = peptides.begin(); a != peptides.end(); ++a) // old loop version

//#ifdef _OPENMP
//#pragma omp parallel for schedule(guided)
//#endif
//    for (Size p1 = 0; p1 < peptide_AASeqs.size(); ++p1)
//    {
//      String seq_first = peptide_AASeqs[p1]->toUnmodifiedString();

//      countA += 1;
//      if (countA % 50 == 0)
//      {
//        LOG_DEBUG << "Enumerating pairs with sequence " << countA << " of " << peptides.size() << ";\t Current pair count: " << mass_to_candidates.size() << endl;
//      }

//      // generate mono-links
//      for (Size i = 0; i < cross_link_mass_mono_link.size(); i++)
//      {
//        double cross_linked_pair_mass = peptide_AASeqs[p1]->getMonoWeight() + cross_link_mass_mono_link[i];
//        // Make sure it is clear this is a monolink, (is a NULL pointer a good idea?)
//#ifdef _OPENMP
//#pragma omp critical
//#endif
//        mass_to_candidates.insert(make_pair(cross_linked_pair_mass, make_pair<const AASequence*, const AASequence*>(peptide_AASeqs[p1], NULL)));
//      }

//      // generate loop-links
//      bool first_res = false;
//      bool second_res = false;
//      for (Size k = 0; k < seq_first.size()-1; ++k)
//      {
//        for (Size i = 0; i < cross_link_residue1.size(); ++i)
//        {
//          if (seq_first.substr(k, 1) == cross_link_residue1[i])
//          {
//            first_res = true;
//          }
//        }
//        for (Size i = 0; i < cross_link_residue2.size(); ++i)
//        {
//          if (seq_first.substr(k, 1) == cross_link_residue2[i])
//          {
//            second_res = true;
//          }
//        }
//      }
//      // If both sides of a homo- or heterobifunctional cross-linker can link to this peptide, generate the loop-link
//      if (first_res && second_res)
//      {
//        double cross_linked_pair_mass = peptide_AASeqs[p1]->getMonoWeight() + cross_link_mass;
//#ifdef _OPENMP
//#pragma omp critical
//#endif
//        mass_to_candidates.insert(make_pair(cross_linked_pair_mass, make_pair<const AASequence*, const AASequence*>(peptide_AASeqs[p1], NULL)));
//      }

//      // Generate cross-link between two peptides
//      //for (map<StringView, AASequence>::const_iterator b = a; b != peptides.end(); ++b)
//      for (Size p2 = p1; p2 < peptide_AASeqs.size(); ++p2)
//      {
//        // mass peptide1 + mass peptide2 + cross linker mass - cross link loss
//        double cross_linked_pair_mass = peptide_AASeqs[p1]->getMonoWeight() + peptide_AASeqs[p2]->getMonoWeight() + cross_link_mass;
//#ifdef _OPENMP
//#pragma omp critical
//#endif
//        mass_to_candidates.insert(make_pair(cross_linked_pair_mass, make_pair<const AASequence*, const AASequence*>(peptide_AASeqs[p1], peptide_AASeqs[p2])));
//      }
//    }

//    return mass_to_candidates;
//  }

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

    // Deisotopes and single charges spectrum, may be more useful with higher resolution data
    // spectrum must not contain 0 intensity peaks and must be sorted by m/z
    RichPeakSpectrum deisotopeAndSingleChargeMSSpectrum(PeakSpectrum& old_spectrum, Int min_charge, Int max_charge, double fragment_tolerance, bool fragment_unit_ppm, bool keep_only_deisotoped = false, Size min_isopeaks = 3, Size max_isopeaks = 10, bool make_single_charged = false)
    {

      // Input Spectrum originally called "in"
      //PeakSpectrum old_spectrum = in;
      RichPeakSpectrum out;
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

// OLD VERSION, creating deisotoped PeakSpectrum
//      in.clear(false);
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
//            in.push_back(old_spectrum[i]);
//          }
//          else
//          {
//            Peak1D p = old_spectrum[i];
//            p.setMZ(p.getMZ() * z - (z - 1) * Constants::PROTON_MASS_U);
//            in.push_back(p);
//          }
//        }
//        else
//        {
//          // keep all unassigned peaks
//          if (features[i] < 0)
//          {
//            in.push_back(old_spectrum[i]);
//            continue;
//          }

//          // convert mono-isotopic peak with charge assigned by deisotoping
//          if (z != 0)
//          {
//            if (!make_single_charged)
//            {
//              in.push_back(old_spectrum[i]);
//            }
//            else
//            {
//              Peak1D p = old_spectrum[i];
//              p.setMZ(p.getMZ() * z - (z - 1) * Constants::PROTON_MASS_U);
//              in.push_back(p);
//            }
//          }
//        }
//      }

      // NEW VERSION, creating RichPeakSpectrum containing charges
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
            p.setMetaValue("z", z);
            out.push_back(p);
          }
          else
          {
            RichPeak1D p;
            p.setIntensity(old_spectrum[i].getIntensity());
            p.setMZ(old_spectrum[i].getMZ() * z - (z - 1) * Constants::PROTON_MASS_U);
            p.setMetaValue("z", 1);
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
            p.setMetaValue("z", 0);
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
              p.setMetaValue("z", z);
              out.push_back(p);
            }
            else
            {
              RichPeak1D p;
              p.setMZ(old_spectrum[i].getMZ());
              p.setIntensity(old_spectrum[i].getIntensity());
              p.setMetaValue("z", z);
              out.push_back(p);
            }
          }
        }
      }
      out.setPrecursors(old_spectrum.getPrecursors());
      out.setRT(old_spectrum.getRT());
      out.sortByPosition();
      return out;
    }

    // Spectrum Alignment function adapted from SpectrumAlignment.h, intensity_cutoff: 0 for not considering, 0.3 = lower intentity has to be at least 30% of higher intensity (< 70% difference)
    template <typename SpectrumType1, typename SpectrumType2>
    void getSpectrumAlignment(std::vector<std::pair<Size, Size> > & alignment, const SpectrumType1 & s1, SpectrumType2 s2, double tolerance, bool relative_tolerance, double intensity_cutoff = 0.0) const
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

              // check for similar intensity values
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

    // Spectrum Alignment function, that considers intensity similarity.  intensity_cutoff: 0 for not considering, 0.3 = lower intentity has to be at least 30% of higher intensity (< 70% difference)
    void getSpectrumIntensityMatching(std::vector<std::pair<Size, Size> > & alignment, PeakSpectrum s1, PeakSpectrum s2, double tolerance, bool relative_tolerance, double intensity_cutoff = 0.0) const
    {
      if (!s2.isSorted())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Input to SpectrumAlignment is not sorted!");
      }
      PeakSpectrum spec1 = s1;

      // clear result
      alignment.clear();

      s1.sortByIntensity(true);

      //TODO s1 spectrum must be prepared more thoroughly for high resolution ms2. Perl code uses deconvolution / deisotoping / single charging, then adding additional peaks
      // add 4 new peaks for each peak in s1: - C13, -2*C13, + C13, + 2 * C13 with C13shift = 1.0033548378
      // Deisotoping not working on lower res data, like the 26S
      // After that look at matching algorithm in bin/compare_peaks3.pl, line 767



      for (Size i = 0; i != s1.size(); ++i)
      {
        const double& s1_mz = s1[i].getMZ();

        PeakSpectrum peaks = getToleranceWindowPeaks(s2, s1_mz, tolerance, relative_tolerance);
        if (peaks.size() > 1)
        {
          LOG_DEBUG << "SpectrumIntensityMatch: Peaks in tolerance window: " << peaks.size() << endl;
        }

        if (peaks.size() > 0)
        {
          peaks.sortByIntensity();
          double intensity1(s1[i].getIntensity());

          for (Size j = 0; j < peaks.size(); ++j)
          {
            // check for similar intensity values
            double intensity2(peaks[j].getIntensity());
            double high_cutoff = 1 / intensity_cutoff;
            //bool diff_int_clear = (min(intensity1, intensity2) / max(intensity1, intensity2)) >= intensity_cutoff;
            bool diff_int_clear = (intensity1 / intensity2) >= intensity_cutoff && (intensity1 / intensity2) <= high_cutoff;

            // found peak match. if intensity similar enough, update alignment and remove peak, so that it will not get matched again to another peak
            if (diff_int_clear)
            {
              Size s2_index = s2.findNearest(peaks[j].getMZ());
              Size s1_index = spec1.findNearest(s1_mz);
              alignment.push_back(std::make_pair(s1_index, s2_index));
              s2.erase(s2.begin() + s2_index);
              break;
            } else
            {
              // erase the most intense peak, because it does not fulfill the similarity constraint and try with the rest of the peaks
              peaks.erase(peaks.begin() + j);
              --j;
            }
          }
        }
      }
    }

    // Write xQuest.xml output
    void writeXQuestXML(const String& out_file, String base_name, const vector< PeptideIdentification >& peptide_ids, const vector< vector< CrossLinkSpectrumMatch > >& all_top_csms, const PeakMap& spectra)
    {
      String spec_xml_name = base_name + "_matched";
      ofstream xml_file;
      xml_file.open(out_file.c_str(), ios::trunc);
      // XML Header
      xml_file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
      xml_file << "<?xml-stylesheet type=\"text/xsl\" href=\"\"?>" << endl;

      // TODO!!! write actual experiment data
      // original date/time: Fri Dec 18 12:28:23 2015
      DateTime time= DateTime::now();
      String timestring = time.getDate() + " " + time.getTime();

      String precursor_mass_tolerance_unit = getStringOption_("precursor:mass_tolerance_unit");
      String fragment_mass_tolerance_unit = getStringOption_("fragment:mass_tolerance_unit");
      double precursor_mass_tolerance = getDoubleOption_("precursor:mass_tolerance");
      double fragment_mass_tolerance = getDoubleOption_("fragment:mass_tolerance");
      double fragment_mass_tolerance_xlinks = getDoubleOption_("fragment:mass_tolerance_xlinks");

      StringList cross_link_names = getStringList_("cross_linker:names");
      double cross_link_mass = getDoubleOption_("cross_linker:mass");
      DoubleList cross_link_mass_mono_link = getDoubleList_("cross_linker:mass_mono_link");
      String mono_masses;
      for (Size k = 0; k < cross_link_mass_mono_link.size()-1; ++k)
      {
        mono_masses += String(cross_link_mass_mono_link[k]) + ", ";
      }
      mono_masses += cross_link_mass_mono_link[cross_link_mass_mono_link.size()-1];

      const string in_fasta(getStringOption_("database"));
      const string in_decoy_fasta(getStringOption_("decoy_database"));
      StringList cross_link_residue1 = getStringList_("cross_linker:residue1");
      StringList cross_link_residue2 = getStringList_("cross_linker:residue2");
      String aarequired1, aarequired2;
      for (Size k= 0; k < cross_link_residue1.size()-1; ++k)
      {
        aarequired1 += cross_link_residue1[k] + ",";
      }
      aarequired1 += cross_link_residue1[cross_link_residue1.size()-1];
      for (Size k= 0; k < cross_link_residue2.size()-1; ++k)
      {
        aarequired2 += cross_link_residue2[k] + ",";
      }
      aarequired2 += cross_link_residue2[cross_link_residue2.size()-1];

      //double cross_link_mass_iso_shift = getDoubleOption_("cross_linker:mass_iso_shift");
      String enzyme_name = getStringOption_("peptide:enzyme");
      Size missed_cleavages = getIntOption_("peptide:missed_cleavages");

      xml_file << "<xquest_results xquest_version=\"openxquest 1.0\" date=\"" << timestring <<
             "\" author=\"Eugen Netz, Timo Sachsenberg\" tolerancemeasure_ms1=\"" << precursor_mass_tolerance_unit  <<
             "\" tolerancemeasure_ms2=\"" << fragment_mass_tolerance_unit << "\" ms1tolerance=\"" << precursor_mass_tolerance <<
             "\" ms2tolerance=\"" << fragment_mass_tolerance << "\" xlink_ms2tolerance=\"" << fragment_mass_tolerance_xlinks <<
             "\" crosslinkername=\"" << cross_link_names[0] << "\" xlinkermw=\"" << cross_link_mass <<
             "\" monolinkmw=\"" << mono_masses << "\" database=\"" << in_fasta << "\" database_dc=\"" << in_decoy_fasta <<
             "\" xlinktypes=\"1111\" AArequired1=\"" << aarequired1 << "\" AArequired2=\"" << aarequired2 <<  "\" cp_isotopediff=\"" << 0 <<
             "\" enzyme_name=\"" << enzyme_name << "\" outputpath=\"" << spec_xml_name <<
             "\" Iontag_charges_for_index=\"1\" missed_cleavages=\"" << missed_cleavages <<
             "\" ntermxlinkable=\"0\" CID_match2ndisotope=\"1" <<
             "\" variable_mod=\"TODO\" nocutatxlink=\"1\" xcorrdelay=\"5\" >" << endl;



      for (vector< vector< CrossLinkSpectrumMatch > >::const_iterator top_csms_spectrum = all_top_csms.begin(); top_csms_spectrum != all_top_csms.end(); ++top_csms_spectrum)
      {
        vector< CrossLinkSpectrumMatch > top_vector = (*top_csms_spectrum);

        if (!top_vector.empty())
        {
          // Spectrum Data, for each spectrum
          Size scan_index = top_vector[0].scan_index_light;
          //Size scan_index_heavy = top_vector[0].scan_index_heavy;
          const PeakSpectrum& spectrum = spectra[scan_index];
          double precursor_charge = spectrum.getPrecursors()[0].getCharge();

          double precursor_mz = spectrum.getPrecursors()[0].getMZ();
          double precursor_rt = spectrum.getRT();
          double precursor_mass = precursor_mz * static_cast<double>(precursor_charge) - static_cast<double>(precursor_charge) * Constants::PROTON_MASS_U;

          // print information about new peak to file (starts with <spectrum_search..., ends with </spectrum_search>
          String spectrum_light_name = base_name + ".light." + scan_index;
          String spectrum_heavy_name = base_name + ".heavy." + scan_index;

          String spectrum_name = spectrum_light_name + "_" + spectrum_heavy_name;
          String rt_scans = String(precursor_rt) + ":" + String(precursor_rt);
          String mz_scans = String(precursor_mz) + ":" + String(precursor_mz);

          // Mean ion intensity (light spectrum, TODO add heavy spectrum?)
          double mean_intensity= 0;
          for (SignedSize j = 0; j < static_cast<SignedSize>(spectrum.size()); ++j) mean_intensity += spectrum[j].getIntensity();
          mean_intensity = mean_intensity / spectrum.size();

          xml_file << "<spectrum_search spectrum=\"" << spectrum_name << "\" mean_ionintensity=\"" << mean_intensity << "\" ionintensity_stdev=\"" << "TODO" << "\" addedMass=\"" << "TODO" << "\" iontag_ncandidates=\"" << "TODO"
            << "\"  apriori_pmatch_common=\"" << "TODO" << "\" apriori_pmatch_xlink=\"" << "TODO" << "\" ncommonions=\"" << "TODO" << "\" nxlinkions=\"" << "TODO" << "\" mz_precursor=\"" << precursor_mz
            << "\" scantype=\"" << "light" << "\" charge_precursor=\"" << precursor_charge << "\" Mr_precursor=\"" << precursor_mass <<  "\" rtsecscans=\"" << rt_scans << "\" mzscans=\"" << mz_scans << "\" >" << endl;


          for (vector< CrossLinkSpectrumMatch>::const_iterator top_csm = top_csms_spectrum->begin(); top_csm != top_csms_spectrum->end(); ++top_csm)
          {
            String xltype = "monolink";
            String structure = top_csm->cross_link.alpha.toUnmodifiedString();
            String letter_first = structure.substr(top_csm->cross_link.cross_link_position.first, 1);


             // TODO track or otherwise find out, which kind of mono-link it was (if there are several possibilities for the weigths)
            double weight = top_csm->cross_link.alpha.getMonoWeight() + top_csm->cross_link.cross_linker_mass;
//            bool is_monolink = (top_csm->cross_link.cross_link_position.second == -1);
            int alpha_pos = top_csm->cross_link.cross_link_position.first + 1;
            int beta_pos = top_csm->cross_link.cross_link_position.second + 1;

            String topology = String("a") + alpha_pos;
            String id = structure + String("-") + letter_first + alpha_pos + String("-") + static_cast<int>(top_csm->cross_link.cross_linker_mass);

            if (top_csm->cross_link.getType() == TheoreticalSpectrumGeneratorXLinks::ProteinProteinCrossLink::CROSS)
            {
              xltype = "xlink";
              structure += "-" + top_csm->cross_link.beta.toUnmodifiedString();
              topology += String("-b") + beta_pos;
              weight += top_csm->cross_link.beta.getMonoWeight();
              id = structure + "-" + topology;
            }
            else if (top_csm->cross_link.getType() == TheoreticalSpectrumGeneratorXLinks::ProteinProteinCrossLink::LOOP)
            {
              xltype = "intralink";
              topology += String("-b") + beta_pos;
              String letter_second = structure.substr(top_csm->cross_link.cross_link_position.second, 1);
              id = structure + String("-") + letter_first + alpha_pos + String("-") + letter_second + beta_pos;
            }

            // Error calculation
            double cl_mz = (weight + (static_cast<double>(precursor_charge) * Constants::PROTON_MASS_U)) / static_cast<double>(precursor_charge);

            double error = precursor_mz - cl_mz;
            //double rel_error = (error / precursor_mz) / 1e-6;
            double rel_error = (error / cl_mz) / 1e-6;
            //cout << "rel_error[ppm]: " << rel_error << "\terror: " << error << "\tcl_mz: " << cl_mz << "\tweight: " << weight << endl;

            PeptideIdentification pep_id = peptide_ids[top_csm->peptide_id_index];
            vector< PeptideHit > pep_hits = pep_id.getHits();

            String prot_alpha = pep_hits[0].getPeptideEvidences()[0].getProteinAccession();
            if (pep_hits[0].getPeptideEvidences().size() > 1)
            {
              for (Size i = 1; i < pep_hits[0].getPeptideEvidences().size(); ++i)
              {
                prot_alpha = prot_alpha + "," + pep_hits[0].getPeptideEvidences()[i].getProteinAccession();
              }
            }

            String prot_beta = "";

            if (pep_hits.size() > 1)
            {
              prot_beta= pep_hits[1].getPeptideEvidences()[0].getProteinAccession();
              if (pep_hits[1].getPeptideEvidences().size() > 1)
              {
                for (Size i = 1; i < pep_hits[1].getPeptideEvidences().size(); ++i)
                {
                  prot_alpha = prot_alpha + "," + pep_hits[1].getPeptideEvidences()[i].getProteinAccession();
                }
              }
            }
            // Hit Data, for each cross-link to Spectrum Hit (e.g. top 5 per spectrum)
            xml_file << "<search_hit search_hit_rank=\"" <<top_csm->rank << "\" id=\"" << id << "\" type=\"" << xltype << "\" structure=\"" << structure << "\" seq1=\"" << top_csm->cross_link.alpha.toUnmodifiedString() << "\" seq2=\"" << top_csm->cross_link.beta.toUnmodifiedString()
                << "\" prot1=\"" << prot_alpha << "\" prot2=\"" << prot_beta << "\" topology=\"" << topology << "\" xlinkposition=\"" << (top_csm->cross_link.cross_link_position.first+1) << "," << (top_csm->cross_link.cross_link_position.second+1)
                << "\" Mr=\"" << weight << "\" mz=\"" << cl_mz << "\" charge=\"" << precursor_charge << "\" xlinkermass=\"" << top_csm->cross_link.cross_linker_mass << "\" measured_mass=\"" << precursor_mass << "\" error=\"" << error
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

    void writeXQuestXMLSpec(String base_name, const vector< vector< CrossLinkSpectrumMatch > >& all_top_csms, const PeakMap& spectra)
    {
      String spec_xml_filename = base_name + "_matched.spec.xml";
      // XML Header
      ofstream spec_xml_file;
      spec_xml_file.open(spec_xml_filename.c_str(), ios::trunc); // ios::app = append to file, ios::trunc = overwrites file
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
          spec_xml_file << TOPPxQuest::getxQuestBase64EncodedSpectrum_(spectra[i], String(""));
          spec_xml_file << "</spectrum>" << endl;

          spec_xml_file << "<spectrum filename=\"" << spectrum_heavy_name << ".dta" << "\" type=\"heavy\">" << endl;
          spec_xml_file << TOPPxQuest::getxQuestBase64EncodedSpectrum_(spectra[i], String(""));
          spec_xml_file << "</spectrum>" << endl;

          String spectrum_common_name = spectrum_name + String("_common.txt");
          spec_xml_file << "<spectrum filename=\"" << spectrum_common_name << "\" type=\"common\">" << endl;
          spec_xml_file << TOPPxQuest::getxQuestBase64EncodedSpectrum_(spectra[i], spectrum_light_name + ".dta," + spectrum_heavy_name + ".dta");
          spec_xml_file << "</spectrum>" << endl;

          String spectrum_xlink_name = spectrum_name + String("_xlinker.txt");
          spec_xml_file << "<spectrum filename=\"" << spectrum_xlink_name << "\" type=\"xlinker\">" << endl;
          spec_xml_file <<TOPPxQuest::getxQuestBase64EncodedSpectrum_(spectra[i], spectrum_light_name + ".dta," + spectrum_heavy_name + ".dta");
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
    StringList cross_link_names = getStringList_("cross_linker:names");

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
    PeakMap spectra;
    MzMLFile f;
    f.setLogType(log_type_);

    PeakFileOptions options;
    options.clearMSLevels();
    options.addMSLevel(2);
    f.getOptions() = options;
    f.load(in_mzml, spectra);
    spectra.sortSpectra(true);

    // preprocess spectra (filter out 0 values, sort by position)
    progresslogger.startProgress(0, 1, "Filtering spectra...");
    preprocessSpectra_(spectra);
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

    progresslogger.endProgress();

    const Size missed_cleavages = getIntOption_("peptide:missed_cleavages");
    EnzymaticDigestion digestor;
    String enzyme_name = getStringOption_("peptide:enzyme");
    digestor.setEnzyme(enzyme_name);
    digestor.setMissedCleavages(missed_cleavages);

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

    // create common peak / shifted peak spectra for all pairs
    progresslogger.startProgress(0, 1, "Preprocessing Spectra Pairs...");
//    PreprocessedPairSpectra_ preprocessed_pair_spectra = preprocessPairs_(spectra, map_light_to_heavy, cross_link_mass, cross_link_mass_heavy);
    //PreprocessedPairSpectra_ preprocessed_pair_spectra = preprocessPairs_(spectra, spectrum_pairs, cross_link_mass_iso_shift);
    progresslogger.endProgress();

    Size count_proteins = 0;
    Size count_peptides = 0;


    // one identification run
    vector<ProteinIdentification> protein_ids(1);
    protein_ids[0].setDateTime(DateTime::now());
    protein_ids[0].setSearchEngine("OpenMSxQuest");
    protein_ids[0].setSearchEngineVersion(VersionInfo::getVersion());
    protein_ids[0].setPrimaryMSRunPath(spectra.getPrimaryMSRunPath());
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
    search_params.fragment_mass_tolerance_ppm =  fragment_mass_tolerance_unit_ppm ? "ppm" : "Da";
    search_params.precursor_mass_tolerance = precursor_mass_tolerance;
    search_params.precursor_mass_tolerance_ppm = precursor_mass_tolerance_unit_ppm ? "ppm" : "Da";

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

//#ifdef _OPENMP
//#pragma omp parallel for
//#endif
    // digest and filter database
    for (SignedSize fasta_index = 0; fasta_index < static_cast<SignedSize>(fasta_db.size()); ++fasta_index)
    {
//#ifdef _OPENMP
//#pragma omp atomic
//#endif
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
        // skip peptides with invalid AAs // TODO is this necessary?
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
        if (cit->getString().find('K') >= cit->getString().size()-1)
        {
          continue;
        }

//#ifdef _OPENMP
//#pragma omp atomic
//#endif
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

//#ifdef _OPENMP
//#pragma omp critical (processed_peptides_access)
//#endif
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

//#ifdef _OPENMP
//#pragma omp parallel for
//#endif
      for (SignedSize fasta_index = 0; fasta_index < static_cast<SignedSize>(fasta_decoys.size()); ++fasta_index)
      {
//  #ifdef _OPENMP
//  #pragma omp atomic
//  #endif
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
//  #ifdef _OPENMP
//  #pragma omp critical (processed_peptides_access)
//  #endif
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

//  #ifdef _OPENMP
//  #pragma omp atomic
//  #endif
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

//  #ifdef _OPENMP
//  #pragma omp critical (processed_peptides_access)
//  #endif
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

    // create spectrum generator
    TheoreticalSpectrumGenerator spectrum_generator;

    TheoreticalSpectrumGeneratorXLinks specGen;
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
    multimap<double, pair<const AASequence*, const AASequence*> > enumerated_cross_link_masses;

    //TODO remove, adapt to ppm
    HashGrid1D hg(0.0, 20000.0, tolerance_binsize);

    map<AASequence*, MSSpectrum<RichPeak1D> > peptide_spectra;

    if (!ion_index_mode)
    {
      progresslogger.startProgress(0, 1, "Enumerating cross-links...");
      enumerated_cross_link_masses = OpenXQuestScores::enumerateCrossLinksAndMasses_(processed_peptides, cross_link_mass, cross_link_mass_mono_link, cross_link_residue1, cross_link_residue2);
      progresslogger.endProgress();
      LOG_DEBUG << "Enumerated cross-links: " << enumerated_cross_link_masses.size() << endl;
    }
    else
    {
      LOG_DEBUG << "Adding peaks to hash map ...";
      for (map<StringView, AASequence>::iterator a = processed_peptides.begin(); a != processed_peptides.end(); ++a)
      {
        //create theoretical spectrum
        MSSpectrum<RichPeak1D> theo_spectrum = MSSpectrum<RichPeak1D>();
        // LOG_DEBUG << a->second.toString() << endl;
        AASequence * seq = &(a->second);
        // generate common ions
        spectrum_generator.getSpectrum(theo_spectrum, *seq, 3); // TODO check which charge and which ion series are used for ion index
        peptide_spectra.insert(make_pair(seq, theo_spectrum));

        //sort by mz (is done in getCommonIonSpectrum)
//        theo_spectrum.sortByPosition();

        for (Size i = 0; i != theo_spectrum.size(); ++i)
        {
          hg.insert(theo_spectrum[i].getMZ(), seq);
        }
      }
      LOG_DEBUG << " finished."  << endl;
    }

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

// Multithreading options: schedule: static, dynamic, guided with chunk size
#ifdef _OPENMP
#pragma omp parallel for schedule(guided)
#endif
    for (Size scan_index = 0; scan_index < spectra.size(); ++scan_index)
    {

      // If this spectra pair has less than 15 common peaks, then ignore it.
      //TODO is a xquest.def parameter in perl xQuest, set to 25 usually
//      if (preprocessed_pair_spectra.spectra_common_peaks[pair_index].size() < 15)
//      {
//        continue;
//      }

      //Size scan_index = spectrum_pairs[pair_index].first;
      //Size scan_index_heavy = spectrum_pairs[pair_index].second;
      const PeakSpectrum& spectrum= spectra[scan_index];
      LOG_DEBUG << "################## NEW SPECTRUM ##############################" << endl;
      //LOG_DEBUG << "Scan index: " << scan_index << "\tSpectrum native ID: " << spectrum.getNativeID()  << endl;
      cout << "Processing spectrum index: " << scan_index << " / " << spectra.size() << "\tSpectrum native ID: " << spectrum.getNativeID()  << endl;

      const double precursor_charge = spectrum.getPrecursors()[0].getCharge();
      const double precursor_mz = spectrum.getPrecursors()[0].getMZ();
      const double precursor_mass = precursor_mz * static_cast<double>(precursor_charge) - static_cast<double>(precursor_charge) * Constants::PROTON_MASS_U;

      // Mean ion intensity (light spectrum, TODO add heavy spectrum?)
      double mean_intensity= 0;
      for (SignedSize j = 0; j < static_cast<SignedSize>(spectrum.size()); ++j)
      {
        mean_intensity += spectrum[j].getIntensity();
      }
      mean_intensity = mean_intensity / spectrum.size();

      const PeakSpectrum& common_peaks = spectra[scan_index];

      vector< CrossLinkSpectrumMatch > top_csms_spectrum;

      if (common_peaks.size() > 0) // TODO: check if this is done in xQuest?
      {
        // determine candidates
        vector<pair<const AASequence*, const AASequence*> > candidates;
        double allowed_error = 0;

        if (ion_index_mode)
        {
          LOG_DEBUG << "Ion tag Mode, start collecting candidate peptides" << endl;
          // Use 50 most intense common peaks of exp. spectrum, consider all peptides that produce any of these as theor. common ions
          NLargest nlargest_filter = NLargest(50);
          PeakSpectrum common_peaks_50 = common_peaks;
          nlargest_filter.filterSpectrum(common_peaks_50);
          common_peaks_50.sortByPosition();

          set<AASequence*> ion_tag_candidates;
          for (Size i = 0; i != common_peaks_50.size(); ++i)
          {
            const vector<AASequence*> new_ion_tag_candidates = hg.get(common_peaks_50[i].getMZ(), 5000);
            LOG_DEBUG << "Number of candidate peptides for current peak: " << new_ion_tag_candidates.size();
            ion_tag_candidates.insert(new_ion_tag_candidates.begin(), new_ion_tag_candidates.end());
          }
          LOG_DEBUG << "Number of candidate peptides for all 50 peaks: " << ion_tag_candidates.size();
          //LOG_DEBUG << "Ion tag Mode, start uniquifying" << endl;
          //sort(ion_tag_candidates.begin(), ion_tag_candidates.end());
          //vector<AASequence*>::iterator last_unique = unique(ion_tag_candidates.begin(), ion_tag_candidates.end());
          //ion_tag_candidates.erase(last_unique, ion_tag_candidates.end());
          //LOG_DEBUG << "Ion tag Mode, end uniquifying" << endl;

          // Pre-Score all candidates
          LOG_DEBUG << "Start pre-scoring candidate peptides...." << endl;
          vector<pair<double, AASequence*> > pre_scores;
          //for (Size i = 0; i < ion_tag_candidates.size(); ++i)
          for (set<AASequence*>::iterator candidate = ion_tag_candidates.begin(); candidate != ion_tag_candidates.end(); ++candidate)
          {
            vector< pair< Size, Size > > matched_spec;
            getSpectrumAlignment(matched_spec, peptide_spectra.at(*candidate), spectrum, fragment_mass_tolerance, false);

            double pre_score = 0;
            if (matched_spec.size() > 0)
            {
              pre_score = OpenXQuestScores::preScore(matched_spec.size(), peptide_spectra.at(*candidate).size());
            }
            pre_scores.push_back(make_pair(pre_score, *candidate));
          }
          LOG_DEBUG << "Sorting scored results" << endl;
          sort(pre_scores.begin(), pre_scores.end());
          LOG_DEBUG << "Sorting finished" << endl;


          // Clear candidates and add 50 highest scoring
          ion_tag_candidates.clear();
          Size max_candidates = 500;
          for (Size i = 0; i < max_candidates; ++i)
          {
            ion_tag_candidates.insert(pre_scores[i].second);
          }
          LOG_DEBUG << "Pre-scoring completed." << endl;
          LOG_DEBUG << "Ion tag candidates before mass filtering: " << ion_tag_candidates.size() << endl;


          //vector<pair<AASequence, AASequence> > candidates;
          //for (Size i = 0; i != ion_tag_candidates.size(); ++i)
          for (set<AASequence*>::iterator candidateA = ion_tag_candidates.begin(); candidateA != ion_tag_candidates.end(); ++candidateA)
          {
            const AASequence* peptide_a = *candidateA;

            //for (Size j = i + 1; j < ion_tag_candidates.size(); ++j)
            for (set<AASequence*>::iterator candidateB = ion_tag_candidates.begin(); candidateB != ion_tag_candidates.end(); ++candidateB)
            {
              const AASequence* peptide_b = *candidateB;

              double cross_linked_pair_mass = peptide_a->getMonoWeight() + peptide_b->getMonoWeight() + cross_link_mass; //TODO: find a way to precalculate individual peptide masses
              double error_Da = abs(cross_linked_pair_mass- precursor_mass);
              if (error_Da < precursor_mass_tolerance)
              {
                candidates.push_back(make_pair(peptide_a, peptide_b));
              }
            }
          }
            LOG_DEBUG << "Ion tag candidates after mass filtering: " << candidates.size() << endl;
        } else // enumeration mode
        {
          //LOG_DEBUG << "Number of common peaks, xlink peaks: " << preprocessed_pair_spectra.spectra_common_peaks[scan_index].size() << "\t" << preprocessed_pair_spectra.spectra_xlink_peaks[scan_index].size();

          // determine MS2 precursors that match to the current peptide mass
          multimap<double, pair<const AASequence*, const AASequence*> >::const_iterator low_it;
          multimap<double, pair<const AASequence*, const AASequence*> >::const_iterator up_it;


          if (precursor_mass_tolerance_unit_ppm) // ppm
          {
            allowed_error = precursor_mass * precursor_mass_tolerance * 1e-6;
            //low_it = enumerated_cross_link_masses.lower_bound(precursor_mass - allowed_error);
            //up_it = enumerated_cross_link_masses.upper_bound(precursor_mass + allowed_error);
          }
          else // Dalton
          {
            allowed_error = precursor_mass_tolerance;
          }
#ifdef _OPENMP
#pragma omp critical (enumerated_cross_link_masses_access)
#endif
          {
            low_it = enumerated_cross_link_masses.lower_bound(precursor_mass - allowed_error);
            up_it =  enumerated_cross_link_masses.upper_bound(precursor_mass + allowed_error);
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
          vector <SignedSize> link_pos_first;
          vector <SignedSize> link_pos_second;
          AASequence peptide_first = *candidate.first;
          AASequence peptide_second;
          if (candidate.second)
          {
            peptide_second = *candidate.second;
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
          if (candidate.second)
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
              TheoreticalSpectrumGeneratorXLinks::ProteinProteinCrossLink cross_link_candidate;
              // if loop link, and the positions are the same, then it is linking the same residue with itself,  skip this combination, also pos1 > pos2 would be the same link as pos1 < pos2
              if (((seq_second.size() == 0) && (link_pos_first[x] >= link_pos_second[y])) && (link_pos_second[y] != -1))
              {
                continue;
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
                cross_link_candidate.cross_linker_name = cross_link_names[0];
                cross_link_candidates.push_back(cross_link_candidate);
              }
              else
              {
                for (Size k = 0; k < cross_link_mass_mono_link.size(); ++k)
                {
                  // only use the correct mono-links (at this point we know it is a mono-link, but not which one)
                  if (abs(precursor_mass - (peptide_first.getMonoWeight() + cross_link_mass_mono_link[k])) <= allowed_error)
                  {
                    cross_link_candidate.cross_linker_mass = cross_link_mass_mono_link[k];
                    cross_link_candidate.cross_linker_name = cross_link_names[k+1];
                    cross_link_candidates.push_back(cross_link_candidate);
                  }
                }
              }
            }
          }
        }

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
          TheoreticalSpectrumGeneratorXLinks::ProteinProteinCrossLink cross_link_candidate = cross_link_candidates[i];
          double candidate_mz = (cross_link_candidate.alpha.getMonoWeight() + cross_link_candidate.beta.getMonoWeight() +  cross_link_candidate.cross_linker_mass+ (static_cast<double>(precursor_charge) * Constants::PROTON_MASS_U)) / precursor_charge;

          LOG_DEBUG << "Pair: " << cross_link_candidate.alpha.toString() << "-" << cross_link_candidate.beta.toString() << " matched to light spectrum " << scan_index << "\t and heavy spectrum " << scan_index
              << " with m/z: " << precursor_mz << "\t" << "and candidate m/z: " << candidate_mz << "\tK Positions: " << cross_link_candidate.cross_link_position.first << "\t" << cross_link_candidate.cross_link_position.second << endl;
//          LOG_DEBUG << a->second.getMonoWeight() << ", " << b->second.getMonoWeight() << " cross_link_mass: " <<  cross_link_mass <<  endl;

	  CrossLinkSpectrumMatch csm;
	  csm.cross_link = cross_link_candidate;

	  RichPeakSpectrum theoretical_spec_common_alpha;
	  RichPeakSpectrum theoretical_spec_common_beta;
	  RichPeakSpectrum theoretical_spec_xlinks_alpha;
	  RichPeakSpectrum theoretical_spec_xlinks_beta;

	  bool type_is_cross_link = cross_link_candidate.getType() == TheoreticalSpectrumGeneratorXLinks::ProteinProteinCrossLink::CROSS;

          specGen.getCommonIonSpectrum(theoretical_spec_common_alpha, cross_link_candidate, 2, true);
          if (type_is_cross_link)
          {
            specGen.getCommonIonSpectrum(theoretical_spec_common_beta, cross_link_candidate, 2, false);
            specGen.getXLinkIonSpectrum(theoretical_spec_xlinks_alpha, theoretical_spec_xlinks_beta, cross_link_candidate, 2, precursor_charge);
          } else
          {
            // Function for mono-links or loop-links
            specGen.getXLinkIonSpectrum(theoretical_spec_xlinks_alpha, cross_link_candidate, 2, precursor_charge);
          }

          vector< pair< Size, Size > > matched_spec_common_alpha;
          vector< pair< Size, Size > > matched_spec_common_beta;
          vector< pair< Size, Size > > matched_spec_xlinks_alpha;
          vector< pair< Size, Size > > matched_spec_xlinks_beta;

          if (spectrum.size() > 0)
          {
            getSpectrumAlignment(matched_spec_common_alpha, theoretical_spec_common_alpha, spectrum, fragment_mass_tolerance, false);
            getSpectrumAlignment(matched_spec_common_beta, theoretical_spec_common_beta, spectrum, fragment_mass_tolerance, false);
          }
          if (spectrum.size() > 0)
          {
            getSpectrumAlignment(matched_spec_xlinks_alpha, theoretical_spec_xlinks_alpha, spectrum, fragment_mass_tolerance_xlinks, false);
            getSpectrumAlignment(matched_spec_xlinks_beta, theoretical_spec_xlinks_beta, spectrum, fragment_mass_tolerance_xlinks, false);
          }

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
               pre_score = OpenXQuestScores::preScore(matched_alpha_count, theor_alpha_count, matched_beta_count, theor_beta_count);
             }
             else
             {
              pre_score = OpenXQuestScores::preScore(matched_alpha_count, theor_alpha_count);
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
             double intsum = OpenXQuestScores::total_matched_current(matched_spec_common_alpha, matched_spec_common_beta, matched_spec_xlinks_alpha, matched_spec_xlinks_beta, spectrum, spectrum);


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
              double intsum_alpha = OpenXQuestScores::matched_current_chain(matched_spec_common_alpha, matched_spec_xlinks_alpha, spectrum, spectrum);
              double intsum_beta = 0;
              if (type_is_cross_link)
              {
                intsum_beta = OpenXQuestScores::matched_current_chain(matched_spec_common_beta, matched_spec_xlinks_beta, spectrum, spectrum);
              }

              // normalize TIC_alpha and  _beta
              if ((intsum_alpha + intsum_beta) > 0.0)
              {
                intsum_alpha = intsum_alpha * intsum / (intsum_alpha + intsum_beta);
                intsum_beta = intsum_beta *  intsum / (intsum_alpha + intsum_beta);
              }

              // compute wTIC
              double wTIC = OpenXQuestScores::weighted_TIC_score(cross_link_candidate.alpha.size(), cross_link_candidate.beta.size(), intsum_alpha, intsum_beta, intsum, total_current, type_is_cross_link);

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
              double match_odds_c_alpha = OpenXQuestScores::match_odds_score(theoretical_spec_common_alpha, matched_spec_common_alpha, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, false);
              double match_odds_x_alpha = OpenXQuestScores::match_odds_score(theoretical_spec_xlinks_alpha, matched_spec_xlinks_alpha, fragment_mass_tolerance_xlinks , fragment_mass_tolerance_unit_ppm, true, n_xlink_charges);
              double match_odds = 0;
              if (type_is_cross_link)
              {
                double match_odds_c_beta = OpenXQuestScores::match_odds_score(theoretical_spec_common_beta, matched_spec_common_beta, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, false);
                double match_odds_x_beta = OpenXQuestScores::match_odds_score(theoretical_spec_xlinks_beta, matched_spec_xlinks_beta, fragment_mass_tolerance_xlinks, fragment_mass_tolerance_unit_ppm, true, n_xlink_charges);
                match_odds = (match_odds_c_alpha + match_odds_x_alpha + match_odds_c_beta + match_odds_x_beta) / 4;
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
              RichPeakSpectrum theoretical_spec_common;
              RichPeakSpectrum theoretical_spec_xlinks;

              if (type_is_cross_link)
              {
                theoretical_spec_common.reserve(theoretical_spec_common_alpha.size() + theoretical_spec_common_beta.size());
                theoretical_spec_xlinks.reserve( theoretical_spec_xlinks_alpha.size() + theoretical_spec_xlinks_beta.size());
                theoretical_spec_common.insert(theoretical_spec_common.end(), theoretical_spec_common_alpha.begin(), theoretical_spec_common_alpha.end());
                theoretical_spec_common.insert(theoretical_spec_common.end(), theoretical_spec_common_beta.begin(), theoretical_spec_common_beta.end());
                theoretical_spec_xlinks.insert(theoretical_spec_xlinks.end(), theoretical_spec_xlinks_alpha.begin(), theoretical_spec_xlinks_alpha.end());
                theoretical_spec_xlinks.insert(theoretical_spec_xlinks.end(), theoretical_spec_xlinks_beta.begin(), theoretical_spec_xlinks_beta.end());
                theoretical_spec_common.sortByPosition();
                theoretical_spec_xlinks.sortByPosition();
              }
              else
              {
                theoretical_spec_common = theoretical_spec_common_alpha;
                theoretical_spec_xlinks = theoretical_spec_xlinks_alpha;
              }

               vector< double > xcorrx = OpenXQuestScores::xCorrelation(spectrum, theoretical_spec_xlinks, 5, 0.3);
               vector< double > xcorrc = OpenXQuestScores::xCorrelation(spectrum, theoretical_spec_common, 5, 0.2);

                // TODO save time: only needs to be done once per light spectrum, here it is repeated for cross-link each candidate
                vector< double > aucorrx = OpenXQuestScores::xCorrelation(spectrum, spectrum, 5, 0.3);
                vector< double > aucorrc = OpenXQuestScores::xCorrelation(spectrum, spectrum, 5, 0.2);
                double aucorr_sumx = accumulate(aucorrx.begin(), aucorrx.end(), 0.0);
                double aucorr_sumc = accumulate(aucorrc.begin(), aucorrc.end(), 0.0);
                double xcorrx_max = accumulate(xcorrx.begin(), xcorrx.end(), 0.0) / aucorr_sumx;
                double xcorrc_max = accumulate(xcorrc.begin(), xcorrc.end(), 0.0) / aucorr_sumc;
//                LOG_DEBUG << "xCorrelation X: " << xcorrx << endl;
//                LOG_DEBUG << "xCorrelation C: " << xcorrc << endl;
//                LOG_DEBUG << "Autocorr: " << aucorr << "\t Autocorr_sum: " << aucorr_sum << "\t xcorrx_max: " << xcorrx_max << "\t xcorrc_max: " << xcorrc_max << endl;

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
                double match_odds_weight = 1.973;
                double wTIC_weight = 12.829;
                double intsum_weight = 1.8;

//                double match_odds_weight = 1;  // 1.973;
//                double wTIC_weight = 24;           // 12.829;
//                double intsum_weight = 2.8;       // 1.8;

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
                for (Size k = 0; k < matched_spec_common_alpha.size(); ++k)
                {
                  PeptideHit::FragmentAnnotation frag_anno;
                  frag_anno.charge = static_cast<int>(theoretical_spec_common_alpha[matched_spec_common_alpha[k].first].getMetaValue("z"));
                  frag_anno.mz = spectrum[matched_spec_common_alpha[k].second].getMZ();
                  frag_anno.intensity = spectrum[matched_spec_common_alpha[k].second].getIntensity();
                  frag_anno.annotation = theoretical_spec_common_alpha[matched_spec_common_alpha[k].first].getMetaValue("IonName");
                  frag_annotations.push_back(frag_anno);
                }
                for (Size k = 0; k < matched_spec_common_beta.size(); ++k)
                {
                  PeptideHit::FragmentAnnotation frag_anno;
                  frag_anno.charge = static_cast<int>(theoretical_spec_common_beta[matched_spec_common_beta[k].first].getMetaValue("z"));
                  frag_anno.mz = spectrum[matched_spec_common_beta[k].second].getMZ();
                  frag_anno.intensity = spectrum[matched_spec_common_beta[k].second].getIntensity();
                  frag_anno.annotation = theoretical_spec_common_beta[matched_spec_common_beta[k].first].getMetaValue("IonName");
                  frag_annotations.push_back(frag_anno);
                }
                for (Size k = 0; k < matched_spec_xlinks_alpha.size(); ++k)
                {
                  PeptideHit::FragmentAnnotation frag_anno;
                  frag_anno.charge = static_cast<int>(theoretical_spec_xlinks_alpha[matched_spec_xlinks_alpha[k].first].getMetaValue("z"));
                  frag_anno.mz = spectrum[matched_spec_xlinks_alpha[k].second].getMZ();
                  frag_anno.intensity = spectrum[matched_spec_xlinks_alpha[k].second].getIntensity();
                  frag_anno.annotation = theoretical_spec_xlinks_alpha[matched_spec_xlinks_alpha[k].first].getMetaValue("IonName");
                  frag_annotations.push_back(frag_anno);
                }
                for (Size k = 0; k < matched_spec_xlinks_beta.size(); ++k)
                {
                  PeptideHit::FragmentAnnotation frag_anno;
                  frag_anno.charge = static_cast<int>(theoretical_spec_xlinks_beta[matched_spec_xlinks_beta[k].first].getMetaValue("z"));
                  frag_anno.mz = spectrum[matched_spec_xlinks_beta[k].second].getMZ();
                  frag_anno.intensity = spectrum[matched_spec_xlinks_beta[k].second].getIntensity();
                  frag_anno.annotation = theoretical_spec_xlinks_beta[matched_spec_xlinks_beta[k].first].getMetaValue("IonName");
                  frag_annotations.push_back(frag_anno);
                }
                LOG_DEBUG << "End writing annotations, size: " << frag_annotations.size() << endl;
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

              if (all_csms_spectrum[max_position].cross_link.getType() == TheoreticalSpectrumGeneratorXLinks::ProteinProteinCrossLink::CROSS)
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

              if (top_csms_spectrum[i].cross_link.getType() == TheoreticalSpectrumGeneratorXLinks::ProteinProteinCrossLink::MONO)
              {
                xltype = "mono-link";
              }
              else if (top_csms_spectrum[i].cross_link.getType() == TheoreticalSpectrumGeneratorXLinks::ProteinProteinCrossLink::LOOP)
              {
                xltype = "loop-link";
              }

              PeptideHit ph_alpha, ph_beta;
              // Set monolink as a modification or add MetaValue for cross-link identity and mass
              AASequence seq_alpha = top_csms_spectrum[i].cross_link.alpha;
              if (top_csms_spectrum[i].cross_link.getType() == TheoreticalSpectrumGeneratorXLinks::ProteinProteinCrossLink::MONO)
              {
                //AASequence seq_alpha = top_csms_spectrum[i].cross_link.alpha;
                vector< String > mods;
                const String residue = seq_alpha[alpha_pos].getOneLetterCode();
                ModificationsDB::getInstance()->getModificationsByDiffMonoMass(mods, residue, top_csms_spectrum[i].cross_link.cross_linker_mass, 0.01);
                LOG_DEBUG << "number of modifications fitting the mono-link: " << mods.size() << "\t" << mods << endl;
                if (mods.size() > 0)
                {
                  // TODO only valid if there are no alternatives, but we know we searched for a cross-linking reagent
                  LOG_DEBUG << "applied modification: " << mods[0] << endl;
                  seq_alpha.setModification(alpha_pos, mods[0]);
                }
                else
                {
//                  // TODO hardcoded for DSS, Xlink:DSS-NH2  not found by mass
//                  LOG_DEBUG << "applied modification: " << "Xlink:DSS-NH2" << endl;
//                  ph_alpha.setMetaValue("xl_mod",  "Xlink:DSS-NH2");
                }
              }
              else
              {
                // TODO hardcoded for DSS, make this an input parameter or something, NO UNIMOD ACCESSION AVAILBALE, for now name and mass
              ph_alpha.setMetaValue("xl_mod", top_csms_spectrum[i].cross_link.cross_linker_name);
              ph_alpha.setMetaValue("xl_mass", DataValue(top_csms_spectrum[i].cross_link.cross_linker_mass));
              }


              if (top_csms_spectrum[i].cross_link.getType() == TheoreticalSpectrumGeneratorXLinks::ProteinProteinCrossLink::LOOP)
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

              ph_alpha.setFragmentAnnotations(top_csms_spectrum[i].frag_annotations);
              LOG_DEBUG << "Annotations of size " << ph_alpha.getFragmentAnnotations().size() << endl;
              phs.push_back(ph_alpha);

              if (top_csms_spectrum[i].cross_link.getType() == TheoreticalSpectrumGeneratorXLinks::ProteinProteinCrossLink::CROSS)
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

                phs.push_back(ph_beta);
              }

              peptide_id.setRT(spectrum.getRT());
              peptide_id.setMZ(precursor_mz);
              String specIDs = spectra[scan_index].getNativeID();
              peptide_id.setMetaValue("spectrum_reference", specIDs);
//              peptide_id.setMetaValue("spec_heavy_RT", spectra[scan_index_heavy].getRT());
//              peptide_id.setMetaValue("spec_heavy_MZ", spectra[scan_index_heavy].getPrecursors()[0].getMZ());
//              peptide_id.setMetaValue("spectrum_reference", spectra[scan_index].getNativeID());
//              peptide_id.setMetaValue("spectrum_reference_heavy", spectra[scan_index_heavy].getNativeID());
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

      writeXQuestXML(out_xquest, base_name, peptide_ids, all_top_csms, spectra);
      writeXQuestXMLSpec(base_name, all_top_csms, spectra);
    }
    progresslogger.endProgress();

    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{

  TOPPxQuest tool;

  return tool.main(argc, argv);
}

