// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Authors: Timo Sachsenberg, Samuel Wein, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h> // for "median"

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/METADATA/SpectrumSettings.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>

// file types
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/SVOutStream.h>

// digestion enzymes
#include <OpenMS/CHEMISTRY/RNaseDigestion.h>
#include <OpenMS/CHEMISTRY/RNaseDB.h>

// ribonucleotides
#include <OpenMS/CHEMISTRY/RibonucleotideDB.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/CHEMISTRY/ModifiedNASequenceGenerator.h>
#include <OpenMS/CHEMISTRY/NASequence.h>

// preprocessing and filtering of spectra
#include <OpenMS/FILTERING/TRANSFORMERS/ThresholdMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/NLargest.h>
#include <OpenMS/FILTERING/TRANSFORMERS/WindowMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>

// spectra comparison
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignment.h>
#include <OpenMS/ANALYSIS/ID/MetaboliteSpectralMatching.h>

// post-processing of results
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/FILTERING/ID/IDFilter.h>

#include <boost/regex.hpp>

#include <QtCore/QProcess>

#include <iostream>
#include <algorithm>
#include <vector>
#include <map>

// multithreading
#ifdef _OPENMP
#include <omp.h>
#define NUMBER_OF_THREADS (omp_get_num_threads())
#else
#define NUMBER_OF_THREADS (1)
#endif

//#define DEBUG_NASEARCH

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_NucleicAcidSearchEngine NucleicAcidSearchEngine

    @brief Searches an mzML file for specified nucleic acid sequences.

    Given a FASTA file containing oligonucleotide sequences (and optionally decoys) and a mzML file from a nucleic acid mass spec experiment:
    - Generate a list of digest fragments from the FASTA file with a specified RNAse
    - Search the mzML input for MS2 spectra with parent masses corresponding to any of these digests
    - Match the MS2 spectra to theoretically generated spectra
    - Score the resulting matches

    Output is in the form of an mzTab file containing the search results and an optional idXML file containing identifications.

    Modified nucleic acids can either be included in the FASTA input file, or set as @p variable or @p fixed modifications in the tool options.
    All modification syntax is taken from the Modomics database (http://modomics.genesilico.pl/)


    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_NucleicAcidSearchEngine.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_NucleicAcidSearchEngine.html
*/

class NucleicAcidSearchEngine :
  public TOPPBase
{
  using ConstRibonucleotidePtr = const Ribonucleotide*;

public:
  NucleicAcidSearchEngine() :
    TOPPBase("NucleicAcidSearchEngine", "Annotate nucleic acid identifications to MS/MS spectra.", false),
    fragment_ion_codes_(ListUtils::create<String>("a-B,a,b,c,d,w,x,y,z"))
  {
  }

protected:
  vector<String> fragment_ion_codes_;

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "Input file: spectra");
    setValidFormats_("in", ListUtils::create<String>("mzML"));

    registerInputFile_("database", "<file>", "", "Input file: sequence database");
    setValidFormats_("database", ListUtils::create<String>("fasta"));

    registerOutputFile_("out", "<file>", "", "Output file: mzTab");
    setValidFormats_("out", ListUtils::create<String>("mzTab"));

    registerOutputFile_("id_out", "<file>", "", "Output file: idXML (for visualization in TOPPView)", false);
    setValidFormats_("id_out", ListUtils::create<String>("idXML"));

    registerOutputFile_("lfq_out", "<file>", "", "Output file: Targets for label-free quantification using FeatureFinderMetaboIdent ('id' input)", false);
    setValidFormats_("lfq_out", vector<String>(1, "tsv"));

    registerOutputFile_("theo_ms2_out", "<file>", "", "Output file: theoretical MS2 spectra for precursor mass matches", false, true);
    setValidFormats_("theo_ms2_out", ListUtils::create<String>("mzML"));
    registerOutputFile_("exp_ms2_out", "<file>", "", "Output file: experimental MS2 spectra for precursor mass matches", false, true);
    setValidFormats_("exp_ms2_out", ListUtils::create<String>("mzML"));

    registerFlag_("decharge_ms2", "Decharge the MS2 spectra for scoring", true);

    registerTOPPSubsection_("precursor", "Precursor (Parent Ion) Options");
    registerDoubleOption_("precursor:mass_tolerance", "<tolerance>", 10.0, "Precursor mass tolerance (+/- around uncharged precursor mass)", false);

    registerStringOption_("precursor:mass_tolerance_unit", "<unit>", "ppm", "Unit of precursor mass tolerance", false, false);
    setValidStrings_("precursor:mass_tolerance_unit", ListUtils::create<String>("Da,ppm"));

    registerIntOption_("precursor:min_charge", "<num>", -1, "Minimum precursor charge to be considered", false, false);
    registerIntOption_("precursor:max_charge", "<num>", -20, "Maximum precursor charge to be considered", false, false);

    registerFlag_("precursor:include_unknown_charge", "Include MS2 spectra with unknown precursor charge - try to match them in any possible charge between 'min_charge' and 'max_charge', at the risk of a higher error rate", false);

    registerFlag_("precursor:use_avg_mass", "Use average instead of monoisotopic precursor masses (appropriate for low-resolution instruments)", false);

    // Whether to look for precursors with salt adducts
    registerFlag_("precursor:use_adducts", "Consider possible salt adducts (see 'precursor:potential_adducts') when matching precursor masses", false);
    registerStringList_("precursor:potential_adducts", "<list>", ListUtils::create<String>("Na:+"), "Adducts considered to explain mass differences. Format: 'Element:Charge(+/-)', i.e. the number of '+' or '-' indicates the charge, e.g. 'Ca:++' indicates +2. Only used if 'precursor:use_adducts' is set.", false, false);

    IntList isotopes = {0, 1, 2, 3, 4};
    registerIntList_("precursor:isotopes", "<list>", isotopes, "Correct for mono-isotopic peak misassignments. E.g.: 1 = precursor may be misassigned to the first isotopic peak. Ignored if 'use_avg_mass' is set.", false, false);

    registerTOPPSubsection_("fragment", "Fragment (Product Ion) Options");
    registerDoubleOption_("fragment:mass_tolerance", "<tolerance>", 10.0, "Fragment mass tolerance (+/- around fragment m/z)", false);

    registerStringOption_("fragment:mass_tolerance_unit", "<unit>", "ppm", "Unit of fragment mass tolerance", false, false);
    setValidStrings_("fragment:mass_tolerance_unit", ListUtils::create<String>("Da,ppm"));

    registerStringList_("fragment:ions", "<choice>", fragment_ion_codes_, "Fragment ions to include in theoretical spectra", false);
    setValidStrings_("fragment:ions", fragment_ion_codes_);

    registerTOPPSubsection_("modifications", "Modifications Options");

    // add modified ribos from database
    vector<String> all_mods;
    for (auto r : *RibonucleotideDB::getInstance())
    {
      if (r->isModified())
      {
        String code = r->getCode();
        // commas aren't allowed in parameter string restrictions:
        all_mods.push_back(code.substitute(',', '_'));
      }
    }

    registerStringList_("modifications:fixed", "<mods>", ListUtils::create<String>(""), "Fixed modifications", false);
    setValidStrings_("modifications:fixed", all_mods);
    registerStringList_("modifications:variable", "<mods>", ListUtils::create<String>(""), "Variable modifications", false);
    setValidStrings_("modifications:variable", all_mods);
    registerIntOption_("modifications:variable_max_per_oligo", "<num>", 2, "Maximum number of residues carrying a variable modification per candidate oligonucleotide", false, false);

    registerTOPPSubsection_("oligo", "Oligonucleotide Options");
    registerIntOption_("oligo:min_size", "<num>", 5, "Minimum size an oligonucleotide must have after digestion to be considered in the search", false);
    registerIntOption_("oligo:missed_cleavages", "<num>", 1, "Number of missed cleavages", false, false);

    StringList all_enzymes;
    RNaseDB::getInstance()->getAllNames(all_enzymes);
    registerStringOption_("oligo:enzyme", "<choice>", "no cleavage", "The enzyme used for RNA digestion", false);
    setValidStrings_("oligo:enzyme", all_enzymes);

    registerTOPPSubsection_("report", "Reporting Options");
    registerIntOption_("report:top_hits", "<num>", 1, "Maximum number of top-scoring hits per spectrum that are reported ('0' for all hits)", false, true);
    setMinInt_("report:top_hits", 0);

    registerTOPPSubsection_("fdr", "False Discovery Rate Options");
    registerStringOption_("fdr:decoy_pattern", "<string>", "", "String used as part of the accession to annotate decoy sequences (e.g. 'DECOY_'). Leave empty to skip the FDR/q-value calculation.", false);
    registerDoubleOption_("fdr:cutoff", "<value>", 1.0, "Cut-off for FDR filtering; search hits with higher q-values will be removed", false);
    setMinFloat_("fdr:cutoff", 0.0);
    setMaxFloat_("fdr:cutoff", 1.0);
    registerFlag_("fdr:remove_decoys", "Do not score hits to decoy sequences and remove them when filtering");
  }


  struct PrecursorInfo
  {
    Size scan_index;
    Int charge;
    Size isotope;
    String adduct;

    PrecursorInfo(Size scan_index, Int charge, Size isotope, const String& adduct):
      scan_index(scan_index), charge(charge), isotope(isotope), adduct(adduct)
    {
    }
  };


  // slimmer structure to store basic hit information
  struct AnnotatedHit
  {
    IdentificationData::IdentifiedOligoRef oligo_ref;
    SignedSize mod_index; // enumeration index of the modification
    double precursor_error_ppm; // precursor mass error in ppm
    vector<PeptideHit::PeakAnnotation> annotations; // peak/ion annotations
    Int charge;
    String adduct;
  };

  typedef multimap<double, AnnotatedHit, greater<double>> HitsByScore;

  // query modified residues from database
  vector<ConstRibonucleotidePtr> getModifications_(const set<String>& mod_names)
  {
    vector<ConstRibonucleotidePtr> modifications;
    for (String m : mod_names)
    {
      m.substitute('_', ',');
      ConstRibonucleotidePtr rm = RibonucleotideDB::getInstance()->getRibonucleotide(m);
      modifications.push_back(rm);
    }
    return modifications;
  }

  // check for minimum size
  class HasInvalidLength
  {
    Size min_size_;
  public:
    explicit HasInvalidLength(Size min_size)
      : min_size_(min_size)
    {
    }
    bool operator()(const NASequence& s) { return s.size() < min_size_; }
  };

  // turn an adduct string (param. "precursor:potential_adducts") into a formula
  EmpiricalFormula parseAdduct_(const String& adduct)
  {
    StringList parts;
    adduct.split(':', parts);
    if (parts.size() != 2)
    {
      String error = "entry in parameter 'precursor:potential_adducts' does not have two parts separated by ':'";
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    error, adduct);
    }

    // determine charge of adduct (by number of '+' or '-')
    Int pos_charge = parts[1].size() - parts[1].remove('+').size();
    Int neg_charge = parts[1].size() - parts[1].remove('-').size();
    LOG_DEBUG << ": " << pos_charge - neg_charge << endl;
    if (pos_charge > 0 && neg_charge > 0)
    {
      String error = "entry in parameter 'precursor:potential_adducts' mixes positive and negative charges";
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    error, adduct);
    }

    // generate the formula for the adduct at neutral charge (!) -
    // compensate for intrinsic charge by adding/removing hydrogens:
    EmpiricalFormula ef(parts[0]);
    if (pos_charge > 0)
    {
      ef -= EmpiricalFormula("H" + String(pos_charge));
    }
    else if ((neg_charge > 0) && (parts[0] != "H-1"))
    {
      ef += EmpiricalFormula("H" + String(neg_charge));
    }
    return ef;
  }

  // spectrum must not contain 0 intensity peaks and must be sorted by m/z
  void deisotopeAndSingleChargeMSSpectrum_(
    MSSpectrum& in,
    Int min_charge,
    Int max_charge,
    double fragment_tolerance,
    bool fragment_unit_ppm,
    bool keep_only_deisotoped = false,
    Size min_isopeaks = 3,
    Size max_isopeaks = 10,
    bool make_single_charged = true)
  {
    if (in.empty()) return;

    MSSpectrum old_spectrum = in;

    // determine charge seeds and extend them
    vector<Size> mono_isotopic_peak(old_spectrum.size(), 0);
    vector<Int> features(old_spectrum.size(), -1);
    Int feature_number = 0;

    bool negative_mode = (max_charge < 0);
    Int step = negative_mode ? -1 : 1;

    for (Size current_peak = 0; current_peak != old_spectrum.size(); ++current_peak)
    {
      double current_mz = old_spectrum[current_peak].getPosition()[0];

      for (Int q = max_charge; abs(q) >= abs(min_charge); q -= step) // important: test charge hypothesis from high to low (in terms of absolute values)
      {
        // try to extend isotopes from mono-isotopic peak
        // if extension larger then min_isopeaks possible:
        //   - save charge q in mono_isotopic_peak[]
        //   - annotate all isotopic peaks with feature number
        if (features[current_peak] == -1) // only process peaks which have no assigned feature number
        {
          bool has_min_isopeaks = true;
          vector<Size> extensions;
          for (Size i = 0; i < max_isopeaks; ++i)
          {
            double expected_mz = current_mz + i * Constants::C13C12_MASSDIFF_U / abs(q);
            Size p = old_spectrum.findNearest(expected_mz);
            double tolerance_dalton = fragment_unit_ppm ? fragment_tolerance * old_spectrum[p].getPosition()[0] * 1e-6 : fragment_tolerance;
            if (fabs(old_spectrum[p].getPosition()[0] - expected_mz) > tolerance_dalton) // test for missing peak
            {
              if (i < min_isopeaks)
              {
                has_min_isopeaks = false;
              }
              break;
            }
            else
            {
/*
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
*/
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
        else // make singly charged
        {
          Peak1D p = old_spectrum[i];
          if (negative_mode) // z < 0 in this case
          {
            z = abs(z);
            p.setMZ(p.getMZ() * z + (z - 1) * Constants::PROTON_MASS_U);
          }
          else
          {
            p.setMZ(p.getMZ() * z - (z - 1) * Constants::PROTON_MASS_U);
          }
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
          else // make singly charged
          {
            Peak1D p = old_spectrum[i];
            if (negative_mode) // z < 0 in this case
            {
              z = abs(z);
              p.setMZ(p.getMZ() * z + (z - 1) * Constants::PROTON_MASS_U);
            }
            else
            {
              p.setMZ(p.getMZ() * z - (z - 1) * Constants::PROTON_MASS_U);
            }
            in.push_back(p);
          }
        }
      }
    }

    in.sortByPosition();
  }


  void preprocessSpectra_(PeakMap& exp, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, bool single_charge_spectra, bool negative_mode, Int min_charge, Int max_charge, bool include_unknown_charge)
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

    // Note: we expect a higher number for NA than e.g., for peptides
    filter_param.setValue("peakcount", 50, "The number of peaks that should be kept.");
    filter_param.setValue("movetype", "jump", "Whether sliding window (one peak steps) or jumping window (window size steps) should be used.");
    window_mower_filter.setParameters(filter_param);

    // Note: we expect a higher number for NA than e.g., for peptides
    NLargest nlargest_filter = NLargest(1000);

    Size n_zero_charge = 0, n_inferred_charge = 0;

#pragma omp parallel for
    for (SignedSize exp_index = 0; exp_index < (SignedSize)exp.size();
         ++exp_index)
    {
      MSSpectrum& spec = exp[exp_index];

      // sort by mz
      spec.sortByPosition();

      if (spec.getPrecursors().empty()) continue; // this shouldn't happen
      Int precursor_charge = spec.getPrecursors()[0].getCharge();
      if (precursor_charge == 0) // no charge information
      {
        n_zero_charge++;
        // maybe we are able to infer the charge state:
        if (spec.getPrecursors().size() > 1) // multiplexed PRM experiment
        {
          // all precursors belong to the same parent, but with different charge
          // states; we want to find the precursor with highest charge and infer
          // its charge state:
          map<double, Size> precursors; // precursor: m/z -> index
          for (Size i = 0; i < spec.getPrecursors().size(); ++i)
          {
            precursors[spec.getPrecursors()[i].getMZ()] = i;
          }
          double mz1 = precursors.begin()->first;
          double mz2 = (++precursors.begin())->first;
          double mz_ratio = mz1 / mz2;

          Int step = negative_mode ? -1 : 1;
          Int inferred_charge = 0;
          for (Int charge = max_charge; abs(charge) > abs(min_charge);
               charge -= step)
          {
            double charge_ratio = (abs(charge) - 1.0) / abs(charge);
            double ratios_ratio = mz_ratio / charge_ratio;
            if ((ratios_ratio > 0.99) && (ratios_ratio < 1.01))
            {
              inferred_charge = charge;
              break;
            }
          }
          if (inferred_charge == 0)
          {
            LOG_ERROR << "Error: unable to determine charge state for spectrum '" << spec.getNativeID() << "' based on precursor m/z values " << mz1 << " and " << mz2 << endl;
          }
          else
          {
            ++n_inferred_charge;
            LOG_DEBUG << "Inferred charge state " << inferred_charge << " for spectrum '" << spec.getNativeID() << "'" << endl;
            // keep only precursor with highest charge, set inferred charge:
            Precursor prec = spec.getPrecursors()[precursors.begin()->second];
            prec.setCharge(abs(inferred_charge));
            spec.setPrecursors(vector<Precursor>(1, prec));
          }
        }
      }

      // deisotope
      Int coef = negative_mode ? -1 : 1;
      // @TODO: what happens here if "precursor_charge" is zero?
      deisotopeAndSingleChargeMSSpectrum_(spec, coef, coef * precursor_charge, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, false, 3, 20, single_charge_spectra);

      // remove noise
      window_mower_filter.filterPeakSpectrum(spec);
      nlargest_filter.filterPeakSpectrum(spec);

      // sort (nlargest changes order)
      spec.sortByPosition();
    }

    if (n_zero_charge)
    {
      LOG_WARN << "Warning: no charge state information available for " << n_zero_charge << " out of " << exp.size() << " spectra." << endl;
      if (n_inferred_charge)
      {
        LOG_INFO << "Inferred charge states for " << n_inferred_charge << " spectra." << endl;
      }
      if (n_zero_charge - n_inferred_charge > 0)
      {
        LOG_INFO << "Spectra without charge information will be " << (include_unknown_charge ? "included in the processing" : "skipped") << " (see parameter 'precursor:include_unknown_charge')" << endl;
      }
    }
  }


  void insertSpectrumPrecursorMass_(multimap<double, PrecursorInfo>& precursor_mass_map, Size scan_index, double mz, Int charge, Size isotope, double adduct_mass, const String& adduct, bool negative_mode)
  {
    // we want to calculate the unadducted (!) precursor mass at neutral charge:
    double mass = mz * charge - adduct_mass;
    // compensate for loss or gain of protons that confer the charge:
    if (negative_mode)
    {
      mass += Constants::PROTON_MASS_U * charge;
    }
    else
    {
      mass -= Constants::PROTON_MASS_U * charge;
    }
    // correct for precursor not being the monoisotopic peak:
    if (isotope > 0)
    {
      mass -= isotope * Constants::C13C12_MASSDIFF_U;
    }

    PrecursorInfo info(scan_index, charge, isotope, adduct);
    precursor_mass_map.insert(make_pair(mass, info));
  }


  void postProcessHits_(const PeakMap& exp,
                        vector<HitsByScore>& annotated_hits,
                        IdentificationData& id_data,
                        const vector<ConstRibonucleotidePtr>& fixed_mods,
                        const vector<ConstRibonucleotidePtr>& variable_mods,
                        Size max_variable_mods_per_oligo,
                        bool negative_mode)
  {
    IdentificationData::InputFileRef file_ref = id_data.getInputFiles().begin();
    IdentificationData::ScoreTypeRef score_ref =
      id_data.getScoreTypes().begin();

// @TODO: change OpenMP schedule from default ("static") to "dynamic"/"guided"?
#pragma omp parallel for
    for (SignedSize scan_index = 0;
         scan_index < (SignedSize)annotated_hits.size(); ++scan_index)
    {
      if (annotated_hits[scan_index].empty()) continue;

      const MSSpectrum& spectrum = exp[scan_index];
      IdentificationData::DataQuery query(spectrum.getNativeID(), file_ref,
                                          spectrum.getRT(),
                                          spectrum.getPrecursors()[0].getMZ());
      query.setMetaValue("scan_index", static_cast<unsigned int>(scan_index));
      IdentificationData::DataQueryRef query_ref;
#pragma omp critical (id_data_access)
      query_ref = id_data.registerDataQuery(query);

      // create full oligo hit structure from annotated hits
      for (const auto& pair : annotated_hits[scan_index])
      {
        const AnnotatedHit& hit = pair.second;
        double score = pair.first;
        // get unmodified string
        NASequence seq = hit.oligo_ref->sequence;
        LOG_DEBUG << "Hit sequence: " << seq << endl;
        // reapply modifications (because for memory reasons we only stored the index and recreation is fast)
        vector<NASequence> all_modified_oligos;
        ModifiedNASequenceGenerator::applyFixedModifications(
          fixed_mods.begin(), fixed_mods.end(), seq);
        ModifiedNASequenceGenerator::applyVariableModifications(
          variable_mods.begin(), variable_mods.end(), seq,
          max_variable_mods_per_oligo, all_modified_oligos);

        // transfer parent matches from unmodified oligo:
        IdentificationData::IdentifiedOligo oligo = *hit.oligo_ref;
        oligo.sequence = all_modified_oligos[hit.mod_index];
        IdentificationData::IdentifiedOligoRef oligo_ref;
#pragma omp critical (id_data_access)
        oligo_ref = id_data.registerIdentifiedOligo(oligo);

        Int charge = hit.charge;
        if ((charge > 0) && negative_mode) charge = -charge;
        IdentificationData::MoleculeQueryMatch match(oligo_ref, query_ref,
                                                     charge);
        match.addScore(score_ref, score, id_data.getCurrentProcessingStep());
        match.peak_annotations[id_data.getCurrentProcessingStep()] =
          hit.annotations;
        // @TODO: add a field for this to "IdentificationData::MoleculeQueryMatch"?
        match.setMetaValue(Constants::PRECURSOR_ERROR_PPM_USERPARAM,
                           hit.precursor_error_ppm);
        if (!hit.adduct.empty()) match.setMetaValue("adduct", hit.adduct);
#pragma omp critical (id_data_access)
        id_data.registerMoleculeQueryMatch(match);
      }
    }
    id_data.cleanup();
  }


  void registerIDMetaData_(
    IdentificationData& id_data, const String& in_mzml,
    const vector<String>& primary_files,
    const IdentificationData::DBSearchParam& search_param)
  {
    IdentificationData::InputFileRef file_ref =
      id_data.registerInputFile(in_mzml);
    IdentificationData::ScoreType score("hyperscore", true);
    IdentificationData::ScoreTypeRef score_ref =
      id_data.registerScoreType(score);
    IdentificationData::DataProcessingSoftware software(toolName_(), version_);
    // if we are in test mode just overwrite with a generic version
    if (test_mode_)
    {
      software.setVersion("test");
    }
    software.assigned_scores.push_back(score_ref);

    IdentificationData::ProcessingSoftwareRef software_ref =
      id_data.registerDataProcessingSoftware(software);
    IdentificationData::SearchParamRef search_ref =
      id_data.registerDBSearchParam(search_param);
    // @TODO: add suitable data processing action
    IdentificationData::DataProcessingStep step(
      software_ref, vector<IdentificationData::InputFileRef>(1, file_ref),
      primary_files);
    IdentificationData::ProcessingStepRef step_ref =
      id_data.registerDataProcessingStep(step, search_ref);
    // reference this step in all following ID data items, if applicable:
    id_data.setCurrentProcessingStep(step_ref);
  }


  void calculateAndFilterFDR_(IdentificationData& id_data, bool only_top_hits)
  {
    IdentificationData::ScoreTypeRef score_ref =
      id_data.findScoreType("hyperscore").first;
    FalseDiscoveryRate fdr;
    Param fdr_params = fdr.getDefaults();
    fdr_params.setValue("use_all_hits", only_top_hits ? "true" : "false");
    bool remove_decoys = getFlag_("fdr:remove_decoys");
    fdr_params.setValue("add_decoy_peptides", remove_decoys ? "false" : "true");
    fdr.setParameters(fdr_params);
    IdentificationData::ScoreTypeRef fdr_ref =
      fdr.applyToQueryMatches(id_data, score_ref);
    double fdr_cutoff = getDoubleOption_("fdr:cutoff");
    if (remove_decoys) // remove references to decoys from shared oligos
    {
      IDFilter::removeDecoys(id_data);
    }
    if (fdr_cutoff < 1.0)
    {
      IDFilter::filterQueryMatchesByScore(id_data, fdr_ref, fdr_cutoff);
      LOG_INFO << "Search hits after FDR filtering: "
               << id_data.getMoleculeQueryMatches().size() << endl;
    }
  }


  void generateLFQInput_(IdentificationData& id_data, const String& out_file)
  {
    // mapping: (oligo, adduct) -> charge -> RTs
    map<pair<NASequence, String>, map<Int, vector<double>>> rt_info;
    for (const IdentificationData::MoleculeQueryMatch& match :
           id_data.getMoleculeQueryMatches())
    {
      auto key = make_pair(match.getIdentifiedOligoRef()->sequence,
                           match.getMetaValue("adduct"));
      rt_info[key][match.charge].push_back(match.data_query_ref->rt);
    }

    SVOutStream tsv(out_file);
    tsv.modifyStrings(false);
    tsv << "CompoundName" << "SumFormula" << "Mass" << "Charge"
        << "RetentionTime" << "RetentionTimeRange" << "IsoDistribution" << endl;
    for (const auto& entry : rt_info)
    {
      String name = entry.first.first.toString();
      EmpiricalFormula ef = entry.first.first.getFormula();
      if (!entry.first.second.empty()) // adduct given
      {
        name += "+[" + entry.first.second + "]";
        ef += parseAdduct_(entry.first.second);
      }
      // @TODO: use charge-specific RTs?
      vector<Int> charges;
      vector<double> rts;
      for (const auto& pair : entry.second)
      {
        charges.push_back(pair.first);
        rts.insert(rts.end(), pair.second.begin(), pair.second.end());
      }
      tsv << name << ef << 0 << ListUtils::concatenate(charges, ",");
      // @TODO: do something more sophisticated for the RT?
      tsv << Math::median(rts.begin(), rts.end(), false) << 0 << 0 << endl;
    }
  }


  ExitCodes main_(int, const char**)
  {
    ProgressLogger progresslogger;
    progresslogger.setLogType(log_type_);
    String in_mzml = getStringOption_("in");
    String in_db = getStringOption_("database");
    String out = getStringOption_("out");
    String id_out = getStringOption_("id_out");
    String lfq_out = getStringOption_("lfq_out");
    String theo_ms2_out = getStringOption_("theo_ms2_out");
    String exp_ms2_out = getStringOption_("exp_ms2_out");
    bool use_avg_mass = getFlag_("precursor:use_avg_mass");

    IdentificationData::DBSearchParam search_param;
    search_param.molecule_type = IdentificationData::MoleculeType::RNA;
    search_param.mass_type = (use_avg_mass ?
                              IdentificationData::MassType::AVERAGE :
                              IdentificationData::MassType::MONOISOTOPIC);
    search_param.database = in_db;
    Int min_charge = getIntOption_("precursor:min_charge");
    Int max_charge = getIntOption_("precursor:max_charge");
    // @TODO: allow zero to mean "any charge state in the data"?
    if ((min_charge == 0) || (max_charge == 0))
    {
      LOG_ERROR << "Error: invalid charge state 0" << endl;
      return ILLEGAL_PARAMETERS;
    }
    // charges can be positive or negative, depending on data acquisition mode:
    if (((min_charge < 0) && (max_charge > 0)) ||
        ((min_charge > 0) && (max_charge < 0)))
    {
      LOG_ERROR << "Error: mixing positive and negative charges is not allowed"
                << endl;
      return ILLEGAL_PARAMETERS;
    }
    // min./max. are based on absolute value:
    if (abs(max_charge) < abs(min_charge)) swap(min_charge, max_charge);
    bool negative_mode = (max_charge < 0);
    Int step = negative_mode ? -1 : 1;
    for (Int charge = min_charge; abs(charge) <= abs(max_charge);
         charge += step)
    {
      search_param.charges.insert(charge);
    }
    search_param.precursor_mass_tolerance =
      getDoubleOption_("precursor:mass_tolerance");
    search_param.precursor_tolerance_ppm =
      (getStringOption_("precursor:mass_tolerance_unit") == "ppm");
    search_param.fragment_mass_tolerance =
      getDoubleOption_("fragment:mass_tolerance");
    search_param.fragment_tolerance_ppm =
      (getStringOption_("fragment:mass_tolerance_unit") == "ppm");
    search_param.min_length = getIntOption_("oligo:min_size");

    StringList fixed_mod_names = getStringList_("modifications:fixed");
    search_param.fixed_mods.insert(fixed_mod_names.begin(),
                                   fixed_mod_names.end());

    StringList var_mod_names = getStringList_("modifications:variable");
    search_param.variable_mods.insert(var_mod_names.begin(),
                                      var_mod_names.end());

    vector<ConstRibonucleotidePtr> fixed_modifications =
      getModifications_(search_param.fixed_mods);
    vector<ConstRibonucleotidePtr> variable_modifications =
      getModifications_(search_param.variable_mods);

    // @TODO: add slots for these to "IdentificationData::DBSearchParam"?
    IntList precursor_isotopes = (use_avg_mass ? vector<Int>(1, 0) :
                                  getIntList_("precursor:isotopes"));
    Size max_variable_mods_per_oligo =
      getIntOption_("modifications:variable_max_per_oligo");
    Size report_top_hits = getIntOption_("report:top_hits");

    StringList potential_adducts =
      getStringList_("precursor:potential_adducts");
    map<double, String> adduct_masses;
    adduct_masses[0.0] = ""; // always consider "no adduct"
    bool use_adducts = getFlag_("precursor:use_adducts");
    bool include_unknown_charge = getFlag_("precursor:include_unknown_charge");
    bool single_charge_spectra = getFlag_("decharge_ms2");
    if (use_adducts)
    {
      for (const String& adduct : potential_adducts)
      {
        EmpiricalFormula ef = parseAdduct_(adduct);
        double mass = use_avg_mass ? ef.getAverageWeight() : ef.getMonoWeight();
        adduct_masses[mass] = adduct;
        LOG_DEBUG << "Added adduct: " << adduct << ", mass: " << mass << endl;
      }
    }

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

    progresslogger.startProgress(0, 1, "filtering spectra...");
    // @TODO: move this into the loop below (run only when checks pass)
    preprocessSpectra_(spectra, search_param.fragment_mass_tolerance,
                       search_param.fragment_tolerance_ppm,
                       single_charge_spectra, negative_mode, min_charge,
                       max_charge, include_unknown_charge);
    progresslogger.endProgress();
    LOG_DEBUG << "preprocessed spectra: " << spectra.getNrSpectra() << endl;

    // build multimap of precursor mass to scan index (and other information):
    multimap<double, PrecursorInfo> precursor_mass_map;
    for (PeakMap::ConstIterator s_it = spectra.begin(); s_it != spectra.end();
         ++s_it)
    {
      int scan_index = s_it - spectra.begin();
      const vector<Precursor>& precursors = s_it->getPrecursors();

      // there should be only one precursor and MS2 should contain at least a
      // few peaks to be considered (at least one per nucleotide in the chain):
      if ((precursors.size() != 1) || (s_it->size() < search_param.min_length))
      {
        continue;
      }

      set<Int> possible_charges;
      Int precursor_charge = precursors[0].getCharge();
      if (precursor_charge == 0) // charge information missing
      {
        if (include_unknown_charge)
        {
          possible_charges = search_param.charges; // try all allowed charges
        }
        else
        {
          continue; // skip
        }
      }
      // compare to charge parameters (the charge value in mzML seems to be
      // always positive, so compare by absolute value in negative mode):
      else if ((negative_mode &&
                ((precursor_charge > abs(*search_param.charges.begin())) ||
                 (precursor_charge < abs(*(--search_param.charges.end()))))) ||
               (!negative_mode &&
                ((precursor_charge < *search_param.charges.begin()) ||
                 (precursor_charge > *(--search_param.charges.end())))))
      {
        continue; // charge not in user-supplied range
      }
      else
      {
        possible_charges.insert(precursor_charge); // only one possibility
      }

      double precursor_mz = precursors[0].getMZ();

      for (Int precursor_charge : possible_charges)
      {
        precursor_charge = abs(precursor_charge); // adjust for neg. mode

        // calculate precursor mass (optionally corrected for adducts and peak
        // misassignment) and map it to MS scan index:
        for (const auto& adduct_pair : adduct_masses)
        {
          for (Int isotope_number : precursor_isotopes)
          {
            insertSpectrumPrecursorMass_(precursor_mass_map, scan_index,
                                         precursor_mz, precursor_charge,
                                         isotope_number, adduct_pair.first,
                                         adduct_pair.second, negative_mode);
          }
        }
      }
    }

    // create spectrum generator
    TheoreticalSpectrumGenerator spectrum_generator;
    Param param = spectrum_generator.getParameters();
    vector<String> temp = getStringList_("fragment:ions");
    set<String> selected_ions(temp.begin(), temp.end());
    for (const auto& code : fragment_ion_codes_)
    {
      String param_name = "add_" + code + "_ions";
      if (selected_ions.count(code))
      {
        param.setValue(param_name, "true");
      }
      else
      {
        param.setValue(param_name, "false");
      }
    }
    param.setValue("add_first_prefix_ion", "true");
    param.setValue("add_metainfo", "true");
    spectrum_generator.setParameters(param);

    vector<HitsByScore> annotated_hits(spectra.size());
    MSExperiment exp_ms2_spectra, theo_ms2_spectra; // debug output

    progresslogger.startProgress(0, 1, "loading database from FASTA file...");
    vector<FASTAFile::FASTAEntry> fasta_db;
    FASTAFile().load(in_db, fasta_db);
    progresslogger.endProgress();

    search_param.missed_cleavages = getIntOption_("oligo:missed_cleavages");
    String enzyme_name = getStringOption_("oligo:enzyme");
    search_param.digestion_enzyme =
      RNaseDB::getInstance()->getEnzyme(enzyme_name);
    RNaseDigestion digestor;
    digestor.setEnzyme(search_param.digestion_enzyme);
    digestor.setMissedCleavages(search_param.missed_cleavages);
    // set minimum size of oligo after digestion
    Size min_oligo_length = getIntOption_("oligo:min_size");

    IdentificationData id_data;
    vector<String> primary_files;
    spectra.getPrimaryMSRunPath(primary_files);
    // this also sets the current processing step in "id_data":
    registerIDMetaData_(id_data, in_mzml, primary_files, search_param);
    String decoy_pattern = getStringOption_("fdr:decoy_pattern");

    LOG_INFO << "Performing in-silico digestion..." << endl;
    IdentificationDataConverter::importSequences(
      id_data, fasta_db, IdentificationData::MoleculeType::RNA, decoy_pattern);
    digestor.digest(id_data, min_oligo_length);

    String msg =  "scoring oligonucleotide models against spectra...";
    progresslogger.startProgress(0, id_data.getIdentifiedOligos().size(), msg);
    Size hit_counter = 0;

    // keep a list of (references to) oligos in the original digest:
    vector<IdentificationData::IdentifiedOligoRef> digest;
    digest.reserve(id_data.getIdentifiedOligos().size());
    for (IdentificationData::IdentifiedOligoRef it =
           id_data.getIdentifiedOligos().begin(); it !=
           id_data.getIdentifiedOligos().end(); ++it)
    {
      digest.push_back(it);
    }

// shorter oligos take (possibly much) less time to process than longer ones;
// due to the sorting order of "NASequence", they also appear earlier in the
// container - therefore use dynamic scheduling to distribute work evenly:
#pragma omp parallel for schedule(dynamic)
    for (SignedSize index = 0; index < SignedSize(digest.size()); ++index)
    {
      IF_MASTERTHREAD
      {
        progresslogger.setProgress(index);
      }

      IdentificationData::IdentifiedOligoRef oligo_ref = digest[index];
      vector<NASequence> all_modified_oligos;
      NASequence ns = oligo_ref->sequence;
      ModifiedNASequenceGenerator::applyFixedModifications(
        fixed_modifications.begin(), fixed_modifications.end(), ns);
      ModifiedNASequenceGenerator::applyVariableModifications(
        variable_modifications.begin(), variable_modifications.end(), ns,
        max_variable_mods_per_oligo, all_modified_oligos, true);

      for (SignedSize mod_idx = 0;
           mod_idx < SignedSize(all_modified_oligos.size()); ++mod_idx)
      {
        const NASequence& candidate = all_modified_oligos[mod_idx];
        double candidate_mass = (use_avg_mass ? candidate.getAverageWeight() :
                                 candidate.getMonoWeight());
        LOG_DEBUG << "Candidate: " << candidate.toString() << " ("
                  << float(candidate_mass) << " Da)" << endl;

        // determine MS2 precursors that match to the current mass
        double tol = search_param.precursor_mass_tolerance;
        if (search_param.precursor_tolerance_ppm)
        {
          tol *= candidate_mass * 1e-6;
        }
        multimap<double, PrecursorInfo>::const_iterator low_it =
          precursor_mass_map.lower_bound(candidate_mass - tol), up_it =
          precursor_mass_map.upper_bound(candidate_mass + tol);

        if (low_it == up_it) continue; // no matching precursor in data

        Int base_charge = negative_mode ? -1 : 1;
        map<Int, PeakSpectrum> theo_spectra_by_charge;
        for (; low_it != up_it; ++low_it)
        {
          LOG_DEBUG << "Matching precursor mass: " << float(low_it->first)
                    << endl;

          Size charge = low_it->second.charge;
          // look up theoretical spectrum for this charge:
          auto pos = theo_spectra_by_charge.find(charge);
          if (pos == theo_spectra_by_charge.end())
          {
            // no spectrum for this charge yet - generate it:
            PeakSpectrum theo_spectrum;
            // @TODO: avoid regenerating spectra with lower charges as part of
            // higher-charged spectra
            spectrum_generator.getSpectrum(theo_spectrum, candidate,
                                           base_charge, charge * base_charge);
            theo_spectrum.setName(candidate.toString());
            pos = theo_spectra_by_charge.insert(make_pair(charge,
                                                          theo_spectrum)).first;
          }
          const PeakSpectrum& theo_spectrum = pos->second;

          Size scan_index = low_it->second.scan_index;
          const PeakSpectrum& exp_spectrum = spectra[scan_index];
          vector<PeptideHit::PeakAnnotation> annotations;
          double score = MetaboliteSpectralMatching::computeHyperScore(
            search_param.fragment_mass_tolerance,
            search_param.fragment_tolerance_ppm, exp_spectrum, theo_spectrum,
            annotations);

          if (!exp_ms2_out.empty())
          {
#pragma omp critical (exp_ms2_out)
            exp_ms2_spectra.addSpectrum(exp_spectrum);
          }
          if (!theo_ms2_out.empty())
          {
#pragma omp critical (theo_ms2_out)
            theo_ms2_spectra.addSpectrum(theo_spectrum);
          }

          if (score < 1e-16) continue; // no hit

          // add oligo hit
          AnnotatedHit ah;
          ah.oligo_ref = oligo_ref;
          ah.mod_index = mod_idx;
          // @TODO: is "observed - calculated" the right way around?
          ah.precursor_error_ppm =
            (low_it->first - candidate_mass) / candidate_mass * 1.0e6;
          ah.annotations = annotations;
          ah.charge = charge;
          ah.adduct = low_it->second.adduct;

#pragma omp atomic
          ++hit_counter;

#ifdef DEBUG_NASEARCH
          cout << "best score in pre-score: " << score << endl;
#endif

#pragma omp critical (annotated_hits_access)
          {
            HitsByScore& scan_hits = annotated_hits[scan_index];
            if ((report_top_hits == 0) || (scan_hits.size() < report_top_hits))
            {
              scan_hits.insert(make_pair(score, ah));
            }
            else // we already have enough hits for this spectrum - replace one?
            {
              double worst_score = (--scan_hits.end())->first;
              if (score >= worst_score)
              {
                scan_hits.insert(make_pair(score, ah));
                // prune list of hits if possible (careful about tied scores!):
                Size n_worst = scan_hits.count(worst_score);
                if (scan_hits.size() - n_worst >= report_top_hits)
                {
                  scan_hits.erase(worst_score);
                }
              }
            }
          }
        }
      }
    }
    progresslogger.endProgress();

    LOG_INFO << "Undigested nucleic acids: " << fasta_db.size() << endl;
    LOG_INFO << "Oligonucleotides: " << id_data.getIdentifiedOligos().size()
             << endl;
    LOG_INFO << "Search hits: " << hit_counter << endl;

    if (!exp_ms2_out.empty())
    {
      MzMLFile().store(exp_ms2_out, exp_ms2_spectra);
    }
    if (!theo_ms2_out.empty())
    {
      MzMLFile().store(theo_ms2_out, theo_ms2_spectra);
    }

    progresslogger.startProgress(0, 1, "post-processing search hits...");
    postProcessHits_(spectra, annotated_hits, id_data, fixed_modifications,
                     variable_modifications, max_variable_mods_per_oligo,
                     negative_mode);
    progresslogger.endProgress();

    // FDR:
    if (!decoy_pattern.empty())
    {
      LOG_INFO << "Performing FDR calculations..." << endl;
      calculateAndFilterFDR_(id_data, report_top_hits == 1);
    }
    id_data.calculateCoverages();

    // store results
    MzTab results = IdentificationDataConverter::exportMzTab(id_data);
    LOG_DEBUG << "Nucleic acid rows: "
              << results.getNucleicAcidSectionRows().size()
              << "\nOligonucleotide rows: "
              << results.getOligonucleotideSectionRows().size()
              << "\nOligo-spectrum match rows: "
              << results.getOSMSectionRows().size() << endl;

    MzTabFile().store(out, results);

    // dummy "peptide" results:
    if (!id_out.empty())
    {
      vector<ProteinIdentification> proteins;
      vector<PeptideIdentification> peptides;
      IdentificationDataConverter::exportIDs(id_data, proteins, peptides, true);
      // proteins[0].setDateTime(DateTime::now());
      // proteins[0].setSearchEngine(toolName_());
      IdXMLFile().store(id_out, proteins, peptides);
    }

    if (!lfq_out.empty())
    {
      generateLFQInput_(id_data, lfq_out);
    }

    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  NucleicAcidSearchEngine tool;
  return tool.main(argc, argv);
}
