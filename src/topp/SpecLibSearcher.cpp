// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: David Wojnar, Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/BinnedSpectrum.h>
#include <OpenMS/COMPARISON/SpectraSTSimilarityScore.h>
#include <OpenMS/COMPARISON/ZhangSimilarityScore.h>
#include <OpenMS/FORMAT/FileHandler.h>
// TODO add ID support to Handler
#include <OpenMS/FORMAT/MSPFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <ctime>
#include <vector>
#include <map>
#include <cmath>
using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_SpecLibSearcher SpecLibSearcher

@brief Identifies peptide MS/MS spectra by spectral matching with a searchable spectral library.

<CENTER>
<table>
    <tr>
        <th ALIGN = "center"> pot. predecessor tools </td>
        <td VALIGN="middle" ROWSPAN=2> &rarr; SpecLibSearcher &rarr;</td>
        <th ALIGN = "center"> pot. successor tools </td>
    </tr>
    <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_SpecLibCreator </td>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter or @n any protein/peptide processing tool</td>
    </tr>
</table>
</CENTER>

@experimental This TOPP-tool is not well tested and not all features might be properly implemented and tested.

@note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_SpecLibSearcher.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_SpecLibSearcher.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPSpecLibSearcher :
  public TOPPBase
{
public:
  TOPPSpecLibSearcher() :
    TOPPBase("SpecLibSearcher", "Identifies peptide MS/MS spectra by spectral matching with a searchable spectral library.")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<files>", ListUtils::create<String>(""), "Input files");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerInputFile_("lib", "<file>", "", "searchable spectral library (MSP format)");
    setValidFormats_("lib", ListUtils::create<String>("msp"));
    registerOutputFileList_("out", "<files>", ListUtils::create<String>(""), "Output files. Have to be as many as input files");
    setValidFormats_("out", ListUtils::create<String>("idXML"));

    registerTOPPSubsection_("precursor", "Precursor (Parent Ion) Options");
    registerDoubleOption_("precursor:mass_tolerance", "<tolerance>", 10.0, "Width of precursor mass tolerance window", false);

    StringList precursor_mass_tolerance_unit_valid_strings;
    precursor_mass_tolerance_unit_valid_strings.push_back("ppm");
    precursor_mass_tolerance_unit_valid_strings.push_back("Da");

    registerStringOption_("precursor:mass_tolerance_unit", "<unit>", "ppm", "Unit of precursor mass tolerance.", false, false);
    setValidStrings_("precursor:mass_tolerance_unit", precursor_mass_tolerance_unit_valid_strings);

    registerIntOption_("precursor:min_charge", "<num>", 2, "Minimum precursor charge to be considered.", false, true);
    registerIntOption_("precursor:max_charge", "<num>", 5, "Maximum precursor charge to be considered.", false, true);

    // consider one before annotated monoisotopic peak and the annotated one
    IntList isotopes = {0, 1};
    registerIntList_("precursor:isotopes", "<num>", isotopes, "Corrects for mono-isotopic peak misassignments. (E.g.: 1 = prec. may be misassigned to first isotopic peak)", false, false);
    
    registerTOPPSubsection_("fragment", "Fragments (Product Ion) Options");
    registerDoubleOption_("fragment:mass_tolerance", "<tolerance>", 10.0, "Fragment mass tolerance", false);
    
//    StringList fragment_mass_tolerance_unit_valid_strings;
//    fragment_mass_tolerance_unit_valid_strings.push_back("ppm");
//    fragment_mass_tolerance_unit_valid_strings.push_back("Da");
    
//    registerStringOption_("fragment:mass_tolerance_unit", "<unit>", "ppm", "Unit of fragment m", false, false);
//    setValidStrings_("fragment:mass_tolerance_unit", fragment_mass_tolerance_unit_valid_strings);

    registerStringOption_("compare_function", "<string>", "ZhangSimilarityScore", "function for similarity comparison", false);
    setValidStrings_("compare_function", {"ZhangSimilarityScore", "SpectraSTSimilarityScore"});

    registerTOPPSubsection_("report", "Reporting Options");
    registerIntOption_("report:top_hits", "<num>", 10, "Maximum number of top scoring hits per spectrum that are reported.", false, true);

    addEmptyLine_();

    registerTOPPSubsection_("filter", "Filtering options. Most are especially useful when the query spectra are raw.");
    registerDoubleOption_("filter:remove_peaks_below_threshold", "<threshold>", 2.01, "All peaks of a query spectrum with intensities below <threshold> will be zeroed.", false);
    registerIntOption_("filter:min_peaks", "<number>", 5, "required minimum number of peaks for a query spectrum", false);
    registerIntOption_("filter:max_peaks", "<number>", 150, "Use only the top <number> of peaks.", false);
    registerIntOption_("filter:cut_peaks_below", "<number>", 1000, "Remove all peaks which are lower than 1/<number> of the highest peaks. Default equals all peaks which are lower than 0.001 of the maximum intensity peak", false);

    registerTOPPSubsection_("modifications", "Modifications Options");
    vector<String> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
    registerStringList_("modifications:fixed", "<mods>", ListUtils::create<String>(""), "Fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)'", false);
    setValidStrings_("modifications:fixed", all_mods);
    registerStringList_("modifications:variable", "<mods>", ListUtils::create<String>(""), "Variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Oxidation (M)'", false);
    setValidStrings_("modifications:variable", all_mods);
    registerIntOption_("modifications:variable_max_per_peptide", "<num>", 2, "Maximum number of residues carrying a variable modification per candidate peptide", false, false);

    addEmptyLine_();
  }

  using MapLibraryPrecursorToLibrarySpectrum = multimap<double, PeakSpectrum>;
    
  MapLibraryPrecursorToLibrarySpectrum annotateIdentificationsToSpectra_(const vector<PeptideIdentification>& ids, 
    const PeakMap& library, 
    StringList variable_modifications, 
    StringList fixed_modifications,
    double remove_peaks_below_threshold)
  {
    MapLibraryPrecursorToLibrarySpectrum annotated_lib;

    ModificationsDB* mdb = ModificationsDB::getInstance();


    // iterate over library spectra and add associated annotations
    PeakMap::const_iterator library_it = library.begin();
    vector<PeptideIdentification>::const_iterator id_it = ids.begin();
    for (; library_it < library.end(); ++library_it, ++id_it)
    {
      const MSSpectrum& lib_spec = *library_it;
      const double& precursor_MZ = lib_spec.getPrecursors()[0].getMZ();

      const PeptideIdentification& id = *id_it;
      const AASequence& aaseq = id.getHits()[0].getSequence();

      PeakSpectrum lib_entry;
      bool variable_modifications_ok(true), fixed_modifications_ok(true);

       // check if each amino acid listed as modified in fixed modifications are modified
       if (!fixed_modifications.empty())
       {
         for (Size j = 0; j < aaseq.size(); ++j)
         {
           const Residue& mod = aaseq.getResidue(j);
           for (Size k = 0; k < fixed_modifications.size(); ++k)
           {
             if (mod.getOneLetterCode()[0] == mdb->getModification(fixed_modifications[k])->getOrigin() && fixed_modifications[k] != mod.getModificationName())
             {
               fixed_modifications_ok = false;
               break;
             }
           }
         }
       }

       // check if each amino acid listed in variable modifications is either unmodified or modified with the corresponding modification
       // Note: this code currently does not allow for multiple variable modifications with same origin
       if (aaseq.isModified() && (!variable_modifications.empty()))
       {
         for (Size j = 0; j < aaseq.size(); ++j)
         {
           if (!aaseq[j].isModified()) { continue; }
          
           const Residue& mod = aaseq.getResidue(j);
           for (Size k = 0; k < variable_modifications.size(); ++k)
           {
             if (mod.getOneLetterCode()[0] == mdb->getModification(variable_modifications[k])->getOrigin() && variable_modifications[k] != mod.getModificationName())
             {
               variable_modifications_ok = false;
               break;
             }
           }
         }
       }

       // TODO: check entries that don't adhere to this rule
       if (!variable_modifications_ok || !fixed_modifications_ok) { continue; }

       // copy peptide identification over to spectrum meta data
       lib_entry.getPeptideIdentifications().push_back(id);
       lib_entry.setPrecursors(lib_spec.getPrecursors());

       // empty array would segfault
       if (id.getHits().empty() || id.getHits()[0].getPeakAnnotations().empty())
       {
         throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Expected StringDataArray of type MSPeakInfo");
       }

       const vector<PeptideHit::PeakAnnotation>& pa = id.getHits()[0].getPeakAnnotations();
       // library entry transformation
       for (UInt l = 0; l < lib_spec.size(); ++l)
       {
         Peak1D peak;
         if (lib_spec[l].getIntensity() > remove_peaks_below_threshold)
         {
           // this is the "MSPPeakInfo" array, see MSPFile which creates a single StringDataArray
           const String& sa = pa[l].annotation;

           // TODO: check why this scaling is done for ? peaks (dubious peaks?)
           if (sa[0] == '?')
           {
             peak.setIntensity(sqrt(0.2 * lib_spec[l].getIntensity()));
           }
           else
           {
             peak.setIntensity(sqrt(lib_spec[l].getIntensity()));
           }

           peak.setMZ(lib_spec[l].getMZ());
           lib_entry.push_back(peak);
         }
       }
       annotated_lib.insert(make_pair(precursor_MZ, lib_entry));
     }
    return annotated_lib;
  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------
    StringList in_spec = getStringList_("in");
    StringList out = getStringList_("out");
    String in_lib = getStringOption_("lib");
    String compare_function = getStringOption_("compare_function");
 
    float precursor_mass_tolerance = getDoubleOption_("precursor:mass_tolerance");
    bool precursor_mass_tolerance_unit_ppm = getStringOption_("precursor:mass_tolerance_unit") == "ppm" ? true : false;

    int pc_min_charge = getIntOption_("precursor:min_charge");
    int pc_max_charge = getIntOption_("precursor:max_charge");

    // consider one before annotated monoisotopic peak and the annotated one
    IntList isotopes = getIntList_("precursor:isotopes");
   
//    float fragment_mass_tolerance = getDoubleOption_("fragment:mass_tolerance");
//    bool fragment_mass_tolerance_unit_ppm = getStringOption_("fragment:mass_tolerance_unit") == "ppm" ? true : false;

    int top_hits = getIntOption_("report:top_hits");

    float remove_peaks_below_threshold = getDoubleOption_("filter:remove_peaks_below_threshold");
    UInt min_peaks = getIntOption_("filter:min_peaks");
    UInt max_peaks = getIntOption_("filter:max_peaks");
    Int cut_peaks_below = getIntOption_("filter:cut_peaks_below");

    StringList fixed_modifications = getStringList_("modifications:fixed");
    StringList variable_modifications = getStringList_("modifications:variable");

    if (top_hits < -1)
    {
      writeLogError_("top_hits (should be  >= -1 )");
      return ILLEGAL_PARAMETERS;
    }

    // -------------------------------------------------------------
    // loading input
    // -------------------------------------------------------------
    if (out.size() != in_spec.size())
    {
      writeLogError_("out (should be as many as input files)");
      return ILLEGAL_PARAMETERS;
    }

    time_t prog_time = time(nullptr);
    MSPFile spectral_library;
    PeakMap query, library;

    time_t start_build_time = time(nullptr);
    // -------------------------------------------------------------
    // building map for faster search
    // -------------------------------------------------------------

    // library containing already identified peptide spectra
    vector<PeptideIdentification> ids;
    spectral_library.load(in_lib, ids, library);

    /*
    // Output bin histogram
    BinnedSpectrum bin_frequency(0.01, 1, PeakSpectrum());
    for (auto const & s : library)
    {
      BinnedSpectrum b(0.01, 1, s);
      // e.g.: bin_frequency.getBins() += b.getBins();  // sum up itensities
      // e.g.: bin_frequency.getBins() += b.getBins().coeffs().cwiseMin(1.0f); // count occupied bins (by truncating intensities >= 1 to 1)
    }

    for (BinnedSpectrum::SparseVectorIteratorType it(bin_frequency.getBins()); it; ++it)
    {
      // output m/z of bin start and average bin intensity
      cout << it.index() * bin_frequency.getBinSize()  << "\t" << static_cast<float>(it.value()/library.size()) << "\n";
      cout << static_cast<float>(it.value()) << "\n";
      cout << static_cast<float>(library.size()) << "\n";
    }
    cout << endl;
    */

    MapLibraryPrecursorToLibrarySpectrum mslib = annotateIdentificationsToSpectra_(ids, library, variable_modifications, fixed_modifications, remove_peaks_below_threshold);

    time_t end_build_time = time(nullptr);
    OPENMS_LOG_INFO << "Time needed for preprocessing data: " << (end_build_time - start_build_time) << "\n";

    //compare function
    std::unique_ptr<PeakSpectrumCompareFunctor> comparator;
    if (compare_function == "SpectraSTSimilarityScore")
    {
      comparator.reset(new SpectraSTSimilarityScore());
    }
    else if (compare_function == "ZhangSimilarityScore")
    {
      comparator.reset(new ZhangSimilarityScore());
    }
    else 
    {
      writeLogError_("Unknown compare function");
      return ILLEGAL_PARAMETERS;
    }
 
   //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    double score;
    StringList::iterator in, out_file;
    for (in  = in_spec.begin(), out_file  = out.begin(); in < in_spec.end(); ++in, ++out_file)
    {
      time_t start_time = time(nullptr);
      FileHandler().loadExperiment(*in, query, {FileTypes::MZML}, log_type_);

      // results
      vector<PeptideIdentification> peptide_ids;
      vector<ProteinIdentification> protein_ids;
      ProteinIdentification prot_id;
 
      //Parameters of identification
      prot_id.setIdentifier("test");
      prot_id.setSearchEngineVersion("SpecLibSearcher");
      prot_id.setDateTime(DateTime::now());
      prot_id.setScoreType(compare_function);

      ProteinIdentification::SearchParameters search_parameters;
      search_parameters.db = getStringOption_("lib");
      search_parameters.charges = String(getIntOption_("precursor:min_charge")) + ":" + String(getIntOption_("precursor:max_charge"));

      ProteinIdentification::PeakMassType mass_type = ProteinIdentification::MONOISOTOPIC;
      search_parameters.mass_type = mass_type;
      search_parameters.fixed_modifications = getStringList_("modifications:fixed");
      search_parameters.variable_modifications = getStringList_("modifications:variable");
      // search_parameters.missed_cleavages = getIntOption_("peptide:missed_cleavages");
      search_parameters.precursor_mass_tolerance = getDoubleOption_("precursor:mass_tolerance");
      search_parameters.precursor_mass_tolerance_ppm = getStringOption_("precursor:mass_tolerance_unit") == "ppm" ? true : false;
//      search_parameters.fragment_mass_tolerance = getDoubleOption_("fragment:mass_tolerance");
//      search_parameters.fragment_mass_tolerance_ppm = getStringOption_("fragment:mass_tolerance_unit") == "ppm" ? true : false;

//TODO: report an Enzyme?

      prot_id.setSearchParameters(search_parameters);


      /***********SEARCH**********/
      for (UInt j = 0; j < query.size(); ++j)
      {
        //Set identifier for each identifications
        PeptideIdentification pid;
        pid.setIdentifier("test");
        pid.setScoreType(compare_function);
        ProteinHit pr_hit;
        pr_hit.setAccession(j);
        prot_id.insertHit(pr_hit);

        // proper MS2?
        if (query[j].empty() || query[j].getMSLevel() != 2)
        {
          continue;
        }

        if (query[j].getPrecursors().empty())
        {
          writeLogWarn_("Warning MS2 spectrum without precursor information");
          continue;
        }

        // filter query spectrum
        double max_intensity = std::max_element(query[j].begin(), query[j].end(), 
                                [](const Peak1D& l, const Peak1D& r) 
                                { 
                                  return (l.getIntensity() < r.getIntensity()); 
                                })->getIntensity();

        double min_high_intensity = max_intensity / cut_peaks_below;

        PeakSpectrum filtered_query;
        for (UInt k = 0; k < query[j].size(); ++k)
        {
          if (query[j][k].getIntensity() >= remove_peaks_below_threshold 
           && query[j][k].getIntensity() >= min_high_intensity)
          {
            Peak1D peak;
            peak.setIntensity(sqrt(query[j][k].getIntensity()));
            peak.setMZ(query[j][k].getMZ());
            filtered_query.push_back(peak);
          }
        }

        // retain only top N peaks
        if (filtered_query.size() > max_peaks)
        {
          filtered_query.sortByIntensity(true);
          filtered_query.resize(max_peaks);
          filtered_query.sortByPosition();
        }

        if (filtered_query.size() < min_peaks)
        { 
          continue;
        }

        const double& query_rt = query[j].getRT();
        const int& query_charge = query[j].getPrecursors()[0].getCharge();
        const double query_mz = query[j].getPrecursors()[0].getMZ();
        
        if (query_charge > 0 && (query_charge < pc_min_charge || query_charge > pc_max_charge))
        { 
          continue;
        } 

        for (auto const & iso : isotopes)
        {
          // isotopic misassignment corrected query
          const double ic_query_mz = query_mz - iso * Constants::C13C12_MASSDIFF_U;

          // if tolerance unit is ppm convert to m/z
          const double precursor_mass_tolerance_mz = precursor_mass_tolerance_unit_ppm ? ic_query_mz * precursor_mass_tolerance * 1e-6 : precursor_mass_tolerance;

          // skip matching of isotopic misassignments if charge not annotated
          if (iso != 0 && query_charge == 0)
          {
            continue;
          }

          // skip matching of isotopic misassignments if search windows around isotopic peaks would overlap (resulting in more than one report of the same hit)
          const double isotopic_peak_distance_mz = Constants::C13C12_MASSDIFF_U / query_charge;
          if (iso != 0 && precursor_mass_tolerance_mz >= 0.5 * isotopic_peak_distance_mz)
          { 
            continue;
          }

          /* TODO: remove old code for charge estimation?
          bool charge_one = false;
          Int percent = (Int) Math::round((query[j].size() / 100.0) * 3.0);
          Int margin  = (Int) Math::round((query[j].size() / 100.0) * 1.0);
          for (vector<Peak1D>::iterator peak = query[j].end() - 1; percent >= 0; --peak, --percent)
          {
            if (peak->getMZ() < query_MZ)
            {
              break;
            }
          }
          if (percent > margin)
          {
            charge_one = true;
          }
          */


          // determine MS2 precursors that match to the current peptide mass
          MapLibraryPrecursorToLibrarySpectrum::const_iterator low_it, up_it;
        
          low_it = mslib.lower_bound(ic_query_mz - 0.5 * precursor_mass_tolerance_mz);
          up_it = mslib.upper_bound(ic_query_mz + 0.5 * precursor_mass_tolerance_mz);
        
          // no matching precursor in data
          if (low_it == up_it)
          { 
            continue;
          }
       
          for (; low_it != up_it; ++low_it)
          {
            const PeakSpectrum& lib_spec = low_it->second;;
            PeptideHit hit = lib_spec.getPeptideIdentifications()[0].getHits()[0];
            const int& lib_charge = hit.getCharge();  

            // check if charge state between library and experimental spectrum match
            if (query_charge > 0 && lib_charge != query_charge)
            {
              continue;
            }

            // Special treatment for SpectraST score as it computes a score based on the whole library
            if (compare_function == "SpectraSTSimilarityScore")
            {
              auto& sp = dynamic_cast<SpectraSTSimilarityScore&>(*comparator);
              BinnedSpectrum quer_bin_spec = sp.transform(filtered_query);
              BinnedSpectrum lib_bin_spec = sp.transform(lib_spec);
              score = sp(filtered_query, lib_spec); //(*sp)(quer_bin,librar_bin);
              double dot_bias = sp.dot_bias(quer_bin_spec, lib_bin_spec, score);
              hit.setMetaValue("DOTBIAS", dot_bias);
            }
            else
            {
              score = (*comparator)(filtered_query, lib_spec);
            }

            DataValue RT(lib_spec.getRT());
            DataValue MZ(lib_spec.getPrecursors()[0].getMZ());
            hit.setMetaValue("lib:RT", RT);
            hit.setMetaValue("lib:MZ", MZ);
            hit.setMetaValue(Constants::UserParam::ISOTOPE_ERROR, iso);
            hit.setScore(score);
            PeptideEvidence pe;
            pe.setProteinAccession(pr_hit.getAccession());
            hit.addPeptideEvidence(pe);
            pid.insertHit(hit);
          }
        }

        pid.setHigherScoreBetter(true);
        pid.sort();

        if (compare_function == "SpectraSTSimilarityScore")
        {
          if (!pid.empty() && !pid.getHits().empty())
          {
            vector<PeptideHit> final_hits;
            final_hits.resize(pid.getHits().size());
            auto& sp = dynamic_cast<SpectraSTSimilarityScore&>(*comparator);
            Size runner_up = 1;
            for (; runner_up < pid.getHits().size(); ++runner_up)
            {
              if (pid.getHits()[0].getSequence().toUnmodifiedString() != pid.getHits()[runner_up].getSequence().toUnmodifiedString() 
               || runner_up > 5)
              {
                break;
              }
            }
            double delta_D = sp.delta_D(pid.getHits()[0].getScore(), pid.getHits()[runner_up].getScore());
            for (Size s = 0; s < pid.getHits().size(); ++s)
            {
              final_hits[s] = pid.getHits()[s];
              final_hits[s].setMetaValue("delta D", delta_D);
              final_hits[s].setMetaValue("dot product", pid.getHits()[s].getScore());
              final_hits[s].setScore(sp.compute_F(pid.getHits()[s].getScore(), delta_D, pid.getHits()[s].getMetaValue("DOTBIAS")));
            }
            pid.setHits(final_hits);
            pid.sort();
            pid.setMZ(query[j].getPrecursors()[0].getMZ());
            pid.setRT(query_rt);
          }
        }

        if (top_hits != -1 && (UInt)top_hits < pid.getHits().size())
        {
          pid.getHits().resize(top_hits);
        }
        peptide_ids.push_back(pid);
      }
      protein_ids.push_back(prot_id);

      //-------------------------------------------------------------
      // writing output
      //-------------------------------------------------------------
      FileHandler().storeIdentifications(*out_file, protein_ids, peptide_ids, {FileTypes::IDXML});
      time_t end_time = time(nullptr);
      OPENMS_LOG_INFO << "Search time: " << difftime(end_time, start_time) << " seconds for " << *in << "\n";
    }
    time_t end_time = time(nullptr);
    OPENMS_LOG_INFO << "Total time: " << difftime(end_time, prog_time) << " seconds\n";
    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPSpecLibSearcher tool;
  return tool.main(argc, argv);
}

/// @endcond
