// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/SearchEngineBase.h>

#include <OpenMS/ANALYSIS/ID/PeptideIndexing.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/PepXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <OpenMS/FORMAT/PercolatorInfile.h>
#include <OpenMS/FORMAT/HANDLERS/IndexedMzMLDecoder.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>

#include <OpenMS/SYSTEM/File.h>

#include <fstream>
#include <regex>

#include <QStringList>
#include <chrono>
#include <map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>

#include <boost/math/distributions/normal.hpp>

using namespace OpenMS;
using namespace std;
using boost::math::normal;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_SageAdapter SageAdapter

@brief Identifies peptides in MS/MS spectra via sage.

<CENTER>
    <table>
        <tr>
            <th ALIGN = "center"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> &rarr; SageAdapter &rarr;</td>
            <th ALIGN = "center"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any signal-/preprocessing tool @n (in mzML format)</td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter or @n any protein/peptide processing tool</td>
        </tr>
    </table>
</CENTER>

@em Sage must be installed before this wrapper can be used.

Only the closed-search identification mode of Sage is supported by this adapter.
Currently, also neither "wide window" (= open or DIA) mode, nor "chimeric" mode is supported,
because of limitations in OpenMS' data structures and file formats.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_SageAdapter.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_SageAdapter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


#define CHRONOSET

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif



class TOPPSageAdapter :
  public SearchEngineBase
{
public: 
  TOPPSageAdapter() :
    SearchEngineBase("SageAdapter", "Annotates MS/MS spectra using Sage.", true,
             {
                 {"Michael Lazear",
                 "Sage: An Open-Source Tool for Fast Proteomics Searching and Quantification at Scale",
                 "J. Proteome Res. 2023, 22, 11, 3652–3659",
                 "https://doi.org/10.1021/acs.jproteome.3c00486"}
             })
  {
  }

  // Saves details of PTMs as well, useful for if more than one PTM is mapped to a given mass 
  struct modification
  {
    double count = 0; 
    vector<double> mass; 
    int numcharges = 0; 
  }; 

    // Define a struct to hold modification data
  struct ModData
  {
      int count;            // Modification rate
      String name;            // Modification name
      int numcharges;      // Number of charges
      vector<double> masses;  // Masses associated with the modification
  };

  // Comparator for approximate comparison of double values
  struct FuzzyDoubleComparator {
      double epsilon;
      FuzzyDoubleComparator(double eps = 1e-9) : epsilon(eps) {}
      bool operator()(const double& a, const double& b) const 
      {
          return std::fabs(a - b) >= epsilon && a < b;
      }
  };

// delta mass counts to delta masses
//typedef map<double, double, FuzzyDoubleComparator> CountToDeltaMass;

typedef map<double, double, FuzzyDoubleComparator> DeltaMassHistogram; // maps delta mass to count
typedef map<double, int, FuzzyDoubleComparator> DeltaMasstoCharge; // maps delta mass to count

// Gaussian function
static double gaussian(double x, double sigma) {
    return exp(-(x*x) / (2 * sigma*sigma)) / (sigma * sqrt(2 * M_PI));
}

// Smooths the PTM-mass histogram , uses a Kernel Density Estimation on top of the histogram. 
// Smooths the PTM-mass histogram using Gaussian Kernel Density Estimation (KDE).
static DeltaMassHistogram smoothDeltaMassHist(const DeltaMassHistogram& hist, double sigma = 0.001)
{
    if (hist.size() < 3)
    {
      return hist; //Not enough data points for smoothing 
    }
    // Create a smoothed histogram with a fuzzy comparator for floating-point keys
    DeltaMassHistogram smoothed_hist(FuzzyDoubleComparator(1e-9));

    // Extract delta masses and counts into vectors for efficient access
    std::vector<double> deltas;
    std::vector<double> counts;
    deltas.reserve(hist.size());
    counts.reserve(hist.size());

    for (const auto& [delta, count] : hist)
    {
        deltas.push_back(delta);
        counts.push_back(count);
    }

    const size_t n = deltas.size();
    std::vector<double> smoothed_counts(n, 0.0);

    // Perform Gaussian smoothing
    for (size_t i = 0; i < n; ++i)
    {
        double weight_sum = 0.0;

        for (size_t j = 0; j < n; ++j)
        {
            double mz_diff = deltas[i] - deltas[j];

            // Ignore points beyond 3 standard deviations
            if (std::abs(mz_diff) > 3.0 * sigma)
                continue;

            double weight = gaussian(mz_diff, sigma);
            smoothed_counts[i] += weight * counts[j];
            weight_sum += weight;
        }

        if (weight_sum != 0.0)
        {
            smoothed_counts[i] /= weight_sum;
        }
    }

    // Populate the smoothed histogram
    for (size_t i = 0; i < n; ++i)
    {
        smoothed_hist[deltas[i]] = smoothed_counts[i];
    }

    return smoothed_hist;
}

// Identifies local maxima in the delta mass histogram based on count threshold and SNR.
static DeltaMassHistogram findPeaksInDeltaMassHistogram(const DeltaMassHistogram& hist, double count_threshold = 0.0, double SNR = 2.0)
{
    if (hist.size() < 3)
    {
        return hist;  // Not enough data points to find peaks
    }

    DeltaMassHistogram peaks(FuzzyDoubleComparator(1e-9));

    // Extract counts to compute noise level (median count)
    std::vector<double> counts;
    counts.reserve(hist.size());

    for (const auto& [_, count] : hist)
    {
        counts.push_back(count);
    }

    // Calculate median as noise level
    std::nth_element(counts.begin(), counts.begin() + counts.size() / 2, counts.end());
    double noise_level = counts[counts.size() / 2];

    // Convert histogram to vector for indexed access
    std::vector<std::pair<double, double>> hist_vector(hist.begin(), hist.end());

    // Check each point except the first and last for local maxima
    for (size_t i = 1; i < hist_vector.size() - 1; ++i)
    {
        double prev_count = hist_vector[i - 1].second;
        double curr_count = hist_vector[i].second;
        double next_count = hist_vector[i + 1].second;

        // Check if current point is a local maximum
        if (curr_count >= prev_count && curr_count >= next_count &&
            curr_count > count_threshold &&
            curr_count / noise_level > SNR)
        {
            peaks[hist_vector[i].first] = curr_count;
        }
    }

    return peaks;
}

// Returns the maxima of a histogram from the delta masses of each peptide.
std::pair<DeltaMassHistogram, DeltaMasstoCharge> getDeltaClusterCenter(const std::vector<PeptideIdentification>& pips, bool smoothing = false, bool debug = false)
{
    // Constants
    constexpr double deltamass_tolerance = 0.0005;
    constexpr double delta_mass_zero_treshold = 0.05;

    // Lambda to round values to the specified tolerance
    auto roundToTolerance = [deltamass_tolerance](double value) {
        return std::round(value / deltamass_tolerance) * deltamass_tolerance;
    };

    // Data structures to store histogram and charge states
    DeltaMassHistogram hist(FuzzyDoubleComparator(1e-9));
    DeltaMasstoCharge num_charges_at_mass(FuzzyDoubleComparator(1e-9));
    std::unordered_map<double, std::unordered_set<int>> charge_states;

    // Process each peptide identification
    for (const auto& id : pips)
    {
        const auto& hits = id.getHits();
        for (const auto& hit : hits)
        {
            // Retrieve delta mass and charge
            double delta_mass = hit.getMetaValue("DeltaMass");
            int charge = hit.getCharge();

            // Ignore delta masses close to zero
            if (std::abs(delta_mass) <= delta_mass_zero_treshold)
                continue;

            // Round delta mass to bin similar values
            double rounded_mass = roundToTolerance(delta_mass);

            // Update histogram count
            hist[rounded_mass] += 1.0;

            // Update unique charge count
            if (charge_states[rounded_mass].insert(charge).second)
            {
                num_charges_at_mass[rounded_mass] += 1.0;
            }
        }
    }

    // Prepare results
    std::pair<DeltaMassHistogram, DeltaMasstoCharge> results;   
    results = { hist, num_charges_at_mass };

    // Apply smoothing if requested
    if (smoothing)
    {
        DeltaMassHistogram smoothed_hist = smoothDeltaMassHist(hist, 0.0001);
        DeltaMassHistogram hist_maxima = findPeaksInDeltaMassHistogram(smoothed_hist, 0.0, 3.0);

        // Update charge counts for the smoothed maxima
        DeltaMasstoCharge num_charges_at_mass_smoothed(FuzzyDoubleComparator(1e-9));
        for (const auto& [mass, _] : hist_maxima)
        {
            num_charges_at_mass_smoothed[mass] = num_charges_at_mass[mass];
        }

        // Update results with smoothed data
        results = { hist_maxima, num_charges_at_mass_smoothed };
    }

    return results;
}

//Fucntion that maps a selection of masses to certain PTMs and returns a summary of said PTMs. Also adds PTM for each petide without in-peptide localization. 
vector<PeptideIdentification> mapDifftoMods(DeltaMassHistogram hist, DeltaMasstoCharge charge_hist, vector<PeptideIdentification>& pips, double precursor_mass_tolerance_ = 5, bool precursor_mass_tolerance_unit_ppm = true, String outfile = "")
{
  vector<vector<PeptideIdentification>> clusters(hist.size(), vector<PeptideIdentification>());
  map<double, String, FuzzyDoubleComparator> mass_of_mods(FuzzyDoubleComparator(1e-9));
  vector<pair<double, String>> mass_of_mods_vec;

  // Load modifications from the database
  vector<String> searchmodifications_names;
  ModificationsDB* mod_db = ModificationsDB::getInstance();
  mod_db->getAllSearchModifications(searchmodifications_names);
  for (const String& m : searchmodifications_names)
  {
      const ResidueModification* residue = mod_db->getModification(m);
      String res_name = residue->getFullName();
      double res_diffmonoMass = residue->getDiffMonoMass();
      if (res_name.find("substitution") == string::npos)
          mass_of_mods[res_diffmonoMass] = res_name;
  }

  // Generate combinations of modifications
  map<double, String, FuzzyDoubleComparator> combo_mods(FuzzyDoubleComparator(1e-9));
  for (auto mit = mass_of_mods.begin(); mit != mass_of_mods.end(); ++mit)
  {
      for (auto mit2 = mit; mit2 != mass_of_mods.end(); ++mit2)
      {
          combo_mods[mit->first + mit2->first] = mit->second + "++" + mit2->second;
      }
  }

  // Variables for mapping
  StringList modnames;
  map<String, modification> modifications;
  map<double, String> hist_found;

  // Helper function to add or update modifications
  auto addOrUpdateModification = [&](const String& mod_name, double mass, double count, int numcharges)
  {
      if (modifications.find(mod_name) == modifications.end())
      {
          modification modi{};
          modi.mass.push_back(mass);
          modi.count = count;
          modi.numcharges = numcharges;
          modifications[mod_name] = modi;
      }
      else
      {
          modifications[mod_name].count += count;
          modifications[mod_name].numcharges = max(numcharges, modifications[mod_name].numcharges);
      }
  };

  // Mapping with tolerances //TODO: fix code again, add back high_it 
  for (const auto& hist_entry : hist)
  {
    //Values from the histogram 
    double current_cluster_mass = hist_entry.first;
    double count = hist_entry.second;
    
    double lowerbound, upperbound; 

    const double epsilon = 1e-8;

    if (precursor_mass_tolerance_unit_ppm) // ppm
    {
        double tolerance = current_cluster_mass * precursor_mass_tolerance_ * 1e-6;
        lowerbound = current_cluster_mass - tolerance;
        upperbound = current_cluster_mass + tolerance;
    }
    else // Dalton
    {
        lowerbound = current_cluster_mass - precursor_mass_tolerance_;
        upperbound = current_cluster_mass + precursor_mass_tolerance_;
    }

    // Search for modifications within bounds
    bool mapping_found = false;
    String mod_name;
    double mod_mass = 0.0;

    // Search in single modifications using lower_bound
    auto it_lower = mass_of_mods.lower_bound(lowerbound - epsilon);
    bool found_lower = false;
    if (it_lower != mass_of_mods.end() && fabs(it_lower->first - current_cluster_mass) <= precursor_mass_tolerance_)
    {
        found_lower = true;
    }

    // Search in single modifications using upper_bound
    auto it_upper = mass_of_mods.upper_bound(upperbound + epsilon);
    bool found_upper = false;
    if (it_upper != mass_of_mods.begin())
    {
        --it_upper; // Move to the largest element <= upperbound
        if (fabs(it_upper->first - current_cluster_mass) <= precursor_mass_tolerance_)
        {
            found_upper = true;
        }
    }

    // Compare results from lower_bound and upper_bound
    if (found_lower && found_upper)
    {
        if (it_lower->first == it_upper->first && it_lower->second == it_upper->second)
        {
            // Both methods found the same modification
            mod_name = it_lower->second;
            mod_mass = it_lower->first;
            hist_found[mod_mass] = mod_name;
            mapping_found = true;
        }
        else
        {
            // Different results from lower_bound and upper_bound
            // Choose the closer one
            mod_name = it_lower->second + "//" + it_upper->second; 
            mod_mass = current_cluster_mass;
            hist_found[it_lower->first] = it_lower->second;
            hist_found[it_upper->first] = it_upper->second;
            mapping_found = true;
        }
    }
    else
    {
      // Check if modification can be explained by known modifications
      for (const auto& hit : hist_found)
      {
          if (fabs(hit.first - current_cluster_mass) < precursor_mass_tolerance_)
          {
              addOrUpdateModification(hit.second, hit.first, count, charge_hist[current_cluster_mass]); 
              mapping_found = true;
              break;
          } // Check if modification can be explained by a +1 Isotope variant of a known modification 
          else if (fabs((hit.first + 1) - current_cluster_mass) < precursor_mass_tolerance_)
          {
              String temp_mod_name = hit.second + "+1Da";
              addOrUpdateModification(temp_mod_name, hit.first + 1, count, charge_hist[current_cluster_mass]);
              hist_found[hit.first + 1] = temp_mod_name;
              mapping_found = true;
              break;
          }
      }
      // Search in combination modifications
      if (!mapping_found)
      {
        auto it = combo_mods.lower_bound(current_cluster_mass - epsilon);
        if (it != combo_mods.end() && fabs(it->first - current_cluster_mass) <= precursor_mass_tolerance_ / 10)
        {
            mod_name = it->second;
            mod_mass = it->first;
            mapping_found = true;
        }
      }
    }
    if (fabs(mod_mass) <  precursor_mass_tolerance_) continue; //If the closest mod_mass is too close to 0, continue

    if (mapping_found)
    {
        modnames.push_back(mod_name);
        addOrUpdateModification(mod_name, mod_mass, count, charge_hist[current_cluster_mass]);
    }
    else
    {
        // Unknown modification
        String unknown_mod_name = "Unknown" + std::to_string(std::round(current_cluster_mass));
        addOrUpdateModification(unknown_mod_name, current_cluster_mass, count, charge_hist[current_cluster_mass]);
    }
  }

  // Collect all modification data into a vector
  vector<ModData> mods_by_count;

  //Fill vetcor 
  for (const auto& mod_pair : modifications)
  {
      ModData mod_data;
      mod_data.count = std::round(mod_pair.second.count);
      mod_data.name = mod_pair.first;
      mod_data.numcharges = mod_pair.second.numcharges;
      mod_data.masses = mod_pair.second.mass;

      mods_by_count.push_back(mod_data);
  }

  // Sort the modifications based on (numcharges + rate) in descending order
  sort(mods_by_count.begin(), mods_by_count.end(),
      [](const ModData& a, const ModData& b)
      {
          return (a.numcharges + a.count) > (b.numcharges + b.count);
      });

  // Add the modifications to the output for each peptide
  for (auto& id : pips)
  {
      auto& hits = id.getHits();
      for (auto& h : hits)
      {
          double deltamass = h.getMetaValue("DeltaMass");
          String PTM = "";

          // Check if too close to zero 
          if (fabs(deltamass) < 0.05)
          {
              h.setMetaValue("PTM", PTM);
              continue;
          }

          bool found = false;
          // Check with error tolerance if already present in histogram
          for (const auto& mit : hist_found)
          {
              if (fabs(deltamass - mit.first) < precursor_mass_tolerance_)
              {
                  PTM = mit.second;
                  found = true;
                  break;
              }
          }
          //Otherwise assign unkwown 
          if (!found)
          {
              PTM = "Unknown" + String(deltamass);
          }
          h.setMetaValue("PTM", PTM);
      }
  }
  // Remove 'idxml' from output file name and write the table
  String output_tab = outfile.substr(0, outfile.size() - 5) + "_OutputTable.tsv";
  std::ofstream outfile_stream(output_tab);

  // Check if the file was opened successfully
  if (!outfile_stream.is_open())
  {
      std::cerr << "Error opening file: " << output_tab << std::endl;
      // Handle the error appropriately, e.g., return or exit
      return pips; // Assuming pips is the default return value
  }

  outfile_stream << "Name\tMass\tModified Peptides (incl. charge variants)\tModified Peptides\n";

  // Iterate over the data and write to the file
  for (const auto& mod_data : mods_by_count)
  {
      outfile_stream << mod_data.name << '\t';

      // Output mass or masses
      if (mod_data.masses.size() < 2)
      {
          outfile_stream << mod_data.masses.at(0) << '\t';
      }
      else
      {
          outfile_stream << mod_data.masses.at(0) << "/" << mod_data.masses.at(1) << '\t';
      }

      // Output rounded values
      outfile_stream << mod_data.numcharges + mod_data.count << '\t'
                    << mod_data.count << '\n';
  }

  // Close the file
  outfile_stream.close();

  //Return the peptides with the additional PTM column 
  return pips; 
} 


protected:
  // create a template-based configuration file for sage
  // variable values correspond to sage parameter that can be configured via TOPP tool parameter.
  // values will be pasted into the config_template at the corresponding tag. E.g. bucket_size at tag ##bucket_size##
  static constexpr size_t bucket_size = 8192;
  static constexpr size_t min_len = 5; 
  static constexpr size_t max_len = 50; 
  static constexpr size_t missed_cleavages = 2;
  static constexpr double fragment_min_mz = 200.0;
  static constexpr double fragment_max_mz = 2000.0;
  static constexpr double peptide_min_mass = 500.0;
  static constexpr double peptide_max_mass = 5000.0;
  static constexpr size_t min_ion_index = 2;
  static constexpr size_t max_variable_mods = 2;
  const std::string precursor_tol_unit = "ppm";
  static constexpr double precursor_tol_left = -6.0;
  static constexpr double precursor_tol_right = 6.0;
  const std::string fragment_tol_unit = "ppm";
  static constexpr double fragment_tol_left = -10.0;
  static constexpr double fragment_tol_right = 10.0;
  const std::string isotope_errors = "-1, 3";
  const std::string charges_if_not_annotated = "2, 5";
  static constexpr size_t min_matched_peaks = 6;
  static constexpr size_t report_psms = 1;
  static constexpr size_t min_peaks = 15;
  static constexpr size_t max_peaks = 150;

  std::string config_template = R"(
{
  "database": {
    "bucket_size": ##bucket_size##,
    "enzyme": {
      "missed_cleavages": ##missed_cleavages##,
      "min_len": ##min_len##,
      "max_len": ##max_len##,
      ##enzyme_details##
    },
    "fragment_min_mz": ##fragment_min_mz##,
    "fragment_max_mz": ##fragment_max_mz##,
    "peptide_min_mass": ##peptide_min_mass##,
    "peptide_max_mass": ##peptide_max_mass##,
    "ion_kinds": ["b", "y"],
    "min_ion_index": ##min_ion_index##,
    "static_mods": {
      ##static_mods##
    },
    "variable_mods": {
      ##variable_mods##
    },
    "max_variable_mods": ##max_variable_mods##,
    "generate_decoys": false,
    "decoy_tag": "##decoy_prefix##"
  },
  "precursor_tol": {
    "##precursor_tol_unit##": [
      ##precursor_tol_left##,
      ##precursor_tol_right##
    ]
  },
  "fragment_tol": {
    "##fragment_tol_unit##": [
    ##fragment_tol_left##,
    ##fragment_tol_right##
    ]
  },
  "precursor_charge": [
    ##charges_if_not_annotated##
  ],
  "isotope_errors": [
    ##isotope_errors##
  ],
  "deisotope": ##deisotope##,
  "chimera": ##chimera##,
  "predict_rt": ##predict_rt##,
  "min_peaks": ##min_peaks##,
  "max_peaks": ##max_peaks##,
  "min_matched_peaks": ##min_matched_peaks##,
  "report_psms": ##report_psms##, 
  "wide_window": ##wide_window##
}
)";

  // formats a single mod entry as sage json entry
  String getModDetails(const ResidueModification* mod, const Residue* res)
  {
    String origin;
    if (mod->getTermSpecificity() == ResidueModification::N_TERM)
    { 
      origin += "^";
    }
    else if (mod->getTermSpecificity() == ResidueModification::C_TERM)
    {
      origin += "$";
    }
    else if (mod->getTermSpecificity() == ResidueModification::PROTEIN_N_TERM)
    {
      origin += "[";
    }
    else if (mod->getTermSpecificity() == ResidueModification::PROTEIN_C_TERM)
    {
      origin += "]";
    }
   if (res != nullptr && res->getOneLetterCode() != "X") // omit letter for "any AA"
   {
     origin += res->getOneLetterCode();
   }

    return String("\"") + origin + "\": " + String(mod->getDiffMonoMass());
  }

  // formats all mod entries into a single multi-line json string
  String getModDetailsString(const OpenMS::ModifiedPeptideGenerator::MapToResidueType& mod_map)
  {
    String mod_details;
    for (auto it = mod_map.val.begin(); it != mod_map.val.end(); ++it)
    {
      const auto& mod = it->first;
      const auto& res = it->second;
      mod_details += getModDetails(mod, res);     
      if (std::next(it) != mod_map.val.end())
      {
        mod_details += ",\n";
      }
    }
    return mod_details;
  }

  // impute values into config_template
  // TODO just iterate over all options??
  String imputeConfigIntoTemplate()
  {
    String config_file = config_template;
    config_file.substitute("##bucket_size##", String(getIntOption_("bucket_size")));
    config_file.substitute("##min_len##", String(getIntOption_("min_len")));
    config_file.substitute("##max_len##", String(getIntOption_("max_len")));
    config_file.substitute("##missed_cleavages##", String(getIntOption_("missed_cleavages")));
    config_file.substitute("##fragment_min_mz##", String(getDoubleOption_("fragment_min_mz")));
    config_file.substitute("##fragment_max_mz##", String(getDoubleOption_("fragment_max_mz")));
    config_file.substitute("##peptide_min_mass##", String(getDoubleOption_("peptide_min_mass")));
    config_file.substitute("##peptide_max_mass##", String(getDoubleOption_("peptide_max_mass")));
    config_file.substitute("##min_ion_index##", String(getIntOption_("min_ion_index")));
    config_file.substitute("##max_variable_mods##", String(getIntOption_("max_variable_mods")));
    config_file.substitute("##precursor_tol_unit##", getStringOption_("precursor_tol_unit") == "Da" ? "da" : "ppm"); // sage might expect lower-case "da"
    config_file.substitute("##precursor_tol_left##", String(getDoubleOption_("precursor_tol_left")));
    config_file.substitute("##precursor_tol_right##", String(getDoubleOption_("precursor_tol_right")));
    config_file.substitute("##fragment_tol_unit##", getStringOption_("fragment_tol_unit") == "Da" ? "da" : "ppm"); // sage might expect lower-case "da"
    config_file.substitute("##fragment_tol_left##", String(getDoubleOption_("fragment_tol_left")));
    config_file.substitute("##fragment_tol_right##", String(getDoubleOption_("fragment_tol_right")));
    config_file.substitute("##isotope_errors##", getStringOption_("isotope_error_range"));
    config_file.substitute("##charges_if_not_annotated##", getStringOption_("charges"));
    config_file.substitute("##min_matched_peaks##", String(getIntOption_("min_matched_peaks")));
    config_file.substitute("##min_peaks##", String(getIntOption_("min_peaks")));
    config_file.substitute("##max_peaks##", String(getIntOption_("max_peaks")));
    config_file.substitute("##report_psms##", String(getIntOption_("report_psms")));
    config_file.substitute("##deisotope##", getStringOption_("deisotope")); 
    config_file.substitute("##chimera##", getStringOption_("chimera")); 
    config_file.substitute("##predict_rt##", getStringOption_("predict_rt")); 
    config_file.substitute("##decoy_prefix##", getStringOption_("decoy_prefix")); 
    config_file.substitute("##wide_window##", getStringOption_("wide_window")); 

    
    //Look at decoy handling 

    String enzyme = getStringOption_("enzyme");
    String enzyme_details;
    if (enzyme == "Trypsin")
    {
      enzyme_details = 
   R"("cleave_at": "KR",
      "restrict": "P",
      "c_terminal": true)";
    }
    else if (enzyme == "Trypsin/P")
    {
      enzyme_details = 
   R"("cleave_at": "KR",
      "restrict": null,
      "c_terminal": true)";
    }
    else if (enzyme == "Chymotrypsin")
    {
      enzyme_details = 
   R"("cleave_at": "FWYL",
      "restrict": "P",
      "c_terminal": true)";
    }
    else if (enzyme == "Chymotrypsin/P")
    {
      enzyme_details = 
   R"("cleave_at": "FWYL",
      "restrict": null,
      "c_terminal": true)";
    }
    else if (enzyme == "Arg-C")
    {
      enzyme_details = 
   R"("cleave_at": "R",
      "restrict": "P",
      "c_terminal": true)";
    }
    else if (enzyme == "Arg-C/P")
    {
      enzyme_details = 
   R"("cleave_at": "R",
      "restrict": null,
      "c_terminal": true)";
    }
    else if (enzyme == "Lys-C")
    {
      enzyme_details = 
   R"("cleave_at": "K",
      "restrict": "P",
      "c_terminal": true)";
    }
    else if (enzyme == "Lys-C/P")
    {
      enzyme_details = 
   R"("cleave_at": "K",
      "restrict": null,
      "c_terminal": true)";
    }    
    else if (enzyme == "Lys-N")
    {
      enzyme_details = 
   R"("cleave_at": "K",
      "restrict": null,
      "c_terminal": false)";
    }
    else if (enzyme == "no cleavage")
    {
      enzyme_details = 
   R"("cleave_at": "$")";
    }    
    else if (enzyme == "unspecific cleavage")
    {
      enzyme_details = 
   R"("cleave_at": "")";
    }
    else if (enzyme == "glutamyl endopeptidase")
    {
      enzyme_details =
   R"("cleave_at": "E",
      "restrict": "E",
      "c_terminal":true)";
    }
    else if (enzyme == "leukocyte elastase")
    {
      enzyme_details =
   R"("cleave_at": "ALIV",
      "restrict": null,
      "c_terminal":true)";
    }

    config_file.substitute("##enzyme_details##", enzyme_details);

    
    auto fixed_mods = getStringList_("fixed_modifications");
    set<String> fixed_unique(fixed_mods.begin(), fixed_mods.end());
    fixed_mods.assign(fixed_unique.begin(), fixed_unique.end());   
    ModifiedPeptideGenerator::MapToResidueType fixed_mod_map = ModifiedPeptideGenerator::getModifications(fixed_mods); // std::unordered_map<const ResidueModification*, const Residue*> val;
    String static_mods_details = getModDetailsString(fixed_mod_map);

    auto variable_mods = getStringList_("variable_modifications");
    set<String> variable_unique(variable_mods.begin(), variable_mods.end());
    variable_mods.assign(variable_unique.begin(), variable_unique.end());
    ModifiedPeptideGenerator::MapToResidueType variable_mod_map = ModifiedPeptideGenerator::getModifications(variable_mods);
    String variable_mods_details = getModDetailsString(variable_mod_map);

    //Treat variables as list for sage v0.15 and beyond 
    StringList static_mods_details_list; 
    StringList variable_mods_details_list; 

    String static_mods_details_split = static_mods_details; 
    String variable_mods_details_split = variable_mods_details; 
    static_mods_details_split.split(",", static_mods_details_list); 
    variable_mods_details_split.split(",", variable_mods_details_list); 

    String temp_String_var; 
    for (auto& x : variable_mods_details_list)
    {
      StringList temp_split; 
      x.split(":", temp_split); 
      
      temp_split.insert(temp_split.begin()+1, ":["); 
      temp_split.insert(temp_split.end(), "]"); 
      String temp_split_Str = ""; 

      for (auto& y : temp_split)
      {
        temp_split_Str = temp_split_Str + y; 
      } 
      temp_String_var = temp_String_var + "," + temp_split_Str ; 
    } 
    String temp_String_var_Fin = temp_String_var.substr(1, temp_String_var.size()-1); 
    config_file.substitute("##static_mods##", static_mods_details);
    config_file.substitute("##variable_mods##", temp_String_var_Fin);

    return config_file;
  }

  std::tuple<std::string, std::string, std::string> getVersionNumber_(const std::string& multi_line_input)
  {
      std::regex version_regex("Version ([0-9]+)\\.([0-9]+)\\.([0-9]+)");

      std::sregex_iterator it(multi_line_input.begin(), multi_line_input.end(), version_regex);
      std::smatch match = *it;
      std::cout << "Found Sage version string: " << match.str() << std::endl;      
          
      return make_tuple(it->str(1), it->str(2), it->str(3)); // major, minor, patch
  }

  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<files>", StringList(), "Input files separated by blank");
    setValidFormats_("in", { "mzML" } );

    registerOutputFile_("out", "<file>", "", "Single output file containing all search results.", true, false);
    setValidFormats_("out", { "idXML" } );

    registerInputFile_("database", "<file>", "", "FASTA file", true, false, {"skipexists"});
    setValidFormats_("database", { "FASTA" } );

    registerInputFile_("sage_executable", "<executable>",
      // choose the default value according to the platform where it will be executed
      #ifdef OPENMS_WINDOWSPLATFORM
        "sage.exe",
      #else
        "sage",
      #endif
      "The Sage executable. Provide a full or relative path, or make sure it can be found in your PATH environment.", true, false, {"is_executable"});

    registerStringOption_("decoy_prefix", "<prefix>", "DECOY_", "Prefix on protein accession used to distinguish decoy from target proteins. NOTE: Decoy suffix is currently not supported by sage.", false, false);
    registerIntOption_("batch_size", "<int>", 0, "Number of files to load and search in parallel (default = # of CPUs/2)", false, false);
    
    registerDoubleOption_("precursor_tol_left", "<double>", -6.0, "Start (left side) of the precursor tolerance window w.r.t. precursor location. Usually used with negative values smaller or equal to the 'right' counterpart.", false, false);
    registerDoubleOption_("precursor_tol_right", "<double>", 6.0, "End (right side) of the precursor tolerance window w.r.t. precursor location. Usually used with positive values larger or equal to the 'left' counterpart.", false, false);
    registerStringOption_("precursor_tol_unit", "<unit>", "ppm", "Unit of precursor tolerance (ppm or Da)", false, false);
    setValidStrings_("precursor_tol_unit", ListUtils::create<String>("ppm,Da"));

    registerDoubleOption_("fragment_tol_left", "<double>", -20.0, "Start (left side) of the fragment tolerance window w.r.t. precursor location. Usually used with negative values smaller or equal to the 'right' counterpart.", false, false);
    registerDoubleOption_("fragment_tol_right", "<double>", 20.0, "End (right side) of the fragment tolerance window w.r.t. precursor location. Usually used with positive values larger or equal to the 'left' counterpart.", false, false);
    registerStringOption_("fragment_tol_unit", "<unit>", "ppm", "Unit of fragment tolerance (ppm or Da)", false, false);
    setValidStrings_("fragment_tol_unit", ListUtils::create<String>("ppm,Da"));

    // add advanced options
    registerIntOption_("min_matched_peaks", "<int>", min_matched_peaks, "Minimum number of b+y ions required to match for PSM to be reported", false, true);
    registerIntOption_("min_peaks", "<int>", min_peaks, "Minimum number of peaks required for a spectrum to be considered", false, true);
    registerIntOption_("max_peaks", "<int>", max_peaks, "Take the top N most intense MS2 peaks only for matching", false, true);
    registerIntOption_("report_psms", "<int>", report_psms, "Number of hits (PSMs) to report for each spectrum", false, true);
    registerIntOption_("bucket_size", "<int>", bucket_size, "How many fragments are in each internal mass bucket (default: 8192 for hi-res data). Try increasing it to 32k or 64k for low-res. See also: fragment_tol_*", false, true);
    registerIntOption_("min_len", "<int>", min_len, "Minimum peptide length", false, true);
    registerIntOption_("max_len", "<int>", max_len, "Maximum peptide length", false, true);
    registerIntOption_("missed_cleavages", "<int>", missed_cleavages, "Number of missed cleavages", false, true);
    registerDoubleOption_("fragment_min_mz", "<double>", fragment_min_mz, "Minimum fragment m/z", false, true);
    registerDoubleOption_("fragment_max_mz", "<double>", fragment_max_mz, "Maximum fragment m/z", false, true);
    registerDoubleOption_("peptide_min_mass", "<double>", peptide_min_mass, "Minimum monoisotopic peptide mass to consider a peptide from the DB", false, true);
    registerDoubleOption_("peptide_max_mass", "<double>", peptide_max_mass, "Maximum monoisotopic peptide mass to consider a peptide from the DB", false, true);
    registerIntOption_("min_ion_index", "<int>", min_ion_index, "Minimum ion index to consider for preliminary scoring. Default = 2 to skip b1/y1 AND (sic) b2/y2 ions that are often missing.", false, true);
    registerIntOption_("max_variable_mods", "<int>", max_variable_mods, "Maximum number of variable modifications", false, true);  
    registerStringOption_("isotope_error_range", "<start,end>", isotope_errors, "Range of (C13) isotope errors to consider for precursor."
      "Can be negative. E.g. '-1,3' for considering '-1/0/1/2/3'", false, true);
    registerStringOption_("charges", "<start,end>", charges_if_not_annotated, "Range of precursor charges to consider if not annotated in the file."
      , false, true);
    

    //Search Enzyme
    vector<String> all_enzymes;
    ProteaseDB::getInstance()->getAllNames(all_enzymes);
    registerStringOption_("enzyme", "<cleavage site>", "Trypsin", "The enzyme used for peptide digestion.", false, false);
    setValidStrings_("enzyme", all_enzymes);

    //Modifications
    vector<String> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
    registerStringList_("fixed_modifications", "<mods>", ListUtils::create<String>("Carbamidomethyl (C)", ','), "Fixed modifications, specified using Unimod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'", false);
    setValidStrings_("fixed_modifications", all_mods);
    registerStringList_("variable_modifications", "<mods>", ListUtils::create<String>("Oxidation (M)", ','), "Variable modifications, specified using Unimod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'", false);
    setValidStrings_("variable_modifications", all_mods);

    //FDR and misc 

    registerDoubleOption_("q_value_threshold", "<double>", 1, "The FDR threshhold for filtering peptides", false, false); 
    registerStringOption_("annotate_matches", "<bool>", "true", "If the matches should be annotated (default: false),", false, false); 
    registerStringOption_("deisotope", "<bool>", "false", "Sets deisotope option (true or false), default: false", false, false ); 
    registerStringOption_("chimera", "<bool>", "false", "Sets chimera option (true or false), default: false", false, false  ); 
    registerStringOption_("predict_rt",  "<bool>", "false", "Sets predict_rt option (true or false), default: false", false, false ); 
    registerStringOption_("wide_window", "<bool>", "false", "Sets wide_window option (true or false), default: false", false, false);
    registerStringOption_("smoothing", "<bool>", "true", "Should the PTM histogram be smoothed and local maxima be picked. If false, uses raw data, default: false", false, false);  
    registerIntOption_("threads", "<int>", 1, "Amount of threads available to the program", false, false); 

    // register peptide indexing parameter (with defaults for this search engine)
    registerPeptideIndexingParameter_(PeptideIndexing().getParameters());
  }


  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------

    // do this early, to see if Sage is installed
    String sage_executable = getStringOption_("sage_executable");
    std::cout << sage_executable << " sage executable" << std::endl; 
    String proc_stdout, proc_stderr;
    TOPPBase::ExitCodes exit_code = runExternalProcess_(sage_executable.toQString(), QStringList() << "--help", proc_stdout, proc_stderr, "");
    if (exit_code != EXECUTION_OK)
    {
      return exit_code;
    }

    auto major_minor_patch = getVersionNumber_(proc_stdout);
    String sage_version = std::get<0>(major_minor_patch) + "." + std::get<1>(major_minor_patch) + "." + std::get<2>(major_minor_patch);
    
    //-------------------------------------------------------------
    // run sage
    //-------------------------------------------------------------
    StringList input_files = getStringList_("in");
    String output_file = getStringOption_("out");
    String output_folder = File::path(output_file);
    String fasta_file = getStringOption_("database");
    int batch = getIntOption_("batch_size");
    String decoy_prefix = getStringOption_("decoy_prefix");

    // create config
    String config = imputeConfigIntoTemplate();

    // store config in config_file
    OPENMS_LOG_INFO << "Creating temp file name..." << std::endl;
    String config_file = File::getTempDirectory() + "/" + File::getUniqueName() + ".json";
    OPENMS_LOG_INFO << "Creating Sage config file..." << config_file << std::endl;
    ofstream config_stream(config_file.c_str());
    config_stream << config;
    config_stream.close();

    // keep config file if debug mode is set
    if (getIntOption_("debug") > 1)
    {
      String debug_config_file = output_folder + "/" + File::getUniqueName() + ".json";
      ofstream debug_config_stream(debug_config_file.c_str());
      debug_config_stream << config;
      debug_config_stream.close();     
    }

    String annotation_check;    

    QStringList arguments;

  if ( (getStringOption_("annotate_matches").compare("true")) == 0)
  {
    arguments << config_file.toQString() 
              << "-f" << fasta_file.toQString() 
              << "-o" << output_folder.toQString() 
              << "--annotate-matches"
              << "--write-pin"; 
  }
  else
  {
    arguments << config_file.toQString() 
              << "-f" << fasta_file.toQString() 
              << "-o" << output_folder.toQString() 
              << "--write-pin"; 
  }

    if (batch >= 1) arguments << "--batch-size" << String(batch).toQString();
    
    for (auto s : input_files) arguments << s.toQString();

    OPENMS_LOG_INFO << "Sage command line: " << sage_executable << " " << arguments.join(' ').toStdString() << std::endl;
    
    //std::chrono lines for testing/writing purposes only! 

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    // Sage execution with the executable and the arguments StringList
    exit_code = runExternalProcess_(sage_executable.toQString(), arguments);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    #ifdef CHRONOSET
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;
    #endif

    if (exit_code != EXECUTION_OK)
    {
      return exit_code;
    }

    //-------------------------------------------------------------
    // writing IdXML output
    //-------------------------------------------------------------

    // read the sage output
    OPENMS_LOG_INFO << "Reading sage output..." << std::endl;
    StringList filenames;
    StringList extra_scores = {"ln(-poisson)", "ln(delta_best)", "ln(delta_next)", 
      "ln(matched_intensity_pct)", "longest_b", "longest_y", 
      "longest_y_pct", "matched_peaks", "scored_candidates"}; 
    double FDR_threshhold = getDoubleOption_("q_value_threshold"); 

    vector<PeptideIdentification> peptide_identifications = PercolatorInfile::load(
      output_folder + "/results.sage.pin",
      true,
      "ln(hyperscore)", //TODO can we get sage's "sage_discriminant_score" out of the pin? Probably not. Suboptimal!
      extra_scores,
      filenames,
      decoy_prefix, 
      FDR_threshhold, 
      true);

    for (auto& id : peptide_identifications)
    {
      auto& hits = id.getHits();
      for (auto& h : hits)
      {
        for (const auto& meta : extra_scores)
        {
          if (h.metaValueExists(meta))
          {
            h.setMetaValue("SAGE:" + meta, h.getMetaValue(meta));
            h.removeMetaValue(meta);    
          }
        }
      }
    }
    
    String smoothing_string = getStringOption_("smoothing"); 
    bool smoothing = !(smoothing_string.compare("true")); 

    const  pair<DeltaMassHistogram, DeltaMasstoCharge> resultsClus =  getDeltaClusterCenter(peptide_identifications, smoothing, false); 
    vector<PeptideIdentification> mapD = mapDifftoMods(resultsClus.first, resultsClus.second, peptide_identifications, 0.01, false, output_file); //peptide_identifications; 
    // remove hits without charge state assigned or charge outside of default range (fix for downstream bugs). TODO: remove if all charges annotated in sage
    IDFilter::filterPeptidesByCharge(peptide_identifications, 2, numeric_limits<int>::max());
    
    if (filenames.empty()) filenames = getStringList_("in");

    // TODO: allow optional split and create multiple idXMLs one per input file
    vector<ProteinIdentification> protein_identifications(1, ProteinIdentification());

    writeDebug_("write idXMLFile", 1);    
    
    protein_identifications[0].setPrimaryMSRunPath(filenames);  
    protein_identifications[0].setDateTime(DateTime::now());
    protein_identifications[0].setSearchEngine("Sage");
    protein_identifications[0].setSearchEngineVersion(sage_version);

    DateTime now = DateTime::now();
    String identifier("Sage_" + now.get());
    protein_identifications[0].setIdentifier(identifier);
    for (auto & pid : peptide_identifications) 
    { 
      pid.setIdentifier(identifier);
      pid.setScoreType("ln(hyperscore)");
      pid.setHigherScoreBetter(true);
    }

    auto& search_parameters = protein_identifications[0].getSearchParameters();
    // protein_identifications[0].getSearchParameters().enzyme_term_specificity = static_cast<EnzymaticDigestion::Specificity>(num_enzyme_termini[getStringOption_("num_enzyme_termini")]);
    protein_identifications[0].getSearchParameters().db = getStringOption_("database");
    
    // add extra scores for percolator rescoring
    vector<String> percolator_features = { "score" };
    for (auto s : extra_scores) percolator_features.push_back("SAGE:" + s);
    search_parameters.setMetaValue("extra_features",  ListUtils::concatenate(percolator_features, ","));
    auto enzyme = *ProteaseDB::getInstance()->getEnzyme(getStringOption_("enzyme"));
    search_parameters.digestion_enzyme = enzyme; // needed for indexing
    search_parameters.enzyme_term_specificity = EnzymaticDigestion::SPEC_FULL;

    search_parameters.charges = "2:5"; // probably hard-coded in sage https://github.com/lazear/sage/blob/master/crates/sage/src/scoring.rs#L301

    search_parameters.mass_type = ProteinIdentification::MONOISOTOPIC;
    search_parameters.fixed_modifications = getStringList_("fixed_modifications");
    search_parameters.variable_modifications = getStringList_("variable_modifications");
    search_parameters.missed_cleavages = getIntOption_("missed_cleavages");
    search_parameters.fragment_mass_tolerance = (std::fabs(getDoubleOption_("fragment_tol_left")) + std::fabs(getDoubleOption_("fragment_tol_right"))) * 0.5;
    search_parameters.precursor_mass_tolerance = (std::fabs(getDoubleOption_("precursor_tol_left")) + std::fabs(getDoubleOption_("precursor_tol_right"))) * 0.5;
    search_parameters.precursor_mass_tolerance_ppm = getStringOption_("precursor_tol_unit") == "ppm";
    search_parameters.fragment_mass_tolerance_ppm = getStringOption_("fragment_tol_unit") == "ppm";

    // write all (!) parameters as metavalues to the search parameters
    if (!protein_identifications.empty())
    {
      DefaultParamHandler::writeParametersToMetaValues(this->getParam_(), protein_identifications[0].getSearchParameters(), this->getToolPrefix());
    }

    // if "reindex" parameter is set to true: will perform reindexing
    if (auto ret = reindex_(protein_identifications, peptide_identifications); ret != EXECUTION_OK) return ret;

    map<String,unordered_map<int,String>> file2specnr2nativeid;
    for (const auto& mzml : input_files)
    {
      // TODO stream mzml?
      MzMLFile m;
      MSExperiment exp;
      auto opts = m.getOptions();
      opts.setMSLevels({2,3});
      opts.setFillData(false);
      //opts.setMetadataOnly(true);
      m.setOptions(opts);
      m.load(mzml, exp);
      String nIDType = "";
      if (!exp.getSourceFiles().empty())
      {
        // TODO we could also guess the regex from the first nativeID if it is not stored here
        //  but I refuse to link to Boost::regex just for this
        //  Someone has to rework the API first!
        nIDType = exp.getSourceFiles()[0].getNativeIDTypeAccession();
      }

      for (const auto& spec : exp)
      {
        const String& nID = spec.getNativeID();
        int nr = SpectrumLookup::extractScanNumber(nID, nIDType);
        if (nr >= 0)
        {
          auto [it, inserted] = file2specnr2nativeid.emplace(File::basename(mzml), unordered_map<int,String>({{nr,nID}}));
          if (!inserted)
          {
            it->second.emplace(nr,nID);
          }
        }
      }
    }

    map<Size, String> idxToFile;
    StringList fnInRun;
    protein_identifications[0].getPrimaryMSRunPath(fnInRun);
    Size cnt = 0;
    for (const auto& f : fnInRun)
    {
      idxToFile.emplace(cnt, f);
      ++cnt;
    }

    for (auto& id : peptide_identifications)
    {
      Int64 scanNrAsInt = 0;
      
      try
      { // check if spectrum reference is a string that just contains a number        
        scanNrAsInt = id.getSpectrumReference().toInt64();
        // no exception -> conversion to int was successful. Now lookup full native ID in corresponding file for given spectrum number.
        id.setSpectrumReference( file2specnr2nativeid[idxToFile[id.getMetaValue(Constants::UserParam::ID_MERGE_INDEX)]].at(scanNrAsInt) );                              
      }
      catch (...)
      {
      }
    }
    IdXMLFile().store(output_file, protein_identifications, peptide_identifications);
    return EXECUTION_OK;
  }
};


int main(int argc, const char** argv)
{
  TOPPSageAdapter tool;
  return tool.main(argc, argv);
}

/// @endcond
