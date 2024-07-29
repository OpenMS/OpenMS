// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow, Philipp Wang $
// $Authors: Philipp Wang $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/ANALYSIS/ID/NeighborSeq.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/KERNEL/MSSpectrum.h>




using namespace OpenMS;
using namespace std;




// Function to generate the theoretical spectrum for a given peptide sequence
MSSpectrum NeighborSeq::generateSpectrum(const String& peptide_sequence)
{
AASequence peptide = AASequence::fromString(peptide_sequence);
TheoreticalSpectrumGenerator spec_gen;
// Set parameters for the spectrum generator
Param params;
params.setValue("add_b_ions", "true");
params.setValue("add_y_ions", "true");
spec_gen.setParameters(params);
// Generate the spectrum and store it in an MSSpectrum object
MSSpectrum spectrum;
spec_gen.getSpectrum(spectrum, peptide, 1, 1);
return spectrum;
}

// Function to compare two spectra and determine if they are similar
bool NeighborSeq::compareSpectra(const MSSpectrum& spec1, const MSSpectrum& spec2, const double& min_shared_ion_fraction, const double& mz_bin_size)
{

  // Calculate the number of shared bins considering the bin frequencies
  int B12 = compareShareSpectra(spec1,  spec2, mz_bin_size);

  // Calculate the fraction of shared bins
  double fraction_shared = static_cast<double>(2 * B12) / (spec1.size() + spec2.size());

  return fraction_shared > min_shared_ion_fraction;
}

//Finds candidate positions based on a given mono-isotopic weight and mass tolerance.
vector<int> NeighborSeq::findCandidatePositions(const map<double, vector<int>>& mass_position_map,const double& mono_weight, const double& mass_tolerance)
{
  // Vector to store the candidate positions
  vector<int> candidate_position;
  // Calculate the lower and upper bounds for the mass tolerance range
  double lower_bound_key = mono_weight;
  double upper_bound_key = mono_weight;
  if (mass_tolerance < 0) 
  {
     lower_bound_key += mass_tolerance;
     upper_bound_key -= mass_tolerance;
  }
  else
  {
    lower_bound_key -= mass_tolerance;
    upper_bound_key += mass_tolerance;
  }
 
  // Find the lower bound iterator in the map
  auto lower = mass_position_map.lower_bound(lower_bound_key);

  // Find the upper bound iterator in the map
  auto upper = mass_position_map.upper_bound(upper_bound_key);

  // Iterate through the map from the lower bound to the upper bound
  for (auto it = lower; it != upper; ++it)
  {
    // Insert the positions from the current map entry into the candidate positions vector
    candidate_position.insert(candidate_position.end(), it->second.begin(), it->second.end());
  }
  return candidate_position;
}


//Creates a map of masses to positions from a vector of peptides.(not protein)
map<double, vector<int>> NeighborSeq::createMassPositionMap(const vector<AASequence>& candidates)
{
  // Map to store the mass and corresponding positions
  map<double, vector<int>> mass_position_map;

  // Iterate through the vector of AASequence objects
  for (size_t i = 0; i < candidates.size(); ++i)
  {
    if (candidates[i].toString().find('X') == String::npos)
    {
      // Calculate the mono-isotopic mass of the sequence
      double mass = candidates[i].getMonoWeight();

      // Insert the mass and the position into the map
      mass_position_map[mass].push_back(i);
    }
   // else
   // {
          //throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Cannot get weight of sequence with unknown AA 'X' with unknown mass.", toString());
    //}
  }
  return mass_position_map;
}

// Method to find neighbor peptides in a given FASTA file
vector<int> NeighborSeq::findNeighborPeptides(const AASequence& peptide,
                                              const vector<AASequence>& neighbor_candidates,
                                              const vector<int>& candidates_position,
                                              const double& min_shared_ion_fraction,
                                              const double& mz_bin_size)

{
    vector<int> result_entries;
    MSSpectrum spec = generateSpectrum(peptide.toString());
    for (size_t i = 0; i < candidates_position.size(); i++)
    {
          
          // Check if the sequence contains an 'X'
          if (neighbor_candidates[candidates_position[i]].toString().find('X') == String::npos)
              {     
                  MSSpectrum neighbor_spec = generateSpectrum(neighbor_candidates[candidates_position[i]].toString());
               // Compare the spectra and add to results if they are similar
                  if (compareSpectra(spec, neighbor_spec, min_shared_ion_fraction, mz_bin_size))
                  {      
                      result_entries.push_back(candidates_position[i]);
                  }
              }
             // else
            //  {
                //throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Cannot get peaks of sequence with unknown AA 'X'.", toString());
            //  }
    }
return result_entries;
}




int NeighborSeq::compareShareSpectra(const MSSpectrum& spec1, const MSSpectrum& spec2, const double& mz_bin_size)
{
  std::map<int, int> bins1, bins2;

  // Lambda function to bin a spectrum into a map
  auto bin_spectrum = [&](const MSSpectrum& spec, std::map<int, int>& bins) {
    for (const auto& peak : spec)
      bins[static_cast<int>(peak.getMZ() / mz_bin_size)]++;
  };

  // Bin both spectra
  bin_spectrum(spec1, bins1);
  bin_spectrum(spec2, bins2);

  // Extract bin keys as vectors
  std::vector<int> vec_bins1, vec_bins2;
  for (const auto& bin : bins1)
    vec_bins1.push_back(bin.first);
  for (const auto& bin : bins2)
    vec_bins2.push_back(bin.first);

  // Calculate the intersection of the bin vectors
  std::vector<int> intersection;
  std::set_intersection(vec_bins1.begin(), vec_bins1.end(), vec_bins2.begin(), vec_bins2.end(), std::back_inserter(intersection));
  // Calculate the number of shared bins considering the bin frequencies
  int shared_peaks = 0;
  for (int bin : intersection)
    shared_peaks += min(bins1[bin], bins2[bin]);

  return shared_peaks;
}